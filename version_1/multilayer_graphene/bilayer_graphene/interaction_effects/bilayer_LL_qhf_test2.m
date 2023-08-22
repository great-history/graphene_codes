%% test
% addpath('.\interaction_effects\')
format long
%% 构造S矩阵
S0000 = 0.64;
S0011 = 0.47;
S0212 = 0.29;
S0102 = -0.29;

S_tensor = zeros(3,3,3,3);
for alpha = 1:3
    for beta = 1:3
        for lambda = 1:3
            for sigma = 1:3
                % 非零必须至少有两个指标一样
                if (alpha == beta) && (alpha == lambda) && (alpha == sigma)
                    S_tensor(alpha, beta, lambda, sigma) = S0000;
                    continue % 一定要加continue!!!
                end
                
                if (alpha == beta) && (lambda == sigma)
                    S_tensor(alpha, beta, lambda, sigma) = S0011;
                    continue % 一定要加continue!!!
                end
                
                if (alpha == sigma) && (lambda == beta)
                    S_tensor(alpha, beta, lambda, sigma) = S0011;
                    continue % 一定要加continue!!!
                end
                
                if (beta == sigma)
                    if ~(alpha == lambda) && ~(alpha == beta) && ~(lambda == beta)
                        if beta == 3
                            S_tensor(alpha, beta, lambda, sigma) = S0212;
                        else
                            S_tensor(alpha, beta, lambda, sigma) = S0102;
                        end
                    end
                    
                    continue
                end
                
                if (alpha == lambda)
                    if ~(beta == sigma) && ~(beta == alpha) && ~(sigma == alpha)
                        if alpha == 3
                            S_tensor(alpha, beta, lambda, sigma) = S0212;
                        else
                            S_tensor(alpha, beta, lambda, sigma) = S0102;
                        end
                    end
                    
                    continue
                end
            end
        end
    end
end

%% energy scale
B_field = 1;
% E_zeeman = 0.289 * B_field / 2; % 不要忘记除以2
E_zeeman = 0.116 * B_field;
E_ee = 140 * sqrt(B_field) / 6;  % 6为相对介电常数
filling_factor = 1;

% 构造一个随机矩阵
density_matrix_temp = construct_random_density_matrix(filling_factor, 3, "complex");
% density_matrix_temp = construct_random_density_matrix(filling_factor, 3, "real");  % 如果是实数的话，很容易进入local minima，并且没有复数那样严谨

% eigvec1 = [-0.2600; 0.8011; 0.5391];
% eigvec2 = [0.7734; -0.1614; 0.6130];
% density_matrix_temp = eigvec1 * eigvec1' + eigvec2 * eigvec2';

% 只能是哈密顿量的构造存在问题
H_hf = construct_H_hf_spin(density_matrix_temp, S_tensor, E_ee, 0);
H_hf = refine_H_hf(H_hf, 3);
[error, density_matrix_temp] = Hartree_Fock_iteration_test(H_hf, density_matrix_temp, filling_factor);
% % density_matrix_last = density_matrix_temp;
% disp(["误差为", error])

%% start self-consistent calculation
% tic
steps = 0;
while error >= 1e-8
    H_hf = construct_H_hf_spin(density_matrix_temp, S_tensor, E_ee, 0);
    H_hf = refine_H_hf(H_hf, 3);
    [error, density_matrix_temp] = Hartree_Fock_iteration_test(H_hf, density_matrix_temp, filling_factor);
    steps = steps + 1;
    % disp(["误差为", error])
end

if steps == 0
    steps = 10;
end

for ii = 1:steps*100
    H_hf = construct_H_hf_spin(density_matrix_temp, S_tensor, E_ee, 0);
    H_hf = refine_H_hf(H_hf, 3);
    [error, density_matrix_temp] = Hartree_Fock_iteration_test(H_hf, density_matrix_temp, filling_factor);
    if ~mod(ii,100)
        disp(["误差为", error])
        % [eigvecs_temp, eigvals_temp] = eig(H_hf);
        % disp(["本征值为：", eigvals_temp(1,1), eigvals_temp(2,2), eigvals_temp(3,3)])
        disp("密度矩阵为：")
        density_matrix_temp
    end
end

[eigvecs_temp, eigvals_temp] = eig(H_hf);
[~, index_list] = sort(real(diag(eigvals_temp)), 'ascend');

% disp(["误差为", error])
% [eigvecs_temp, eigvals_temp] = eig(H_hf);
% toc

%% 设置参数：
v = 1 * 10^6; % 来自gamma0，单位是m/s
gamma1 = 0.40;
% gamma3 = 0.1 * gamma0; % 0.0 / 0.1 * gamma0
v3 = 0.1 * v; % 0.0 * v / 0.1 * v
v4 = 0;
delta = 0;
u = 0.2 * gamma1; % 0.0 / 0.04 / 0.065 / 0.07 / 0.08
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
single_electron = 1.6021892 * 10^(-19); % 以C为单位

B_field = 1.0; % 1.0 / 2.2 / 3.0 / 3.5
mag_length = 25.66 / sqrt(B_field); % 以nm为单位
x0 = h_bar * v / mag_length * 10^9;
x3 = h_bar * v3 / mag_length * 10^9;
x4 = h_bar * v4 / mag_length * 10^9;

LL_index_cutoff = 12;
dims = 2 * LL_index_cutoff;
eigvals_LL_K = zeros(1, dims);
eigvals_LL_Kp = zeros(1, dims);

%% construct Hamiltonian @ valley K
Ham_LL_K = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, +1, LL_index_cutoff, dims);
Ham_LL_Kp = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, -1, LL_index_cutoff, dims);

% helper_check_hermite(Ham_LL_K, 1e-8);

%% 对哈密顿量进行对角化
% call the eig sovler
[eigvec_HK_now, eigval_HK] = eig(Ham_LL_K);
eigval_HK_diag_now = diag(eigval_HK);

[eigvec_HKp_now, eigval_HKp] = eig(Ham_LL_Kp);
eigval_HKp_diag_now = diag(eigval_HKp);

% push into the LLs
eigvals_LL_K(:) = eigval_HK_diag_now;
eigvals_LL_Kp(:) = eigval_HKp_diag_now;

%% 选取最靠近CNP的6条朗道能级(3条位于导带，3条位于价带) ———— triplets states
% K valley 取价带能量绝对值最小的三个, Kp valley 取导带能量绝对值最小的三个
negative_index_list = find(eigvals_LL_K < 0);
negative_ene_vals = eigvals_LL_K(eigvals_LL_K < 0);
[~, sort_index_list] = sort(negative_ene_vals, 'descend');
eigvals_LLL_K = eigvals_LL_K(negative_index_list(sort_index_list(1:3)));
eigvecs_LLL_K = eigvec_HK_now(:, negative_index_list(sort_index_list(1:3)));

positive_index_list = find(eigvals_LL_Kp > 0);
positive_ene_vals = eigvals_LL_Kp(eigvals_LL_Kp > 0);
[~, sort_index_list] = sort(positive_ene_vals, 'ascend');
eigvals_LLL_Kp = eigvals_LL_Kp(positive_index_list(sort_index_list(1:3)));
eigvecs_LLL_Kp = eigvec_HKp_now(:, positive_index_list(sort_index_list(1:3)));

% 这里只考虑Kp的三个LL能级
eigvecs_LLL_Kp_A_component = zeros(LL_index_cutoff + 1, 3);
eigvecs_LLL_Kp_B_component = zeros(LL_index_cutoff + 1, 3);
for ii = 1:3
    [eigvec_A_component, eigvec_B_component] = get_LL_components_each_sublattice_two_bands(eigvecs_LLL_Kp(:, ii), LL_index_cutoff, -1);
    eigvecs_LLL_Kp_A_component(:, ii) = eigvec_A_component;
    eigvecs_LLL_Kp_B_component(:, ii) = eigvec_B_component;
end

eigvecs_LLL_Kp_triplet = zeros(dims, filling_factor);
eigvecs_LLL_Kp_A_component_triplet = zeros(LL_index_cutoff + 1, filling_factor);
eigvecs_LLL_Kp_B_component_triplet = zeros(LL_index_cutoff + 1, filling_factor);
for lambda = 1:filling_factor
    index = index_list(lambda);
    S = eigvecs_temp(:, index);
    
    eigvecs_LLL_Kp_triplet(:, lambda) = eigvecs_LLL_Kp(:, 1) * S(1) + eigvecs_LLL_Kp(:, 2) * S(2) + eigvecs_LLL_Kp(:, 3) * S(3);
    eigvecs_LLL_Kp_A_component_triplet(:, lambda) = eigvecs_LLL_Kp_A_component(:, 1) * S(1) + eigvecs_LLL_Kp_A_component(:, 2) * S(2) + eigvecs_LLL_Kp_A_component(:, 3) * S(3);
    eigvecs_LLL_Kp_B_component_triplet(:, lambda) = eigvecs_LLL_Kp_B_component(:, 1) * S(1) + eigvecs_LLL_Kp_B_component(:, 2) * S(2) + eigvecs_LLL_Kp_B_component(:, 3) * S(3);
end

%% 计算每个(kx,ky)处的权重:|\phi(kx,ky)|^2
% k mesh grid
num_points = 301;
kx_list = gamma1 / (h_bar * v * 10^9) * linspace(-0.3, 0.3, num_points);
ky_list = gamma1 / (h_bar * v * 10^9) * linspace(-0.3, 0.3, num_points);
[kx_mesh, ky_mesh] = meshgrid(kx_list, ky_list);
phi_k_mesh_triplet = zeros(filling_factor, num_points, num_points); % triplet * kx * ky

for kx_index = 1:num_points
    kx = kx_list(kx_index);
    for ky_index = 1:num_points
        ky = ky_list(ky_index);
        alpha_eig = 1j * mag_length * (kx + 1j * ky) / sqrt(2);
        
        for lambda = 1:filling_factor
            T_coeff = 0.0;
            norm = 0.0;
            for n = 0:LL_index_cutoff
                T_coeff = T_coeff + eigvecs_LLL_Kp_A_component_triplet(n + 1, lambda) * (alpha_eig)^n / sqrt((factorial(n)));
                norm = norm + (abs(eigvecs_LLL_Kp_A_component_triplet(n + 1, lambda)))^2;
            end
            T_coeff = (abs(T_coeff))^2;
            T_coeff = T_coeff * exp(-(abs(alpha_eig))^2);
            T_coeff = T_coeff / norm;

            phi_k_mesh_triplet(lambda, ky_index, kx_index) = T_coeff;
        end
        
    end
end

addpath('.\plot_funcs\')
% rescale
save_path = ['D:\matlab\graphene-package\multilayer_graphene\bilayer_graphene\bilayer_figs\bilayer_LLLs_qhf_k_dist_B_', ...
                num2str(roundn(B_field, -1)), '_v3_', num2str(roundn(v3 / v, -2)), '_filling_', num2str(filling_factor), '.jpg'];
fig1 = plot_LLs_k_dist(phi_k_mesh_triplet, kx_list * (h_bar * v * 10^9) / gamma1, ky_list * (h_bar * v * 10^9) / gamma1, filling_factor, save_path);