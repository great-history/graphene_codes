%% parameters set up
format long
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

addpath('.\interaction_effects\')
% 这里只考虑Kp的三个LL能级
eigvecs_LLL_Kp_A_component = zeros(LL_index_cutoff + 1, 3);
eigvecs_LLL_Kp_B_component = zeros(LL_index_cutoff + 1, 3);
for ii = 1:3
    [eigvec_A_component, eigvec_B_component] = get_LL_components_each_sublattice_two_bands(eigvecs_LLL_Kp(:, ii), LL_index_cutoff, -1);
    eigvecs_LLL_Kp_A_component(:, ii) = eigvec_A_component;
    eigvecs_LLL_Kp_B_component(:, ii) = eigvec_B_component;
end

%% k mesh grid
num_points = 301;
kx_list = gamma1 / (h_bar * v * 10^9) * linspace(-0.3, 0.3, num_points);
ky_list = gamma1 / (h_bar * v * 10^9) * linspace(-0.3, 0.3, num_points);
[kx_mesh, ky_mesh] = meshgrid(kx_list, ky_list);

% 计算每个(kx,ky)处的权重
phi_k_mesh = zeros(3, num_points, num_points); % triplet * kx * ky
for kx_index = 1:num_points
    kx = kx_list(kx_index);
    for ky_index = 1:num_points
        ky = ky_list(ky_index);
        alpha_eig = 1j * mag_length * (kx + 1j * ky) / sqrt(2);
        
        for lambda = 1:3
            T_coeff = 0.0;
            norm = 0.0;
            for n = 0:LL_index_cutoff
                % T_coeff = T_coeff + (eigvecs_LLL_Kp_A_component(n + 1, lambda) + eigvecs_LLL_Kp_B_component(n + 1, lambda)) * (alpha_eig)^n / sqrt((factorial(n)));
                T_coeff = T_coeff + eigvecs_LLL_Kp_A_component(n + 1, lambda) * (alpha_eig)^n / sqrt((factorial(n)));
                norm = norm + (abs(eigvecs_LLL_Kp_A_component(n + 1, lambda)))^2;
                % T_coeff = T_coeff + eigvecs_LLL_Kp_B_component(n + 1, lambda) * (alpha_eig)^n / sqrt((factorial(n)));
                % norm = norm + (abs(eigvecs_LLL_Kp_B_component(n + 1, lambda)))^2;
            end
            T_coeff = (abs(T_coeff))^2;
            T_coeff = T_coeff * exp(-(abs(alpha_eig))^2);
            T_coeff = T_coeff / norm;

            phi_k_mesh(lambda, ky_index, kx_index) = T_coeff;
        end
        
    end
end

%% 寻找三个valley @ pi, pi / 3, - pi / 3
% 划分区域(按照角度2*pi/3, 4*pi/3， 2*pi), (2*pi/3, 4*pi/3)是区域1，(4*pi/3,4*pi/3)是区域1，(2*pi/3, 4*pi/3)是区域3
region_k_mesh = zeros(num_points, num_points); 
for kx_index = 1:num_points
    kx = kx_list(kx_index);
    for ky_index = 1:num_points
        ky = ky_list(ky_index);
        
        % region 1 criterion
        if ky > sqrt(3) * kx && ky <= - sqrt(3) * kx
            region_k_mesh(ky_index, kx_index) = 1;
            continue
        end
        
        % region 2 criterion
        if ky < 0 && ky <= sqrt(3) * kx
            region_k_mesh(ky_index, kx_index) = 2;
            continue
        end
        
        % region 3 criterion
        if ky >= 0 && ky > - sqrt(3) * kx
            region_k_mesh(ky_index, kx_index) = 3;
            continue
        end
    end
end

k_valley_list = zeros(3, 2);
for lambda = 1:3
    region_index_list = find(region_k_mesh == lambda);
    phi_k_mesh_current = reshape(phi_k_mesh(lambda, :, :), num_points, num_points);
    [~, valley_index] = max(phi_k_mesh_current(region_index_list));
    valley_index = region_index_list(valley_index);
    
    kx_valley = kx_mesh(valley_index);
    ky_valley = ky_mesh(valley_index);
    k_valley_list(lambda, 1) = kx_valley;
    k_valley_list(lambda, 2) = ky_valley;
end
clear region_k_mesh

%% build a state from the above three degenerate states that is solely concentrated around a single valley
T_coeff_list = zeros(3, 3); % k_valley * triplet
for lambda = 1:3
    for valley = 1:3
        kx_valley = k_valley_list(valley, 1);
        ky_valley = k_valley_list(valley, 2);
        alpha_eig = 1j * mag_length * (kx_valley + 1j * ky_valley) / sqrt(2);
        T_coeff = 0.0;
        for n = 0:LL_index_cutoff
            % T_coeff = T_coeff + (eigvecs_LLL_Kp_A_component(n + 1, lambda) + eigvecs_LLL_Kp_B_component(n + 1, lambda)) * (alpha_eig)^n / sqrt((factorial(n)));
            T_coeff = T_coeff + eigvecs_LLL_Kp_A_component(n + 1, lambda) * (alpha_eig)^n / sqrt((factorial(n)));
            % T_coeff = T_coeff + eigvecs_LLL_Kp_B_component(n + 1, lambda) * (alpha_eig)^n / sqrt((factorial(n)));
        end
        T_coeff_list(valley, lambda) = T_coeff;
    end
end

eigvecs_LLL_Kp_triplet = zeros(dims, 3);
eigvecs_LLL_Kp_A_component_triplet = zeros(LL_index_cutoff + 1, 3);
eigvecs_LLL_Kp_B_component_triplet = zeros(LL_index_cutoff + 1, 3);
for lambda = 1:3
    B = zeros(3,1);
    B(lambda) = 1;
    S = linsolve(T_coeff_list, B);
    % 归一化
    norm_factor = sqrt(sum((abs(S)).^2));
    S = S ./ norm_factor;
    
    eigvecs_LLL_Kp_triplet(:, lambda) = eigvecs_LLL_Kp(:, 1) * S(1) + eigvecs_LLL_Kp(:, 2) * S(2) + eigvecs_LLL_Kp(:, 3) * S(3);
    eigvecs_LLL_Kp_A_component_triplet(:, lambda) = eigvecs_LLL_Kp_A_component(:, 1) * S(1) + eigvecs_LLL_Kp_A_component(:, 2) * S(2) + eigvecs_LLL_Kp_A_component(:, 3) * S(3);
    eigvecs_LLL_Kp_B_component_triplet(:, lambda) = eigvecs_LLL_Kp_B_component(:, 1) * S(1) + eigvecs_LLL_Kp_B_component(:, 2) * S(2) + eigvecs_LLL_Kp_B_component(:, 3) * S(3);
end

% 计算每个(kx,ky)处的权重
phi_k_mesh_triplet = zeros(3, num_points, num_points); % triplet * kx * ky
for kx_index = 1:num_points
    kx = kx_list(kx_index);
    for ky_index = 1:num_points
        ky = ky_list(ky_index);
        alpha_eig = 1j * mag_length * (kx + 1j * ky) / sqrt(2);
        
        for lambda = 1:3
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

eigvals_LLL_Kp_quasi = sum(eigvals_LLL_Kp) / 3 * ones(3,1);  % quasi-degenerate 近似简并

% 存放LLL信息
LLL_Kp_info_cell = cell(1,8);
LLL_Kp_info_cell{1, 1} = eigvals_LLL_Kp;
LLL_Kp_info_cell{1, 2} = eigvecs_LLL_Kp;
LLL_Kp_info_cell{1, 3} = eigvecs_LLL_Kp_A_component;
LLL_Kp_info_cell{1, 4} = eigvecs_LLL_Kp_B_component;
LLL_Kp_info_cell{1, 5} = eigvals_LLL_Kp_quasi;
LLL_Kp_info_cell{1, 6} = eigvecs_LLL_Kp_triplet;
LLL_Kp_info_cell{1, 7} = eigvecs_LLL_Kp_A_component_triplet;
LLL_Kp_info_cell{1, 8} = eigvecs_LLL_Kp_B_component_triplet;

%% 作图并保存数据
save_path = ['.\bilayer_data\LLL_Kp_info_B_', num2str(roundn(B_field, -1)), 'T.mat'];
save(save_path,'LLL_Kp_info_cell');

addpath('.\plot_funcs\')
% rescale
save_path = ['D:\matlab\graphene-package\multilayer_graphene\bilayer_graphene\bilayer_figs\bilayer_LLLs_k_dist_B_', num2str(roundn(B_field, -1)), '_v3_', num2str(roundn(v3 / v, -2)), '.jpg'];
fig1 = plot_LLs_k_dist(phi_k_mesh, kx_list * (h_bar * v * 10^9) / gamma1, ky_list * (h_bar * v * 10^9) / gamma1, 3, save_path);
save_path = ['D:\matlab\graphene-package\multilayer_graphene\bilayer_graphene\bilayer_figs\bilayer_LLLs_triplet_k_dist_B_', num2str(roundn(B_field, -1)), '_v3_', num2str(roundn(v3 / v, -2)), '.jpg'];
fig2 = plot_LLs_k_dist(phi_k_mesh_triplet, kx_list * (h_bar * v * 10^9) / gamma1, ky_list * (h_bar * v * 10^9) / gamma1, 3, save_path);