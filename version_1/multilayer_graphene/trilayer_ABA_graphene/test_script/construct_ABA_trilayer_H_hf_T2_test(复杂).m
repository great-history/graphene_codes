%%
% 这个脚本用了比较复杂的构造哈密顿量的方法，不太好
% 添加路径

%% 加载数据
B_field = 1.25;
LL_index_max = 13;
% 加载T2 states相关的信息
file_path = ['.\trilayer_ABA_data\T2_states_info_B_', num2str(roundn(B_field, -2)), 'T_LL_index_max_', num2str(LL_index_max), '.mat'];
load(file_path);
% 加载交换积分
file_path = ['.\trilayer_ABA_data\exchange_integrals_B_', num2str(roundn(B_field, -2)), 'T_LL_index_max_', num2str(LL_index_max), '.mat'];
load(file_path);

%% 基本参数
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
one_electron_coulomb = 1.602176634 * 10^(-19); % 单位是C
epsilon_0 = 8.85 * 10^(-12) * 10^(-9); % 单位是F/nm
epsilon_bn = 6.6; % 无量纲，the in-plane dielectric constant of BN
d_interlayer = 0.335; % 单位是nm
a_intralayer = 0.246; % 单位是nm

%% energy scale
num_LL = LL_index_max + 1;
mag_length = 25.6 / sqrt(B_field); % 以nm为单位
E_zeeman = 0.116 * B_field / 1000; % 单位是eV
E_F = one_electron_coulomb / (4 * pi * epsilon_0 * mag_length); % 单位是 eV fock interaction strength
E_exchange = E_F / epsilon_bn; % 由于hBN存在介电常数（这里是in-plane dielectric constant），这里取为6.6，它可以有效减小交换相互作用的强度
E_H = d_interlayer / (2 * mag_length) * E_F; % 单位是 eV hartree interaction strength, 注意在ABC trilayer那里是(2 * d / l_b)

%% 其它矩阵
% 将本征态从phi_k basis转换到alpha basis
% 在eigvec_phi_k_component_Delta1_array中已经将本征态按照
% [|0,phi1>,|1,phi1>,|2,phi1>,|3,phi1>,...] // [|0,phi2>,|1,phi2>,|2,phi2>,|3,phi2>,...] // 
% [|0,phi3>,|1,phi3>,|2,phi3>,|3,phi3>,...] // [|0,phi4>,|1,phi4>,|2,phi4>,|3,phi4>,...] //
% [|0,phi5>,|1,phi5>,|2,phi5>,|3,phi5>,...] // [|0,phi6>,|1,phi6>,|2,phi6>,|3,phi6>,...] //

% transform_mat = zeros(num_LL, 6, 3); % 第二个指标是子格指标，第三个指标是本征态指标
transform_mat_cell = cell(num_LL, 1);
for n = 0:LL_index_max
    transform_mat_cell{n + 1} = zeros(6, 3);
    
    eig_vec_mat = zeros(6,3);
    for ii = 1:3
        eig_vec = T2_eigvecs_cell{n + 1}(:, ii);
        % eig_vec = reshape(eig_vec, [6,1]);
        eig_vec_mat(:, ii) = eig_vec;
    end
    
    transform_mat_cell{n + 1} = T_mat * eig_vec_mat;
end

%% 随机构造一个密度矩阵(in the basis of triplet T2 states)
dim_subspace = 3; % ABA trilayer graphene 的 T2 states只有三个，并且都在K valley处
filling_factor = 1; % (整数)填充数

% 构造方法四（最科学） ： 在构造方法一的基础上添加一个权重(这样的话容易出现局域极小值的解)，权重随机生成
filling_list = randi([0 dim_subspace], 2, 1);
while filling_list == [0;0]
    filling_list = randi([0 dim_subspace], 2, 1);
end

density_matrix_down_temp = construct_random_density_matrix(filling_list(1), dim_subspace, "complex");
% density_matrix_up_temp = construct_random_density_matrix(filling_list(2), dim_subspace, "complex");
density_matrix_up_temp = zeros(dim_subspace);

weight1 = unifrnd(0,1);
weight_list = [weight1; 1 - weight1];
density_matrix_down_temp = weight_list(1) * density_matrix_down_temp;
density_matrix_up_temp = weight_list(2) * density_matrix_up_temp;
disp("自旋向上密度矩阵为：")
density_matrix_up_temp
disp("自旋向下密度矩阵为：")
density_matrix_down_temp

% 计算第二层上的电荷密度Delta_mid
Delta_mid = 0.0;
for a = 1:3
    for b = 1:3
        
        val = 0.0;
        for n = 0:LL_index_max
            val = val + conj(eigvecs_phi_k_component_Delta1(a, n+1, 5)) * eigvecs_phi_k_component_Delta1(b, n+1, 5);
            val = val + conj(eigvecs_phi_k_component_Delta1(a, n+1, 6)) * eigvecs_phi_k_component_Delta1(b, n+1, 6);
        end
        Delta_mid = Delta_mid + (density_matrix_down_temp(a, b) + density_matrix_up_temp(a, b)) * val;
    end
end

% 第二种方式
% Delta_mid = 0.0;
% for n = 0:LL_index_max
%     index_start = 6 * n;
%     Delta_mid = Delta_mid + density_matrix_down_temp_LL(index_start + 5, index_start + 5);
%     Delta_mid = Delta_mid + density_matrix_down_temp_LL(index_start + 6, index_start + 6);
%     Delta_mid = Delta_mid + density_matrix_up_temp_LL(index_start + 5, index_start + 5);
%     Delta_mid = Delta_mid + density_matrix_up_temp_LL(index_start + 6, index_start + 6);
% end

% 将密度矩阵从eigenstate basis转换到LL basis
density_matrix_down_temp_LL = zeros(6*num_LL); % 只有一个谷,并且是spin down，所以是6*num_LL 【|A1,n>，|B1,n>，|A2,n>，|B2,n>，|A3,n>，|B3,n>】
density_matrix_up_temp_LL = zeros(6*num_LL); % 只有一个谷,并且是spin down，所以是6*num_LL 【|A1,n>，|B1,n>，|A2,n>，|B2,n>，|A3,n>，|B3,n>
for n = 0:LL_index_max
    index_left = 6 * n;
    transform_mat_left = conj(transform_mat_cell{n + 1});
    
    for m = 0:LL_index_max
        index_right = 6 * m;
        
        transform_mat_right = transform_mat_cell{m + 1};
        transform_mat_right = transpose(transform_mat_right);
        
        density_matrix_down_temp_LL((index_left + 1):(index_left + 6), (index_right + 1):(index_right + 6)) = ...
                                                        transform_mat_left * density_matrix_down_temp * transform_mat_right;
        density_matrix_up_temp_LL((index_left + 1):(index_left + 6), (index_right + 1):(index_right + 6)) = ...
                                                        transform_mat_left * density_matrix_up_temp * transform_mat_right;
    end
end
density_matrix_down_temp_LL = refine_H_hf(density_matrix_down_temp_LL, 6*num_LL);
density_matrix_up_temp_LL = refine_H_hf(density_matrix_up_temp_LL, 6*num_LL);

% 构造Hartree-Fock Hamiltonian
% eigenvalues
H_hf_down = diag(eigvals_list);
H_hf_up = diag(eigvals_list);

% Zeeman term
H_hf_down = H_hf_down - E_zeeman * diag([1,1,1]);
H_hf_up = H_hf_up + E_zeeman * diag([1,1,1]);

% Hartree term
layer_potential_mat = diag([1,1,-1,-1,1,1]);
for n = 0:LL_index_max
    transform_mat_left = (transform_mat_cell{n + 1})';
    transform_mat_right = transform_mat_cell{n + 1};

    mat_temp = E_H / 2 * Delta_mid * (transform_mat_left * layer_potential_mat * transform_mat_right);
    
    H_hf_down(:, :) = H_hf_down(:, :) + mat_temp;
    H_hf_up(:, :) = H_hf_up(:, :) + mat_temp;
end

% Fock term
% 分别用1，2，3代表intralayer, nearest interlayer, next-nearest interlayer
dist_index_array = zeros(6,6);
for alpha = 1:6
    layer_alpha = mod(alpha, 2) + fix(alpha / 2);
    for beta = 1:6
        layer_beta = mod(beta, 2) + fix(beta / 2);
        if abs(layer_beta - layer_alpha) == 1
            dist_index_array(alpha, beta) = 2;
        elseif abs(layer_beta - layer_alpha) == 2
            dist_index_array(alpha, beta) = 3;
        else % 对应 alpha == beta
            dist_index_array(alpha, beta) = 1;
        end
    end
end

% 开始计算矩阵元<a|U_ex|b>
for n1 = 0:LL_index_max
    for m2 = 0:LL_index_max
        transform_mat_left = (transform_mat_cell{n1 + 1})';
        transform_mat_right = transform_mat_cell{m2 + 1};
        
        X_mat_up = zeros(6,6);
        X_mat_down = zeros(6,6);
        for m1 = 0:LL_index_max
            n2 = n1 + m1 - m2;
            if n2 < 0 || n2 > LL_index_max
                continue
            end
            
            index_left = 6 * m1;
            index_right = 6 * n2;
                
            a = abs(n1 - n2);
            b_prime = min(n1, n2);
            c_prime = min(m1, m2);
            b = max(b_prime, c_prime);
            c = min(b_prime, c_prime);
            ex_index = indice_cell{a + 1}{b + 1}{c + 1};
            
            for alpha = 1:6
                for beta = 1:6
                    dist_index = dist_index_array(alpha, beta);
                    X_mat_down(alpha, beta) = X_mat_down(alpha, beta) ...
                                            + exchange_integrals_list(ex_index, dist_index) * density_matrix_down_temp_LL(index_left + beta, index_right + alpha);
                    X_mat_up(alpha, beta) = X_mat_up(alpha, beta) ...
                                          + exchange_integrals_list(ex_index, dist_index) * density_matrix_up_temp_LL(index_left + beta, index_right + alpha);
                end
            end
            
        end
        X_mat_up = - E_F * X_mat_up;
        X_mat_down = - E_F * X_mat_down;
        
        H_hf_up = H_hf_up + transform_mat_left * X_mat_up * transform_mat_right;
        H_hf_down = H_hf_down + transform_mat_left * X_mat_down * transform_mat_right;
    end
end

helper_check_hermite(H_hf_down, 1e-12)
%% Hartree-Fock 初始化
% % 构造the full Hartree-Fock Hamitonian
% [H_hf_up, H_hf_down] = construct_ABC_trilayer_H_hf(density_matrix_up_temp, density_matrix_down_temp, ...
%                                                    exchange_integrals_intralayer_cell, exchange_integrals_interlayer_cell, E_exchange, E_H, E_zeeman);
% if ~helper_check_hermite(H_hf_down, 1e-3)
%     disp("H_hf_down不是厄米矩阵")
% end
% 
% H_hf_up = refine_H_hf(H_hf_up, dim_subspace);
% H_hf_down = refine_H_hf(H_hf_down, dim_subspace);
% [error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, filling_factor, dim_subspace);
% [eigvecs_down_last, eigvals_down_last] = eig(H_hf_down);
% [eigvecs_up_last, eigvals_up_last] = eig(H_hf_up);
% eigvals_down_last = diag(eigvals_down_last);
% eigvals_up_last = diag(eigvals_up_last);
% disp(["误差为", error])
