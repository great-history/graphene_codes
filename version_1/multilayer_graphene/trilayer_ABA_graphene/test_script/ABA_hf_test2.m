%% 添加路径
addpath('.\interaction_effects\');
addpath('..\hartree_fock_package\');
addpath('..\utils\')

%% 加载数据
B_field = 1.25;
LL_index_max = 12;
% 加载T2 states相关的信息
file_path = ['.\trilayer_ABA_data\T2_states_info_B_', num2str(roundn(B_field, -2)), 'T_LL_index_max_', num2str(LL_index_max), '.mat'];
load(file_path);
% 加载交换积分
file_path = ['.\trilayer_ABA_data\exchange_integrals_B_', num2str(roundn(B_field, -2)), 'T_LL_index_max_', num2str(LL_index_max), '.mat'];
load(file_path);

%% 基本参数
LL_index_max = 12;
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

%% 其它矩阵
% 将本征态从phi_k basis转换到alpha basis
% 在eigvec_phi_k_component_Delta1_array中已经将本征态按照
% [|0,phi1>,|1,phi1>,|2,phi1>,|3,phi1>,...] // [|0,phi2>,|1,phi2>,|2,phi2>,|3,phi2>,...] // 
% [|0,phi3>,|1,phi3>,|2,phi3>,|3,phi3>,...] // [|0,phi4>,|1,phi4>,|2,phi4>,|3,phi4>,...] //
% [|0,phi5>,|1,phi5>,|2,phi5>,|3,phi5>,...] // [|0,phi6>,|1,phi6>,|2,phi6>,|3,phi6>,...] //

% transform_mat = zeros(num_LL, 6, 3); % 第二个指标是子格指标，第三个指标是本征态指标
transform_mat_cell = cell(num_LL, 1);
for n = 0:LL_index_max
    transform_mat_cell{n + 1} = T_mat * T2_eigvecs_cell{n + 1};
end

%% 随机构造一个密度矩阵(in the basis of triplet T2 states)
dim_subspace = 3; % ABA trilayer graphene 的 T2 states只有三个，并且都在K valley处
filling_factor = 1; % (整数)填充数

% density_matrix_down_temp = diag([1,0,0]);
% density_matrix_up_temp = zeros(3);

% 构造方法三
density_matrix_down_temp = construct_random_density_matrix(1, dim_subspace, "complex");
% density_matrix_up_temp = construct_random_density_matrix(filling_list(2), dim_subspace, "complex");
density_matrix_up_temp = zeros(dim_subspace);

% % 构造方法四（最科学） ： 在构造方法一的基础上添加一个权重(这样的话容易出现局域极小值的解)，权重随机生成
% filling_list = randi([0 dim_subspace], 2, 1);
% while filling_list == [0;0]
%     filling_list = randi([0 dim_subspace], 2, 1);
% end
% 
% density_matrix_down_temp = construct_random_density_matrix(filling_list(1), dim_subspace, "complex");
% % density_matrix_up_temp = construct_random_density_matrix(filling_list(2), dim_subspace, "complex");
% density_matrix_up_temp = zeros(dim_subspace);
% 
% weight1 = unifrnd(0,1);
% weight_list = [weight1; 1 - weight1];
% density_matrix_down_temp = weight_list(1) * density_matrix_down_temp;
% density_matrix_up_temp = weight_list(2) * density_matrix_up_temp;


disp("自旋向上密度矩阵为：")
density_matrix_up_temp
disp("自旋向下密度矩阵为：")
density_matrix_down_temp

%% 将密度矩阵从eigenstate basis转换到LL basis
[density_matrix_up_temp_LL, density_matrix_down_temp_LL] = convert_density_matrix(density_matrix_up_temp, density_matrix_down_temp, ...
                                                                                  transform_mat_cell, num_LL, LL_index_max);

% density_matrix_down_temp_LL = refine_H_hf(density_matrix_down_temp_LL, 6*num_LL); % 保证严格的厄密性
% density_matrix_up_temp_LL = refine_H_hf(density_matrix_up_temp_LL, 6*num_LL); % 保证严格的厄密性

%% 计算第二层上的电荷密度Delta_mid
% % 第一种方式
% Delta_mid = get_delta_mid_method2(density_matrix_down_temp, density_matrix_up_temp, T2_eigvecs_cell, LL_index_max);
% 第二种方式 % 与第一种方法得到的结果应该是一致的，推荐使用第二种，算起来会快一些
Delta_mid = get_delta_mid_method1(density_matrix_up_temp_LL, density_matrix_down_temp_LL, LL_index_max);
% 
%% 构造Hartree-Fock Hamiltonian
[H_hf_up, H_hf_down] = construct_ABA_trilayer_H_hf_T2(density_matrix_up_temp_LL, density_matrix_down_temp_LL, exchange_integrals_list, ...
                                                      dist_index_array, indice_cell, transform_mat_cell, ...
                                                      T2_eigvals_list, E_zeeman, E_H, Delta_mid, E_exchange, LL_index_max);
H_hf_down = refine_H_hf(H_hf_down, 3);
H_hf_up = refine_H_hf(H_hf_up, 3);

[error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, ...
                                                                                   filling_factor, dim_subspace);
[eigvecs_down_last, eigvals_down_last] = eig(H_hf_down);
[eigvecs_up_last, eigvals_up_last] = eig(H_hf_up);
eigvals_down_last = diag(eigvals_down_last);
eigvals_up_last = diag(eigvals_up_last);
disp(["误差为", error])

%% self-consistent loop
steps = 0;
for ii = 1:1000
    [density_matrix_up_temp_LL, density_matrix_down_temp_LL] = convert_density_matrix(density_matrix_up_temp, density_matrix_down_temp, ...
                                                                                           transform_mat_cell, num_LL, LL_index_max);
    Delta_mid = get_delta_mid_method1( density_matrix_up_temp_LL, density_matrix_down_temp_LL, LL_index_max);
    [H_hf_up, H_hf_down] = construct_ABA_trilayer_H_hf_T2(density_matrix_up_temp_LL, density_matrix_down_temp_LL, exchange_integrals_list, ...
                                                          dist_index_array, indice_cell, transform_mat_cell, ...
                                                          T2_eigvals_list, E_zeeman, E_H, Delta_mid, E_exchange, LL_index_max);
    if ~helper_check_hermite(H_hf_down, 1e-16)
        disp("H_hf_down没有厄密性")
    end
                                                  
    H_hf_down = refine_H_hf(H_hf_down, 3);
    H_hf_up = refine_H_hf(H_hf_up, 3);
    
    [error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, ...
                                                                                       filling_factor, dim_subspace);
    steps = steps + 1;
    if ~mod(steps, 100)
        fprintf("进行到第%d步，误差为%10.9f\n", steps, error)
    end
    %     [eigvecs_down_temp, eigvals_down_temp] = eig(H_hf_down);
    %     [eigvecs_up_temp, eigvals_up_temp] = eig(H_hf_up);
    %     eigvals_down_temp = diag(eigvals_down_temp);
    %     eigvals_up_temp = diag(eigvals_up_temp);
    %
    %     error_down = max(abs(eigvals_down_temp - eigvals_down_last));
    %     error_up = max(abs(eigvals_up_temp - eigvals_up_last));
    %     % error = max(error_down, error_up);
    %     eigvals_down_last = eigvals_down_temp;
    %     eigvals_up_last = eigvals_up_temp;
end

density_matrix_down_list = zeros(dim_subspace, dim_subspace, 10); % 存放连续的10个密度矩阵
eigvals_down_list = zeros(dim_subspace, 10);
for ii = 1:10
    [density_matrix_up_temp_LL, density_matrix_down_temp_LL] = convert_density_matrix(density_matrix_up_temp, density_matrix_down_temp, ...
                                                                                           transform_mat_cell, num_LL, LL_index_max);
    Delta_mid = get_delta_mid_method1( density_matrix_up_temp_LL, density_matrix_down_temp_LL, LL_index_max);
    [H_hf_up, H_hf_down] = construct_ABA_trilayer_H_hf_T2(density_matrix_up_temp_LL, density_matrix_down_temp_LL, exchange_integrals_list, ...
                                                          dist_index_array, indice_cell, transform_mat_cell, ...
                                                          T2_eigvals_list, E_zeeman, E_H, Delta_mid, E_exchange, LL_index_max);
    % helper_check_hermite(H_hf_down, 1e-16)
                                                  
    H_hf_down = refine_H_hf(H_hf_down, 3);
    H_hf_up = refine_H_hf(H_hf_up, 3);
                                       
    [error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, ...
                                                                                       filling_factor, dim_subspace);
                                                                                   
    disp(["误差为", error])
    % disp("自旋向上密度矩阵为：")
    % density_matrix_up_temp
    % disp("自旋向下密度矩阵为：")
    % density_matrix_down_temp
    
    [eigvecs_down_temp, eigvals_down_temp] = eig(H_hf_down);
    density_matrix_down_list(:, :, ii) = density_matrix_down_temp;
    eigvals_down_list(:, ii) = real(diag(eigvals_down_temp));
    % [eigvecs_temp, eigvals_temp] = eig(H_hf_down);
    % disp(["本征值为：", eigvals_temp(1,1), eigvals_temp(2,2), eigvals_temp(3,3)])
end