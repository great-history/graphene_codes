%% 添加路径
% addpath('../interaction_effects/')
% addpath('..')
addpath('../interaction_effects/')
addpath('../bilayer_hamitonians/')
addpath('../plot_funcs/')

B_field = 6.0;
d_interlayer = 0.335; % 单位是 nm
mag_length = 25.66 / sqrt(B_field); % 单位是 nm
% 基本参数
one_electron_coulomb = 1.602176634 * 10^(-19); % C
epsilon_0 = 8.85 * 10^(-12); % F/m
dielectric_const = 6;
epsilon_now = epsilon_0 * dielectric_const;

C_bg = 6.85 * 10^(14); % 单位是m^(-2)*V^(-1)
% 转换为nm^(-2)*V^(-1)
C_bg = C_bg * 10^(-18);
V_bg = 1; % 单位是V
electric_external = C_bg * V_bg * one_electron_coulomb; % 即电位移场
% electric_internal = 0.0; % 即层间电容效应产生的内部电场，来部分抵消电位移场的影响 % 2来自自旋自由度

% (下面这种表示是错误的，因为单位是焦耳)
% u = one_electron_coulomb * d_interlayer / (epsilon_0 * 10^(-9)); 
% u = u * (electric_external + electric_internal);
% (下面这种表示是正确的，因为单位是eV)
u_factor = d_interlayer / (2 * epsilon_0 * 10^(-9)); 
% u = u_factor * (electric_external - electric_internal);
% conventions : 第二层的电势为 - u / 2 ; 第一层的电势为 + u / 2
% u_external_list = linspace(0.2, 1.2, 11) * d_interlayer;
u_external_list = linspace(0.0, 0.12, 13);

%% 设置参数：
v = 1 * 10^6; % 来自gamma0，单位是m/s
gamma1 = 0.38;
% gamma3 = 0.1 * gamma0; % 0.0 / 0.1 * gamma0
v3 = 0.1 * v; % 0.0 * v / 0.1 * v
v4 = 0;
delta = 0;
% u = 0.2 * gamma1; % 0.0 / 0.04 / 0.065 / 0.07 / 0.08
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
single_electron = 1.6021892 * 10^(-19); % 以C为单位

B_field = 6.0; % 1.0 / 2.2 / 3.0 / 3.5
mag_length = 25.66 / sqrt(B_field); % 以nm为单位
x0 = h_bar * v / mag_length * 10^9;
x3 = h_bar * v3 / mag_length * 10^9;
x4 = h_bar * v4 / mag_length * 10^9;

LL_index_cutoff = 100;
dims = 2 * LL_index_cutoff;
eigvals_LL_K = zeros(1, dims);
eigvals_LL_Kp = zeros(1, dims);

%% 计算加入电位移场后的LL能带结构
u_self_consistent_list = zeros(1,11);
for u_index = 1:11
    u_external = u_external_list(u_index);
    u = u_external;
    
    %% construct Hamiltonian @ valley K
    Ham_LL_K = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, +1, LL_index_cutoff, dims);
    Ham_LL_Kp = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, -1, LL_index_cutoff, dims);

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

    %% 计算每一层上的电荷密度(单位是？？？ nm^(-2)) 
    % [rho_layer1_K, rho_layer2_K] = get_layer_density_two_bands(eigvec_HK_now(:, 1:negative_index_list(end)), 1, LL_index_cutoff);
    % [rho_layer1_Kp, rho_layer2_Kp] = get_layer_density_two_bands(eigvec_HKp_now(:, 1:(positive_index_list(1)-1)), -1, LL_index_cutoff);
    [rho_layer1_K, rho_layer2_K] = get_layer_density_two_bands(eigvec_HK_now(:, 1:LL_index_cutoff), 1, LL_index_cutoff);
    [rho_layer1_Kp, rho_layer2_Kp] = get_layer_density_two_bands(eigvec_HKp_now(:, 1:LL_index_cutoff - 1), -1, LL_index_cutoff);
    rho_layer1_new = rho_layer1_K + rho_layer1_Kp;
    rho_layer2_new = rho_layer2_K + rho_layer2_Kp;

    % ele_density_diff = 2 * one_electron_coulomb / (2 * pi * mag_length^2) * (rho_layer1_new - rho_layer2_new); % 单位是 C * nm^(-2)
    % electric_internal = (rho_layer1_new - rho_layer2_new) * one_electron_coulomb * 2 / (4 * pi * mag_length^2); % 即层间电容效应产生的内部电场，来部分抵消电位移场的影响 % 2来自自旋自由度
    electric_internal = (rho_layer2_new - rho_layer1_new) * one_electron_coulomb * 2 / (2 * pi * mag_length^2); % 即层间电容效应产生的内部电场，来部分抵消电位移场的影响 % 2来自自旋自由度
    % u_new = u_factor * (electric_external - electric_internal);
    u_new = u_external - u_factor * electric_internal;
    u = u_new;
    
    %% 进行自洽计算
    for step = 1:10
        %% construct Hamiltonian @ valley K
        Ham_LL_K = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, +1, LL_index_cutoff, dims);
        Ham_LL_Kp = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, -1, LL_index_cutoff, dims);

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

        %% 计算每一层上的电荷密度(单位是？？？ nm^(-2))
        rho_layer1_last = rho_layer1_new;
        rho_layer2_last = rho_layer2_new;
        % [rho_layer1_K, rho_layer2_K] = get_layer_density_two_bands(eigvec_HK_now(:, 1:negative_index_list(end)), 1, LL_index_cutoff);
        % [rho_layer1_Kp, rho_layer2_Kp] = get_layer_density_two_bands(eigvec_HKp_now(:, 1:(positive_index_list(1)-1)), -1, LL_index_cutoff);
        [rho_layer1_K, rho_layer2_K] = get_layer_density_two_bands(eigvec_HK_now(:, 1:LL_index_cutoff), 1, LL_index_cutoff);
        [rho_layer1_Kp, rho_layer2_Kp] = get_layer_density_two_bands(eigvec_HKp_now(:, 1:LL_index_cutoff), -1, LL_index_cutoff);
        rho_layer1_new = rho_layer1_K + rho_layer1_Kp;
        rho_layer2_new = rho_layer2_K + rho_layer2_Kp;

        % disp("层密度之差为：")
        % abs(rho_layer1_last - (rho_layer1_K + rho_layer1_Kp))
        % abs(rho_layer2_last - (rho_layer2_K + rho_layer2_Kp))

        % ele_density_diff = 2 * one_electron_coulomb / (2 * pi * mag_length^2) * (rho_layer1_new - rho_layer2_new); % 单位是 C * nm^(-2)
        electric_internal = (rho_layer1_new - rho_layer2_new) * one_electron_coulomb * 2 / (2 * pi * mag_length^2); % 即层间电容效应产生的内部电场，来部分抵消电位移场的影响 % 2来自自旋自由度
        % electric_internal = (rho_layer2_new - rho_layer1_new) * one_electron_coulomb * 2 / (2 * pi * mag_length^2); % 即层间电容效应产生的内部电场，来部分抵消电位移场的影响 % 2来自自旋自由度
        % u_new = u_factor * (electric_external - electric_internal);
        u_new = u_external + u_factor * electric_internal;
        % disp("u之差为：")
        % abs(u - u_new)
        u = u_new; 
    end
    
    disp("u之差为：")
    abs(u - u_new)
    
    u_self_consistent_list(u_index) = u_new;
end