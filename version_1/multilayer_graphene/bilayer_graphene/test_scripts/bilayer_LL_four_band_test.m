%% 该脚本计算的是双层石墨烯四能带模型下的朗道能级图
% 该脚本分为两个部分
% 第一部分：Landau level as a function of E & B
% 第二部分：estimation of gap by self-consistent calculation

%% parameters set up
v = 1 * 10^6; % 来自gamma0，单位是m/s
gamma1 = 0.40;
% gamma3 = 0.1 * gamma0; % 0.0 / 0.1 * gamma0
v3 = 0.1 * v; % 0.0 * v / 0.1 * v
v4 = 0;
delta = 0;
%  u = 0.2 * gamma1; % 0.0 / 0.04 / 0.065 / 0.07 / 0.08
% u的作用：high LL中电荷主要集中在第二层(u > 0)或第一层(u < 0)，而LLL中电荷主要集中在第二层(u > 0)或第一层(u <
% 0) ———— 由此也可推测出写代码时假定了u > 0时第一层电势比第二层低
% 当 | u | 比较小时，high LL在两层上的电荷分布几乎是差不多的，而lowest LL会存在层极化
% 当 | u | 比较大时，high LL 和 lowest LL都会存在层极化， LLL的层极化可能会被改变
% 当 | u | 非常非常大时：
u = 0.04; 
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位

B_start = 0.0;
B_end = 6;
B_steps = 601;
B_fields_list = linspace(B_start, B_end, B_steps);

LL_index_cutoff = 30;
dims = 4 * LL_index_cutoff;
eigvals_LL_K = zeros(B_steps, dims);
eigvals_LL_Kp = zeros(B_steps, dims);

%% 计算LL fan diagram
tic
for B_index = 1:B_steps
    %% parameters set up
    B_field = B_fields_list(B_index);
    mag_length = 25.66 / sqrt(B_field); % 以nm为单位
    x0 = h_bar * v / mag_length * 10^9;
    x3 = h_bar * v3 / mag_length * 10^9;
    x4 = h_bar * v4 / mag_length * 10^9;
    
    %% construct Hamiltonian @ valley K
    Ham_LL_K = construct_bilayer_LL_four_bands(x0, x3, x4, gamma1, delta, u, +1, LL_index_cutoff, dims);
    Ham_LL_Kp = construct_bilayer_LL_four_bands(x0, x3, x4, gamma1, delta, u, -1, LL_index_cutoff, dims);
    
    % helper_check_hermite(Ham_LL_K, 1e-8);
    
    %% 对哈密顿量进行对角化
    % call the eig sovler
    [eigvec_HK_now, eigval_HK] = eig(Ham_LL_K);
    eigval_HK_diag_now = diag(eigval_HK);
    
    [eigvec_HKp_now, eigval_HKp] = eig(Ham_LL_Kp);
    eigval_HKp_diag_now = diag(eigval_HKp);

    % push into the LLs
    eigvals_LL_K(B_index, :) = eigval_HK_diag_now;
    eigvals_LL_Kp(B_index, :) = eigval_HKp_diag_now;
    
    [rho_layer1_K, rho_layer2_K] = get_layer_density_four_bands(eigvec_HK_now(:, 1:2 * LL_index_cutoff - 3), 1, LL_index_cutoff);
    [rho_layer1_Kp, rho_layer2_Kp] = get_layer_density_two_bands(eigvec_HKp_now(:, 1:2 * LL_index_cutoff - 3), -1, LL_index_cutoff);
    rho_layer1 = rho_layer1_K + rho_layer1_Kp;
    rho_layer2 = rho_layer2_K + rho_layer2_Kp;
    disp("第一层的电荷密度")
    rho_layer1
    disp("第二层的电荷密度")
    rho_layer2
    
%     [rho_layer1_K, rho_layer2_K] = get_layer_density_four_bands(eigvec_HK_now(:, 2 * LL_index_cutoff - 3:2 * LL_index_cutoff + 3), 1, LL_index_cutoff);
%     [rho_layer1_Kp, rho_layer2_Kp] = get_layer_density_four_bands(eigvec_HKp_now(:, 2 * LL_index_cutoff - 3:2 * LL_index_cutoff + 3), -1, LL_index_cutoff);
%     rho_layer1 = rho_layer1_K + rho_layer1_Kp;
%     rho_layer2 = rho_layer2_K + rho_layer2_Kp;
%     disp("LLL在第一层上的电荷密度")
%     rho_layer1
%     disp("LLL在第二层上的电荷密度")
%     rho_layer2
end
toc

%% 作图
addpath("..\plot_funcs\")
eigvals_LL_cell = cell(2, 4);
eigvals_LL_cell{1,1} = eigvals_LL_K / gamma1;
eigvals_LL_cell{1,2} = dims; % 维数
eigvals_LL_cell{1,3} = 'b'; % line color
eigvals_LL_cell{1,4} = 0.5; % line width

eigvals_LL_cell{2,1} = eigvals_LL_Kp / gamma1;
eigvals_LL_cell{2,2} = dims;
eigvals_LL_cell{2,3} = 'b--';
eigvals_LL_cell{2,4} = 1.0;

save_path = "";
fig0 = plot_LLs(B_start, B_end, -0.125, 0.125, B_fields_list, eigvals_LL_cell, save_path);
fig1 = plot_LLs_weights(eigvec_HK_now(:,2 * LL_index_cutoff - 3 :2 * LL_index_cutoff + 3));

%% estimation of u by self-consistent calculation
d_interlayer = 0.335; % 单位是nm
one_electron_coulomb = 1.602176634 * 10^(-19); % 单位是C
epsilon_0 = 8.85 * 10^(-12); % 单位是F/m
mag_length = 25.66 / sqrt(B_field); % 以nm为单位
u_factor = one_electron_coulomb * d_interlayer / (4 * pi * epsilon_0 * 10^(-9) * mag_length^2); 


u_external = u; % 外部的电位移场导致的电势能
[rho_layer1_K, rho_layer2_K] = get_layer_density_four_bands(eigvec_HK_now(:, 1:2 * LL_index_cutoff), 1, LL_index_cutoff);
[rho_layer1_Kp, rho_layer2_Kp] = get_layer_density_four_bands(eigvec_HKp_now(:, 1:2 * LL_index_cutoff), -1, LL_index_cutoff);
rho_layer1_new = rho_layer1_K + rho_layer1_Kp;
rho_layer2_new = rho_layer2_K + rho_layer2_Kp;

diff_rho = (rho_layer2_new - rho_layer1_new) * 2; % 即层间电容效应产生的内部电场，来部分抵消电位移场的影响 % 2来自自旋自由度
u_new = u_external - u_factor * diff_rho;
% u = u_new;