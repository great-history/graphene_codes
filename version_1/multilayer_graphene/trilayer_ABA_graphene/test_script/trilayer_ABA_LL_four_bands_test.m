%% 该脚本计算的是三层ABA石墨烯四能带模型下的朗道能级图
% 该脚本分为两个部分
% 第一部分：Landau level as a function of E & B
% 第二部分：Landau level as a function of E & D @ fixed B

%% 添加路径
addpath(".\trilayer_ABA_hamiltonian\")
addpath(".\plot_funcs\")

%% parameters set up
%% 基本参数
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
one_electron_coulomb = 1.602176634 * 10^(-19); % 单位是C
epsilon_0 = 8.85 * 10^(-12); % 单位是F/m
epsilon_bn = 6.6;
d_interlayer = 0.335; % 单位是nm
a_intralayer = 0.246; % 单位是nm

%% 能带参数(以eV为单位)
% SWMC parameters from different references
% gamma0  gamma1  gamma2  gamma3  gamma4  gamma5  delta  Delta2
%   3.1    0.39   -0.028   0.315   0.041   0.05   0.046   0.0
%   3.1    0.38   -0.021   0.29    0.141   0.05   0.0355  0.0035   % Andrea Young PRL
%   3.16   0.39   -0.020   0.315   0.044   0.038  0.037   0.0      % Pablo Nature Physics
gamma0 = 3.16; 
gamma1 = 0.39;
gamma2 = -0.020;
gamma3 = 0.315;
gamma4 = 0.044;
gamma5 = 0.038;
delta = 0.037;
Delta2 = 0.0;

%% 注意点:
% 将gamma_i转换为v_i，同时某些gamma_i需要扩大sqrt(2)倍!!!!
gamma1 = gamma1 * sqrt(2);
v0 = sqrt(3) / 2 * gamma0 * (a_intralayer * 10^(-9)) / h_bar;
v3 = sqrt(3) / 2 * sqrt(2) * gamma3 * (a_intralayer * 10^(-9)) / h_bar; % 不要忘了乘以sqrt(2) !!!
v4 = sqrt(3) / 2 * sqrt(2) * gamma4 * (a_intralayer * 10^(-9)) / h_bar; % 不要忘了乘以sqrt(2) !!!

%% 外界参数
%  u = 0.2 * gamma1; % 0.0 / 0.04 / 0.065 / 0.07 / 0.08
% u的作用：high LL中电荷主要集中在第二层(u > 0)或第一层(u < 0)，而LLL中电荷主要集中在第二层(u > 0)或第一层(u <
% 0) ———— 由此也可推测出写代码时假定了u > 0时第一层电势比第二层低
% 当 | u | 比较小时，high LL在两层上的电荷分布几乎是差不多的，而lowest LL会存在层极化
% 当 | u | 比较大时，high LL 和 lowest LL都会存在层极化， LLL的层极化可能会被改变
% 当 | u | 非常非常大时：
Delta1 = 0.0; 

B_start = 0.0;
B_end = 3;
B_steps = 301;
B_fields_list = linspace(B_start, B_end, B_steps);

%% 第一部分：Landau level as a function of E & B
%% 朗道能级参数
LL_index_cutoff = 30;
dims_m = 2 * LL_index_cutoff + 1;
dims_b = 2 * LL_index_cutoff;
dims = dims_b + dims_m;
Ham_LL_K = zeros(dims);
Ham_LL_Kp = zeros(dims);
eigvals_LL_K = zeros(B_steps, dims);
eigvals_LL_Kp = zeros(B_steps, dims);

%% 计算LL fan diagram
tic
for B_index = 1:B_steps
    %% parameters set up
    B_field = B_fields_list(B_index);
    mag_length = 25.66 / sqrt(B_field); % 以nm为单位
    
    x0 = h_bar * v0 / mag_length * 10^9; % 以eV为单位
    x3 = h_bar * v3 / mag_length * 10^9; % 以eV为单位
    x4 = h_bar * v4 / mag_length * 10^9; % 以eV为单位
    
    %% construct Hamiltonian @ valley K
    [Ham_m_LL_K, Ham_b_LL_K] = construct_trilayer_ABA_LL_four_bands(x0, x3, gamma1, delta, gamma2, gamma5, Delta2, +1, LL_index_cutoff, dims_m, dims_b);
    [Ham_m_LL_Kp, Ham_b_LL_Kp] = construct_trilayer_ABA_LL_four_bands(x0, x3, gamma1, delta, gamma2, gamma5, Delta2, -1, LL_index_cutoff, dims_m, dims_b);
    
    % helper_check_hermite(Ham_LL_K, 1e-8);
    
    %% 对哈密顿量进行对角化
    % call the eig sovler
    % 分别对monolayer-like brach和bilayer-like brach进行求解（在这里，四能带模型不考虑加电位移场的情形）
    [eigvec_HK_b_now, eigval_HK_b] = eig(Ham_b_LL_K);
    eigval_HK_b_diag_now = diag(eigval_HK_b);
    [eigvec_HK_m_now, eigval_HK_m] = eig(Ham_m_LL_K);
    eigval_HK_m_diag_now = diag(eigval_HK_m);

    [eigvec_HKp_b_now, eigval_HKp_b] = eig(Ham_b_LL_Kp);
    eigval_HKp_b_diag_now = diag(eigval_HKp_b);
    [eigvec_HKp_m_now, eigval_HKp_m] = eig(Ham_m_LL_Kp);
    eigval_HKp_m_diag_now = diag(eigval_HKp_m);

    % push into the LLs
    eigvals_LL_K(B_index, 1:dims_m) = eigval_HK_m_diag_now;
    eigvals_LL_K(B_index, (dims_m + 1):end) = eigval_HK_b_diag_now;

    eigvals_LL_Kp(B_index, 1:dims_m) = eigval_HKp_m_diag_now;
    eigvals_LL_Kp(B_index, (dims_m + 1):end) = eigval_HKp_b_diag_now;
    
end
toc

%% 作图
eigvals_LL_cell = cell(2, 4);
eigvals_LL_cell{1,1} = eigvals_LL_K;
eigvals_LL_cell{1,2} = dims; % 维数
eigvals_LL_cell{1,3} = 'b--'; % line color
eigvals_LL_cell{1,4} = 0.5; % line width

eigvals_LL_cell{2,1} = eigvals_LL_Kp;
eigvals_LL_cell{2,2} = dims;
eigvals_LL_cell{2,3} = 'b';
eigvals_LL_cell{2,4} = 1.0;

save_path = "";
fig0 = plot_LLs(B_start, B_end, -0.100, 0.100, B_fields_list, eigvals_LL_cell, save_path);
% fig1 = plot_LLs_weights(eigvec_HK_now(:,2 * LL_index_cutoff - 3 :2 * LL_index_cutoff + 3)); % 对本征态在LL basis下的权重作图