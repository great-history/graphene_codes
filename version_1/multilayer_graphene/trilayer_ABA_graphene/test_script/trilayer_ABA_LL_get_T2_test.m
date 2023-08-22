%% 该脚本是对ABC-stacked trilayer graphene的两能带模型的LLL进行Hartree-Fock simulation
% ABC-stacked trilayer graphene的两能带模型与bilayer的两能带模型在数学结构上很相似，这得益于它们相同的手性特征
% 该脚本分为两部分：
% 第一部分：LL能带结构的计算，提取出LLL的本征态和本征能量信息
% 第二部分：考虑LLL之间的相互作用(包括Hartree & Fock)

%% 添加路径
addpath('.\interaction_effects\')
addpath('..\utils\')
addpath('..\hartree_fock_package\')

%% 参数准备
% 基本参数
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
one_electron_coulomb = 1.602176634 * 10^(-19); % 单位是C
epsilon_0 = 8.85 * 10^(-12) * 10^(-9); % 单位是F/nm
epsilon_bn = 6.6; % 无量纲，the in-plane dielectric constant of BN
d_interlayer = 0.335; % 单位是nm
a_intralayer = 0.246; % 单位是nm

B_field = 1.25;
mag_length = 25.6 / sqrt(B_field); % 以nm为单位

% energy scale
E_zeeman = 0.116 * B_field / 1000; % 单位是eV
E_F = one_electron_coulomb / (4 * pi * epsilon_0 * mag_length); % 单位是 eV fock interaction strength
E_exchange = E_F / epsilon_bn; % 由于hBN存在介电常数（这里是in-plane dielectric constant），这里取为6.6，它可以有效减小交换相互作用的强度
E_H = d_interlayer / (2 * mag_length) * E_F; % 单位是 eV hartree interaction strength, 注意在ABC trilayer那里是(2 * d / l_b)

%% 能带参数(以eV为单位)
% SWMC parameters from different references
% gamma0  gamma1  gamma2  gamma3  gamma4  gamma5  delta  Delta2
%   3.1    0.39   -0.028   0.315   0.041   0.05   0.046   0.0
%   3.1    0.38   -0.021   0.29    0.141   0.05   0.0355  0.0035   % Andrea Young PRL
%   3.16   0.39   -0.020   0.315   0.044   0.038  0.037   0.0      % Pablo Nature Physics
gamma0 = 3.1; 
gamma1 = 0.38;
gamma2 = -0.021;
gamma3 = 0.29;
gamma4 = 0.141;
gamma5 = 0.05;
delta = 0.0355;
Delta2 = 0.0035;

%% 注意点:
% 将gamma_i转换为v_i，同时某些gamma_i需要扩大sqrt(2)倍!!!!
gamma1 = gamma1 * sqrt(2);
v0 = sqrt(3) / 2 * gamma0 * (a_intralayer * 10^(-9)) / h_bar;
v3 = sqrt(3) / 2 * sqrt(2) * gamma3 * (a_intralayer * 10^(-9)) / h_bar; % 不要忘了乘以sqrt(2) !!!
v4 = sqrt(3) / 2 * sqrt(2) * gamma4 * (a_intralayer * 10^(-9)) / h_bar; % 不要忘了乘以sqrt(2) !!!

B_field = 1.25;
mag_length = 25.66 / sqrt(B_field); % 以nm为单位
x0 = h_bar * v0 / mag_length * 10^9; % 以eV为单位
x3 = h_bar * v3 / mag_length * 10^9; % 以eV为单位
x4 = h_bar * v4 / mag_length * 10^9; % 以eV为单位

Delta1_steps = 121;
Delta1_start = 0.0;
Delta1_end = 0.12;
Delta1_list = linspace(Delta1_start, Delta1_end, Delta1_steps); % 以eV为单位

%% 朗道能级参数
LL_index_cutoff = 30;
dims_m = 2 * LL_index_cutoff + 1;
dims_b = 4 * LL_index_cutoff;
dims = dims_b + dims_m;
Ham_LL_K = zeros(dims);
Ham_LL_Kp = zeros(dims);
eigvals_LL_K_Delta1 = zeros(Delta1_steps, dims);
eigvals_LL_Kp_Delta1 = zeros(Delta1_steps, dims);

% 用来存放T2 triplet states
T2_LL_index_list = [round(dims / 2) - 2, round(dims / 2) - 1, round(dims / 2)]; % 是人为找出来的
T2_eigvals = zeros(Delta1_steps, 3);
T2_eigvecs = zeros(Delta1_steps, dims, 3);
% T2_eigvals = eigvals_LL_K_Delta1(:, T2_LL_index_list);

%% 第一部分：LL能带结构的计算(LLs as E & D)，提取出LLL的本征态和本征能量信息
tic
for D_index = 1:Delta1_steps
    %% parameters set up
    Delta1 = Delta1_list(D_index);
    
    %% construct Hamiltonian @ valley K
    [Ham_m_LL_K, Ham_b_LL_K, D_K] = construct_trilayer_ABA_LL_six_bands(x0, x3, x4, gamma1, delta, gamma2, gamma5, Delta1, Delta2, +1, LL_index_cutoff, dims_m, dims_b);
    [Ham_m_LL_Kp, Ham_b_LL_Kp, D_Kp] = construct_trilayer_ABA_LL_six_bands(x0, x3, x4, gamma1, delta, gamma2, gamma5, Delta1, Delta2, -1, LL_index_cutoff, dims_m, dims_b);
    
    % helper_check_hermite(Ham_LL_K, 1e-8);
    
    %% 对哈密顿量进行对角化
    % call the eig sovler
    Ham_LL_K(1:dims_m, 1:dims_m) = Ham_m_LL_K;
    Ham_LL_K((dims_m + 1):end, (dims_m + 1):end) = Ham_b_LL_K;
    Ham_LL_K((dims_m + 1):end, 1:dims_m) = D_K;
    Ham_LL_K(1:dims_m, (dims_m + 1):end) = D_K';

    [eigvec_HK_now, eigval_HK] = eig(Ham_LL_K);
    eigval_HK_diag_now = diag(eigval_HK);

    Ham_LL_Kp(1:dims_m, 1:dims_m) = Ham_m_LL_Kp;
    Ham_LL_Kp((dims_m + 1):end, (dims_m + 1):end) = Ham_b_LL_Kp;
    Ham_LL_Kp((dims_m + 1):end, 1:dims_m) = D_Kp;
    Ham_LL_Kp(1:dims_m, (dims_m + 1):end) = D_Kp';

    [eigvec_HKp_now, eigval_HKp] = eig(Ham_LL_Kp);
    eigval_HKp_diag_now = diag(eigval_HKp);

    % push into the LLs
    eigvals_LL_K_Delta1(D_index, :) = eigval_HK_diag_now;
    eigvals_LL_Kp_Delta1(D_index, :) = eigval_HKp_diag_now;
    
    % 提取出T2 triplets states的本征态和本征能量
    T2_eigvals(D_index, :) = eigval_HK_diag_now(T2_LL_index_list);
    T2_eigvecs(D_index, :, :) = eigvec_HK_now(:,T2_LL_index_list);
end
toc

%% 对T2 triplet进行作图
% Extended view of the triplet T2 with the average energy of the triplet substracted
ene_av_list = sum(T2_eigvals, 2) / 3;
fig4 = figure();
axis([Delta1_start Delta1_end -0.0025 0.0015])
hold on
for i = 1:3
    plot(Delta1_list, eigvals_LL_K_Delta1(:, T2_LL_index_list(i)) - ene_av_list, 'r')
end

%% 第二部分：考虑T2 triplet states之间的相互作用(包括Hartree & Fock)
% 只计算Delta1 = 0.08eV 时的对称性破缺

% 找到Delta1 = 0.08eV附近的几个值和相应的index
for D_index = 1:Delta1_steps
    if abs(Delta1_list(D_index) - 0.08) <= 0.0005
        D_select_index = D_index;
    end
end

T2_eigvecs_cell = cell(LL_index_cutoff + 1, 1); % 用来存放T2 states在不同LL_index下的块
T2_eigvals_list = T2_eigvals(D_select_index, :);
for n = 0:LL_index_cutoff
    T2_eigvecs_cell{n + 1} = zeros(6, 3); % 第一个指标是phi_k, 第二个指标是T2 state index
end

for ii = 1:3
    eigvec = T2_eigvecs(D_select_index, :, ii);
    eigvec = squeeze(eigvec);
    eigvec_phi_k_component_array = get_LL_components_each_phi_k(eigvec, LL_index_cutoff, +1, dims_m, dims_b);
    for n = 0:LL_index_cutoff
        T2_eigvecs_cell{n + 1}(:, ii) = eigvec_phi_k_component_array(n + 1, :);
    end
end

% 确定所需考虑的最高LL_index（我们不需要考虑到LL_index_cutoff，因为一计算量大，二没有必要）
threshold = 0.90; % 只有权重大于threshold的LL_index才会被记录
LL_index_max = 0;
for ii = 1:3
    norm_sum = 0.0;
    for n = 0:LL_index_cutoff
        eigvec_weights = T2_eigvecs_cell{n + 1}(:, ii);
        norm_sum = norm_sum + sum(abs(eigvec_weights).^2);
        if sqrt(norm_sum) >= threshold
            if n > LL_index_max
                LL_index_max = n;
            end
            
            break
        end
    end
end

% 变换矩阵
T_mat = zeros(6); % 从{|alpha>}变换到{|phi_k>}的矩阵
T_mat(1,1) = 1 / sqrt(2);
T_mat(1,3) = 1 / sqrt(2);
T_mat(2,2) = 1 / sqrt(2);
T_mat(2,4) = 1 / sqrt(2);
T_mat(3,5) = 1;
T_mat(4,6) = 1;
T_mat(5,1) = - 1 / sqrt(2);
T_mat(5,3) = 1 / sqrt(2);
T_mat(6,2) = - 1 / sqrt(2);
T_mat(6,4) = 1 / sqrt(2);

%% 保存数据(T2 本征能量和本征态)
file_path = ['.\trilayer_ABA_data\T2_states_info_B_', num2str(roundn(B_field, -2)), 'T_LL_index_max_', num2str(LL_index_max), '.mat'];
save(file_path, 'LL_index_max', 'T2_eigvecs_cell', 'T2_eigvals_list', 'T_mat');

% 计算这几个Delta1下T2 states的自发性对称性破缺
