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
%   3.1    0.38   -0.0257  0.288   0.050   0.056  0.0475  0.0     % 我们的0GPa的拟合结果
% gamma0 = 3.1; 
% gamma1 = 0.38;
% gamma2 = -0.0257;
% gamma3 = 0.288;
% gamma4 = 0.050;
% gamma5 = 0.056;
% delta = 0.0475;
% Delta2 = 0.0;

gamma0 = 3.1; 
gamma1 = 0.3892;
gamma2 = -0.0342;
gamma3 = 0.1569;  % 通过测试发现gamma3在[-0.29, 0.49]的范围内对LL_K_b_0和LL_Kp_b_0对应的本征态和本征能量的影响是可以降低的
gamma4 = 0.0440;
gamma5 = 0.0612;
delta = 0.0496;
Delta2 = 0.0020;

%% 注意点:
% 将gamma_i转换为v_i，同时某些gamma_i需要扩大sqrt(2)倍!!!!
gamma1 = gamma1 * sqrt(2);
v0 = sqrt(3) / 2 * gamma0 * (a_intralayer * 10^(-9)) / h_bar;
v3 = sqrt(3) / 2 * sqrt(2) * gamma3 * (a_intralayer * 10^(-9)) / h_bar; % 不要忘了乘以sqrt(2) !!!
v4 = sqrt(3) / 2 * sqrt(2) * gamma4 * (a_intralayer * 10^(-9)) / h_bar; % 不要忘了乘以sqrt(2) !!!

%% 第二部分：Landau level as a function of E & D @ fixed B
% 朗道能级参数
LL_index_cutoff = 30;
dims_m = 2 * LL_index_cutoff + 1;
dims_b = 4 * LL_index_cutoff;
dims = dims_b + dims_m;
Ham_LL_K = zeros(dims);
Ham_LL_Kp = zeros(dims);

B_field = 3.0;
mag_length = 25.66 / sqrt(B_field); % 以nm为单位
x0 = h_bar * v0 / mag_length * 10^9; % 以eV为单位
x3 = h_bar * v3 / mag_length * 10^9; % 以eV为单位
x4 = h_bar * v4 / mag_length * 10^9; % 以eV为单位

Delta1_steps = 121;
Delta1_start = 0.0;
Delta1_end = 0.4 * d_interlayer;
Delta1_list = linspace(Delta1_start, Delta1_end, Delta1_steps); % 以eV为单位
eigvals_LL_K_Delta1 = zeros(Delta1_steps, dims);
eigvals_LL_Kp_Delta1 = zeros(Delta1_steps, dims);

%% 计算LL fan diagram as E & D
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
    
end
toc

%% 作图
eigvals_LL_cell_Delta1 = cell(2, 4);
eigvals_LL_cell_Delta1{1,1} = eigvals_LL_K_Delta1;
eigvals_LL_cell_Delta1{1,2} = dims; % 维数
eigvals_LL_cell_Delta1{1,3} = 'b'; % line color
eigvals_LL_cell_Delta1{1,4} = 0.5; % line width

eigvals_LL_cell_Delta1{2,1} = eigvals_LL_Kp_Delta1;
eigvals_LL_cell_Delta1{2,2} = dims;
eigvals_LL_cell_Delta1{2,3} = 'b--';
eigvals_LL_cell_Delta1{2,4} = 1.0;

save_path = "";
fig2 = plot_LLs(Delta1_start, Delta1_end, - 0.02, 0.01, Delta1_list, eigvals_LL_cell_Delta1, save_path);
% fig3 = plot_LLs_weights(eigvec_HK_now(:,2 * LL_index_cutoff - 3 :2 * LL_index_cutoff + 3)); % 对本征态在LL basis下的权重作图

% 画出cnp曲线
% cnp_LL_index = dims;
% eigvals_LL_Delta1 = [eigvals_LL_K_Delta1, eigvals_LL_Kp_Delta1];
% plot_cnp_line(fig2, eigvals_LL_Delta1, Delta1_list, cnp_LL_index, 0.75);

% 画出triplet T2对应的三条LLs(全部来自于Kp valley)
% T3_LL_index_list = [round(dims / 2), round(dims / 2) + 1, round(dims / 2) + 2];
% figure(fig2)
% hold on
% plot(Delta1_list, eigvals_LL_Kp_Delta1(:, T3_LL_index_list), 'r--')

T2_LL_index_list = [round(dims / 2) - 2, round(dims / 2) - 1, round(dims / 2)];
T2_eigvals = eigvals_LL_K_Delta1(:, T2_LL_index_list);
figure(fig2)
hold on
plot(Delta1_list, T2_eigvals, 'r')

% Extended view of the triplet T2 with the average energy of the triplet substracted
ene_av_list = sum(T2_eigvals, 2) / 3;
fig4 = figure();
axis([Delta1_start Delta1_end -0.0025 0.0015])
hold on
for i = 1:3
    plot(Delta1_list, eigvals_LL_K_Delta1(:, T2_LL_index_list(i)) - ene_av_list)
end

%% plot cnp
flag_plot_cnp = false;
if flag_plot_cnp
    cnp_index1 = 181;
    cnp_index2 = 182;
    energy_cnp_list = zeros(Delta1_steps, 1);
    for D_index = Delta1_steps:-1:1
        eigvals_LL_sort = sort([eigvals_LL_K_Delta1(D_index, :), eigvals_LL_Kp_Delta1(D_index, :)], 'ascend');
        energy_cnp_list(D_index) = (eigvals_LL_sort(cnp_index1) + eigvals_LL_sort(cnp_index2)) / 2;
    end
    figure(fig2);
    hold on
    plot(Delta1_list, energy_cnp_list, 'g*', 'MarkerSize', 0.5);
end