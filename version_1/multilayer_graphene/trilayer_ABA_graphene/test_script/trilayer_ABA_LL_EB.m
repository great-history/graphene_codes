%% 该脚本计算的是三层ABA石墨烯六能带模型下的朗道能级图
% 该脚本分为两个部分
% 第一部分：Landau level as a function of E & B
% 第二部分：Landau level as a function of E & D @ fixed B

%% 添加路径
addpath(".\trilayer_ABA_hamiltonian\")
addpath(".\plot_funcs\")
addpath(".\utilities\")

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
%   2.92   0.27   -0.022   0.15    0.10    0.0063 0.0362  0.0
%   3.0    0.40    0.0     0.3     0.15    0.0    0.018   0.0
%
%   3.1    0.381  -0.032   0.291   0.0491  0.0561 0.0440  0.0     % our fitting results @ P = 1GPa
%   3.1    0.390  -0.0213  0.30    0.0442  0.052  0.0366  0.0028  % our fitting results @ P = 1GPa （似乎更可信）
%
%   3.1    0.383  -0.0234  0.21    0.0750  0.0700 0.0467  0.0008  % our fitting results @ P = 0GPa（似乎不可信）
%   3.1    0.384  -0.0254  0.200   0.0646  0.0600 0.0425  0.000   % our fitting results @ P = 0GPa（似乎更可信）
%   3.1    0.394  -0.0161  0.2087  0.0846  0.0795 0.0478  0.0037   % our fitting results @ P = 0GPa（似乎更可信）
% 
%     3.1000    0.3700   -0.0241    0.2500    0.0850    0.0787    0.0514    0.0013
%     3.1000    0.3700   -0.0185    0.2500    0.0543    0.0470    0.0327    0.0025
%     3.1000    0.3700   -0.0184    0.2500    0.0455    0.0470    0.0327    0.0024
%     3.1000    0.3700   -0.0293    0.3420    0.0850    0.0786    0.0539         0
%     3.1000    0.3700   -0.0219    0.2844    0.0640    0.0435    0.0327    0.0012
%     3.1000    0.3701   -0.0154    0.2500    0.0754    0.0787    0.0471    0.0050
%     3.1000    0.3700   -0.0184    0.2500    0.0455    0.0470    0.0327    0.0024
%     3.1000    0.3700   -0.0241    0.2500    0.0850    0.0787    0.0514    0.0013
%     3.1000    0.3700   -0.0184    0.2500    0.0455    0.0470    0.0327    0.0024
%     3.1000    0.3958   -0.0224    0.3287    0.0719    0.0695    0.0460    0.0010
%     3.1000    0.3858   -0.0292    0.2000    0.0361    0.0541    0.0398    0.0011
gamma0 = 3.1; 
gamma1 = 0.3858;
gamma2 = -0.0292;
gamma3 = 0.2000;  % 通过测试发现gamma3在[-0.29, 0.49]的范围内对LL_K_b_0和LL_Kp_b_0对应的本征态和本征能量的影响是可以降低的
gamma4 = 0.0361;
gamma5 = 0.0541;
delta = 0.0398;
Delta2 = 0.0011;

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

B_start = 1.0;
B_end = 9;
B_steps = 401;
B_fields_list = linspace(B_start, B_end, B_steps);

%% 第一部分：Landau level as a function of E & B
%% 朗道能级参数
LL_index_cutoff = 30;
dims_m = 2 * LL_index_cutoff + 1;
dims_b = 4 * LL_index_cutoff;
dims = dims_b + dims_m;
% energy window
ene_ub = 0.1; % 100meV
ene_lb = - 0.1; % -100meV

%% 计算LL fan diagram
tic
[eigvals_LL_K, eigvals_LL_Kp, eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, eig_info_HK_m_select_cell, eig_info_HKp_m_select_cell] = ...
            trilayer_ABA_LL_solver_nB_without_Delta1_with_ene_win(v0, v3, v4, gamma1, delta, gamma2, gamma5, Delta1, Delta2, B_fields_list, B_steps, LL_index_cutoff, dims_m, dims_b, ene_ub, ene_lb);
toc  

%% 作图
flag_plot_LL = true;
if flag_plot_LL
    % 单层双层全部画(颜色一样)
    % eigvals_LL_cell = cell(2, 4);
    % eigvals_LL_cell{1,1} = eigvals_LL_K;
    % eigvals_LL_cell{1,2} = dims; % 维数 dims : 单层双层全部画
    % eigvals_LL_cell{1,3} = 'b--'; % line color
    % eigvals_LL_cell{1,4} = 0.5; % line width
    % 
    % eigvals_LL_cell{2,1} = eigvals_LL_Kp;
    % eigvals_LL_cell{2,2} = dims;
    % eigvals_LL_cell{2,3} = 'b';
    % eigvals_LL_cell{2,4} = 1.0;

    % 单层双层全部画(颜色不同)
    eigvals_LL_cell = cell(4, 4); % mono_K, bi_K, mono_Kp, bi_Kp
    eigvals_LL_cell{1,1} = eigvals_LL_K(:, 1:dims_m);
    eigvals_LL_cell{1,2} = dims_m; % 维数 dims : 单层双层全部画
    eigvals_LL_cell{1,3} = 'b'; % line color
    eigvals_LL_cell{1,4} = 0.5; % line width

    eigvals_LL_cell{2,1} = eigvals_LL_Kp(:, 1:dims_m);
    eigvals_LL_cell{2,2} = dims_m; % 维数 dims : 单层双层全部画
    eigvals_LL_cell{2,3} = 'b--'; % line color
    eigvals_LL_cell{2,4} = 0.5; % line width

    eigvals_LL_cell{3,1} = eigvals_LL_K(:, (dims_m + 1):dims);
    eigvals_LL_cell{3,2} = dims_b;
    eigvals_LL_cell{3,3} = 'r';
    eigvals_LL_cell{3,4} = 0.5;

    eigvals_LL_cell{4,1} = eigvals_LL_Kp(:, (dims_m + 1):dims);
    eigvals_LL_cell{4,2} = dims_b;
    eigvals_LL_cell{4,3} = 'r--';
    eigvals_LL_cell{4,4} = 0.5;

    % % 只画单层
    % eigvals_LL_cell = cell(2, 4);
    % eigvals_LL_cell{1,1} = eigvals_LL_K(:, 1:dims_m); % eigvals_LL_K(B_index, 1:dims_m) = eigval_HK_m_diag_now;
    % eigvals_LL_cell{1,2} = dims_m;
    % eigvals_LL_cell{1,3} = 'b--'; % line color
    % eigvals_LL_cell{1,4} = 0.5; % line width
    % 
    % eigvals_LL_cell{2,1} = eigvals_LL_Kp(:, 1:dims_m);
    % eigvals_LL_cell{2,2} = dims_m;
    % eigvals_LL_cell{2,3} = 'b';
    % eigvals_LL_cell{2,4} = 1.0;

    % % 只画双层
    % eigvals_LL_cell = cell(1, 4);
    % eigvals_LL_cell{1,1} = eigvals_LL_K(:, (dims_m + 1):dims); 
    % eigvals_LL_cell{1,2} = dims_b;
    % eigvals_LL_cell{1,3} = 'b--'; % line color
    % eigvals_LL_cell{1,4} = 0.5; % line width
    % 
    % eigvals_LL_cell{2,1} = eigvals_LL_Kp(:, (dims_m + 1):dims);
    % eigvals_LL_cell{2,2} = dims_b;
    % eigvals_LL_cell{2,3} = 'b';
    % eigvals_LL_cell{2,4} = 1.0;
    % 
    save_path = "";
    E_bottom = -0.100;
    E_top = 0.100;
    fig0 = plot_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eigvals_LL_cell, save_path);
    % fig1 = plot_LLs_weights(eigvec_HK_now(:,2 * LL_index_cutoff - 3 :2 * LL_index_cutoff + 3)); % 对本征态在LL basis下的权重作图
    % fig2 = plot_select_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eig_info_HK_b_select_cell, "r", save_path);
    % fig3 = plot_select_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eig_info_HKp_b_select_cell, "r", save_path);
    % fig4 = plot_select_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eig_info_HK_m_select_cell, "b", save_path);
    % fig5 = plot_select_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eig_info_HKp_m_select_cell, "b", save_path);
end

%% 对selected_LL eig info进行重新排列（从大磁场开始）
% ene of LLL
[LL_K_m_0, LL_Kp_m_0, LL_K_b_0, LL_Kp_b_0] = trilayer_LLL_ene(gamma2, gamma5, delta, Delta2);
LL_K_b_1 = LL_K_b_0;
LL_Kp_b_1 = LL_Kp_b_0;
% 
% LL_K_m0_index_cell = cell(B_steps, 1);
% LL_Kp_m0_index_cell = cell(B_steps, 1);
% LL_K_b0_index_cell = cell(B_steps, 1);
% LL_K_b1_index_cell = cell(B_steps, 1);
% LL_Kp_b0_index_cell = cell(B_steps, 1);
% LL_Kp_b1_index_cell = cell(B_steps, 1);
% 
% % 如果不考虑gamma3的影响，bilayer-like branch的0th LL和monolayer-like branch的0th LL都是有解析形式的并且是完全极化的
weight = 0.8; % 如果不考虑gamma3的影响，那么weight应该是1，因此我们只需要关心gamma3的大小即可，即weight实际上应该是gamma3的函数, weight(gamma3 = 0) = 1
ene_eps = 0.005; % 5meV的误差

% 先确定单层对应的LL_K_m0, LL_Kp_m0, 
[eig_info_HK_m_select_cell, eig_info_HKp_m_select_cell] = sort_LLs_m(eig_info_HK_m_select_cell, eig_info_HKp_m_select_cell, ...
                                                                     LL_K_m_0, LL_Kp_m_0, LL_index_cutoff, B_steps, ene_eps, weight);

% 寻找LL_K_b0和LL_K_b1的指标
weight_K_b0 = weight;
weight_Kp_b0 = weight;
flag_slope_K_b0 = false;
flag_slope_Kp_b0 = false;
[eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, eigval_LL_K_b0_list, eigval_LL_Kp_b0_list] = ...
                    sort_LLs_b(eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, LL_K_b_0, LL_Kp_b_0, LL_index_cutoff, B_steps, ene_eps, weight_K_b0, weight_Kp_b0);
                
% 以cnp的能量设置为能量原点E=0（energy shifiting）
flag_plot_cnp = true;
if flag_plot_cnp
    % 找到最大磁场下的cnp的指标和对应的能量
    energy_cnp_list = zeros(B_steps, 1);

    index1 = eig_info_HK_b_select_cell{B_steps, 4} == 0;
    eigval_temp1 = eig_info_HK_b_select_cell{B_steps, 2}(index1);
    index2 = eig_info_HK_b_select_cell{B_steps, 4} == 1;
    eigval_temp2 = eig_info_HK_b_select_cell{B_steps, 2}(index2);
    energy_cnp_list(B_steps) = (eigval_temp1 + eigval_temp2) / 2;

    eigvals_LL_sort = sort([eigvals_LL_K(B_steps, :), eigvals_LL_Kp(B_steps, :)], 'ascend');
    index1 = find(eigvals_LL_sort == eigval_temp1);
    index2 = find(eigvals_LL_sort == eigval_temp2);

    for B_index = (B_steps - 1):-1:1
        eigvals_LL_sort = sort([eigvals_LL_K(B_index, :), eigvals_LL_Kp(B_index, :)], 'ascend');
        energy_cnp_list(B_index) = (eigvals_LL_sort(index1) + eigvals_LL_sort(index2)) / 2;
    end

    figure(fig0);
    hold on
    plot(B_fields_list, energy_cnp_list, 'g*', 'MarkerSize', 0.5);

    % 单层双层全部画(颜色不同)
    eigvals_LL_K_new = zeros(B_steps, dims);
    eigvals_LL_Kp_new = zeros(B_steps, dims);
    for B_index = 1:B_steps
        eigvals_LL_K_new(B_index, :) = eigvals_LL_K(B_index, :) - energy_cnp_list(B_index);
        eigvals_LL_Kp_new(B_index, :) = eigvals_LL_Kp(B_index, :) - energy_cnp_list(B_index);
    end
    
    eigvals_LL_cell = cell(4, 4); % mono_K, bi_K, mono_Kp, bi_Kp   
    eigvals_LL_cell{1,1} = eigvals_LL_K_new(:, 1:dims_m);
    eigvals_LL_cell{1,2} = dims_m; % 维数 dims : 单层双层全部画
    eigvals_LL_cell{1,3} = 'b'; % line color
    eigvals_LL_cell{1,4} = 0.5; % line width

    eigvals_LL_cell{2,1} = eigvals_LL_Kp_new(:, 1:dims_m);
    eigvals_LL_cell{2,2} = dims_m; % 维数 dims : 单层双层全部画
    eigvals_LL_cell{2,3} = 'b--'; % line color
    eigvals_LL_cell{2,4} = 0.5; % line width

    eigvals_LL_cell{3,1} = eigvals_LL_K_new(:, (dims_m + 1):dims);
    eigvals_LL_cell{3,2} = dims_b;
    eigvals_LL_cell{3,3} = 'r';
    eigvals_LL_cell{3,4} = 0.5;

    eigvals_LL_cell{4,1} = eigvals_LL_Kp_new(:, (dims_m + 1):dims);
    eigvals_LL_cell{4,2} = dims_b;
    eigvals_LL_cell{4,3} = 'r--';
    eigvals_LL_cell{4,4} = 0.5;
    
    save_path = "";
    E_bottom = -0.100;
    E_top = 0.100;
    fig7 = plot_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eigvals_LL_cell, save_path);
    
    % 计算态密度(DOS)
    ene_width = 0.0008; % 1meV的LL broadening
    Ene_steps = round((ene_ub - ene_lb) / ene_width) * 2 + 1;
    energy_list = linspace(ene_lb, ene_ub, Ene_steps);
    delta_ene = abs(energy_list(1) - energy_list(2));
    
    for B_index = 1:B_steps
        eig_info_HK_m_select_cell{B_index, 5} = eig_info_HK_m_select_cell{B_index, 2} - energy_cnp_list(B_index);
        eig_info_HKp_m_select_cell{B_index, 5} = eig_info_HKp_m_select_cell{B_index, 2} - energy_cnp_list(B_index);
        eig_info_HK_b_select_cell{B_index, 5} = eig_info_HK_b_select_cell{B_index, 2} - energy_cnp_list(B_index);
        eig_info_HKp_b_select_cell{B_index, 5} = eig_info_HKp_b_select_cell{B_index, 2} - energy_cnp_list(B_index);
    end
    
    dos_HK_m_EB_mat = get_dos_asfo_EB(eig_info_HK_m_select_cell, energy_list, Ene_steps, B_fields_list, B_steps, ene_width);
    dos_HKp_m_EB_mat = get_dos_asfo_EB(eig_info_HKp_m_select_cell, energy_list, Ene_steps, B_fields_list, B_steps, ene_width);
    dos_HK_b_EB_mat = get_dos_asfo_EB(eig_info_HK_b_select_cell, energy_list, Ene_steps, B_fields_list, B_steps, ene_width);
    dos_HKp_b_EB_mat = get_dos_asfo_EB(eig_info_HKp_b_select_cell, energy_list, Ene_steps, B_fields_list, B_steps, ene_width);
    
    dos_EB_mat = dos_HK_m_EB_mat + dos_HKp_m_EB_mat + dos_HK_b_EB_mat + dos_HKp_b_EB_mat;
    [density_EB_mat, density_max, density_min] = get_n_asfo_EB(dos_EB_mat, energy_list, Ene_steps, B_steps, delta_ene);
    
    density_steps = 200;
    density_list = linspace(density_min, density_max, density_steps);
    dos_nB_mat = get_dos_asfo_nB(dos_EB_mat, density_EB_mat, density_list, B_steps, density_steps);
    
    figure
    subplot(2,2,1)
    imagesc(B_fields_list, energy_list, dos_HK_m_EB_mat');
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
    subplot(2,2,2)
    imagesc(B_fields_list, energy_list, dos_HKp_m_EB_mat');
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
    subplot(2,2,3)
    imagesc(B_fields_list, energy_list, dos_HK_b_EB_mat');
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
    subplot(2,2,4)
    imagesc(B_fields_list, energy_list, dos_HKp_b_EB_mat');
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
    
    figure
    subplot(2,2,1)
    imagesc(B_fields_list, energy_list, dos_EB_mat');
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
    subplot(2,2,2)
    imagesc(B_fields_list, energy_list, density_EB_mat');
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
    subplot(2,2,3)
    % imagesc(B_fields_list, density_list, dos_nB_mat');
    imagesc(density_list, B_fields_list, dos_nB_mat);
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
end


%% 验证是否找对了
figure
hold on
plot(B_fields_list, eigval_LL_Kp_b0_list)
plot(B_fields_list, eigval_LL_K_b0_list)

%% 确定了所有LL的index之后可以确定交叉点(crossing point)了
% 确定 LL_K_b2, LL_K_b3, LL_K_b4, LL_Kp_b2, LL_Kp_b3, LL_Kp_b4 与 LL_K_m0, LL_Kp_m0 的交叉点的位置
% 其中 LL_K_m0, LL_Kp_m0 是严格的直线，因此在考虑disorder之后就是一个矩形
ene_width_disorder = 0.002; % 假设有2meV的disorder带来的能量展开

[LL_K_m0_crossing_points_list, LL_Kp_m0_crossing_points_list, poly_LL_K_m0, poly_LL_Kp_m0, eigval_LL_K_b_list, eigval_LL_Kp_b_list] = ...
            trilayer_ABA_LLs_find_cross_points(eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, LL_K_m_0, LL_Kp_m_0, ene_width_disorder, B_steps, B_fields_list);
figure
hold on
plot(poly_LL_K_m0)
plot(poly_LL_Kp_m0)
for ii = 1:4
    plot(B_fields_list, eigval_LL_K_b_list(ii, :), 'r');
    plot(B_fields_list, eigval_LL_Kp_b_list(ii, :), 'r--');
    
    plot(LL_K_m0_crossing_points_list(1, ii), LL_K_m_0 - ene_width_disorder, 'g*') % left crossing points
    plot(LL_K_m0_crossing_points_list(2, ii), LL_K_m_0 + ene_width_disorder, 'g*') % left crossing points
    
    plot(LL_Kp_m0_crossing_points_list(1, ii), LL_Kp_m_0 - ene_width_disorder, 'b*') % left crossing points
    plot(LL_Kp_m0_crossing_points_list(2, ii), LL_Kp_m_0 + ene_width_disorder, 'b*') % left crossing points
end