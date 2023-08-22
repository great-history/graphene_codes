%% 添加路径
addpath('.\trilayer_ABA_hamiltonian\')
addpath('.\utilities\')
addpath('.\plot_funcs\')

%% parameters set up
%% 基本参数
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
one_electron_coulomb = 1.602176634 * 10^(-19); % 单位是C
epsilon_0 = 8.85 * 10^(-12); % 单位是F/m
epsilon_bn = 6.6;
% d_interlayer = 0.335; % 单位是nm
d_interlayer = 0.100; % 单位是nm
a_intralayer = 0.246; % 单位是nm

%% 能带参数(以eV为单位)
% SWMC parameters from different references
% gamma0  gamma1  gamma2  gamma3  gamma4  gamma5  delta  Delta2
%   3.1    0.39   -0.028   0.315   0.041   0.05   0.046   0.0
%   3.1    0.38   -0.021   0.29    0.141   0.05   0.0355  0.0035   % Andrea Young PRL
%   3.16   0.39   -0.020   0.315   0.044   0.038  0.037   0.0      % Pablo Nature Physics
%   3.1    0.38   -0.0257  0.288   0.050   0.056  0.0475  0.0     % 我们的0GPa的拟合结果
%   3.1000    0.3858   -0.0292    0.2000    0.0361    0.0541    0.0398    0.0011
%   [3.10000000000000,0.384475107006668,-0.0193696369488659,0.350924488875967,0.0652882505455456,0.0636471402508554,0.0415083869026046,0.00208064832974783]
%   [3.10000000000000,0.384550640770226,-0.0235656429959901,0.332671140262836,0.0653751020545148,0.0623899653524236,0.0429778037649064,0.00252607886890637]
% gamma0 = 3.1; 
% gamma1 = 0.384;
% gamma2 = -0.021;
% gamma3 = 0.29;
% gamma4 = 0.141;
% gamma5 = 0.05;
% delta = 0.0355;
% Delta2 = 0.0035;

gamma0 = 3.1; 
gamma1 = 0.384;
gamma2 = -0.0236;
gamma3 = 0.333;  % 通过测试发现gamma3在[-0.29, 0.49]的范围内对LL_K_b_0和LL_Kp_b_0对应的本征态和本征能量的影响是可以降低的
gamma4 = 0.0654;
gamma5 = 0.0624;
delta = 0.0430;
Delta2 = 0.0025;

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
mag_length = 25.66 / sqrt(B_field);  % 以nm为单位
x0 = h_bar * v0 / mag_length * 10^9; % 以eV为单位
x3 = h_bar * v3 / mag_length * 10^9; % 以eV为单位
x4 = h_bar * v4 / mag_length * 10^9; % 以eV为单位

Delta1_steps = 241;
Delta1_start = 0.0;
Delta1_end = 0.4 * d_interlayer; % 与实验数据匹配
Delta1_list = linspace(Delta1_start, Delta1_end, Delta1_steps); % 以eV为单位

eigvals_LL_K_Delta1 = zeros(Delta1_steps, dims);
eigvals_LL_Kp_Delta1 = zeros(Delta1_steps, dims);

% energy window
ene_ub = 0.05; % 50meV
ene_lb = - 0.05; % -50meV

%% 计算LL fan diagram as E & D
tic
[eigvals_LL_K_ED, eigvals_LL_Kp_ED, eig_info_HK_ED_select_cell, eig_info_HKp_ED_select_cell] = ...
            trilayer_ABA_LL_solver_nD_fixed_B_with_ene_win(v0, v3, v4, gamma1, delta, gamma2, gamma5, Delta2, Delta1_list, Delta1_steps, B_field, LL_index_cutoff, dims_m, dims_b, ene_ub, ene_lb);
toc

%% 作图
eigvals_LL_cell_Delta1 = cell(2, 4);
eigvals_LL_cell_Delta1{1,1} = eigvals_LL_K_ED;
eigvals_LL_cell_Delta1{1,2} = dims; % 维数
eigvals_LL_cell_Delta1{1,3} = 'b'; % line color
eigvals_LL_cell_Delta1{1,4} = 0.5; % line width

eigvals_LL_cell_Delta1{2,1} = eigvals_LL_Kp_ED;
eigvals_LL_cell_Delta1{2,2} = dims;
eigvals_LL_cell_Delta1{2,3} = 'b--';
eigvals_LL_cell_Delta1{2,4} = 1.0;

save_path = "";
fig0 = plot_LLs(Delta1_start, Delta1_end, ene_lb, ene_ub, Delta1_list, eigvals_LL_cell_Delta1, save_path);

%% plot cnp
flag_plot_cnp = true;
if flag_plot_cnp
    cnp_index1 = 181;
    cnp_index2 = 182;
    energy_cnp_list = zeros(Delta1_steps, 1);
    for D_index = Delta1_steps:-1:1
        eigvals_LL_sort = sort([eigvals_LL_K_ED(D_index, :), eigvals_LL_Kp_ED(D_index, :)], 'ascend');
        energy_cnp_list(D_index) = (eigvals_LL_sort(cnp_index1) + eigvals_LL_sort(cnp_index2)) / 2;
    end
    
    figure(fig0);
    hold on
    plot(Delta1_list, energy_cnp_list, 'g*', 'MarkerSize', 0.5);
    
    eigvals_LL_K_new = zeros(Delta1_steps, dims);
    eigvals_LL_Kp_new = zeros(Delta1_steps, dims);
    for D_index = 1:Delta1_steps
        eigvals_LL_K_new(D_index, :) = eigvals_LL_K_ED(D_index, :) - energy_cnp_list(D_index);
        eigvals_LL_Kp_new(D_index, :) = eigvals_LL_Kp_ED(D_index, :) - energy_cnp_list(D_index);
    end
    
    eigvals_LL_cell = cell(2, 4);
    eigvals_LL_cell{1,1} = eigvals_LL_K_new;
    eigvals_LL_cell{1,2} = dims; % 维数 dims : 单层双层全部画
    eigvals_LL_cell{1,3} = 'b'; % line color
    eigvals_LL_cell{1,4} = 0.5; % line width
    
    eigvals_LL_cell{2,1} = eigvals_LL_Kp_new;
    eigvals_LL_cell{2,2} = dims;
    eigvals_LL_cell{2,3} = 'b--';
    eigvals_LL_cell{2,4} = 0.5;
    
    save_path = "";
    E_bottom = ene_lb;
    E_top = ene_ub;
    fig7 = plot_LLs(Delta1_start, Delta1_end, E_bottom, E_top, Delta1_list, eigvals_LL_cell, save_path);
    
    % 计算态密度(DOS)
    ene_width = 0.0008; % 1meV左右的LL broadening
    Ene_steps = round((ene_ub - ene_lb) / ene_width) * 2 + 1;
    energy_list = linspace(ene_lb, ene_ub, Ene_steps);
    delta_ene = abs(energy_list(1) - energy_list(2));
    
    % energy shifting with regard to cnp
    for D_index = 1:Delta1_steps
        eig_info_HK_ED_select_cell{D_index, 5} = eig_info_HK_ED_select_cell{D_index, 2} - energy_cnp_list(D_index);
        eig_info_HKp_ED_select_cell{D_index, 5} = eig_info_HKp_ED_select_cell{D_index, 2} - energy_cnp_list(D_index);
    end
    
    % 计算态密度(DOS) as a function of E & D
    dos_HK_ED_mat = get_dos_asfo_ED(eig_info_HK_ED_select_cell, energy_list, Ene_steps, Delta1_list, Delta1_steps, B_field, ene_width);
    dos_HKp_ED_mat = get_dos_asfo_ED(eig_info_HKp_ED_select_cell, energy_list, Ene_steps, Delta1_list, Delta1_steps, B_field, ene_width);
 
    dos_ED_mat = dos_HK_ED_mat + dos_HKp_ED_mat;
    
    % 计算载流子浓度(n) as a function of E & D
    [density_ED_mat, density_max, density_min] = get_n_asfo_ED(dos_ED_mat, energy_list, Ene_steps, Delta1_steps, delta_ene);
    
    % 计算态密度(DOS) as a function of n & D
    density_steps = 200;
    density_list = linspace(density_min, density_max, density_steps);
    dos_nD_mat = get_dos_asfo_nD(dos_ED_mat, density_ED_mat, density_list, Delta1_steps, density_steps);
    
    % 作图check
    figure
    subplot(1,2,1)
    imagesc(Delta1_list, energy_list, dos_HK_ED_mat');
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
    subplot(1,2,2)
    imagesc(Delta1_list, energy_list, dos_HKp_ED_mat');
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
    
    figure
    subplot(1,3,1)
    imagesc(Delta1_list, energy_list, dos_ED_mat');
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
    subplot(1,3,2)
    imagesc(Delta1_list, energy_list, density_ED_mat');
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
    subplot(1,3,3)
    imagesc(density_list, Delta1_list, dos_nD_mat);
    set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
    colorbar
end