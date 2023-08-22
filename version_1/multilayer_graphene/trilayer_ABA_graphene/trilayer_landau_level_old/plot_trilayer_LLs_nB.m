% gamma0 = 3100;
% gamma1 = 390;
% gamma2 = -28;
% gamma3 = 315;
% gamma4 = 41;
% gamma5 = 50;
% 
% delta = 46;
% Delta1 = 0;
% Delta2 = 0;

gamma0 = 3100;    
gamma1 = 390;
gamma2 = -21.6;  
gamma3 = 315;    
gamma4 = 46.5;   
gamma5 = 44.1; 

delta = 35.8;
Delta1 = 0;
Delta2 = 0;

N_LL = 30;

B_start = 0.01;
B_end = 14;
B_steps = 1000;
B_fields = linspace(B_start, B_end, B_steps);

eps = 0.005 * gamma0;  % 取一个能够大致描述能隙大小的值
tic
[LL_K_m, LL_Kp_m, LL_K_b, LL_Kp_b] = trilayer_LLs_solver_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_fields, B_steps, eps);
toc

figure
% subplot(2,1,1)
axis([B_start B_end -80 80])
% axis([B_start B_end -0.08 0.08])
hold on

dims_m = 2 * N_LL + 3;
dims_b = 4 * N_LL + 6;

% plot bilayer-like LL
for i = 1:dims_b
    plot(B_fields, LL_K_b(:,i), 'b', 'LineWidth', 0.5)
    plot(B_fields, LL_Kp_b(:,i), 'b--', 'LineWidth', 1.0)
end

% plot monolayer-like LL
for i = 1:dims_m
    plot(B_fields, LL_K_m(:,i), 'r', 'LineWidth', 0.5)
    plot(B_fields, LL_Kp_m(:,i), 'r--', 'LineWidth', 1.0)
end

%% convert the E-B to n-B
% % 寻找charge neutrality point
[LL_K_m_0, LL_Kp_m_0, LL_K_b_0, LL_Kp_b_0, LL_K_b_1_slope, LL_K_b_1_intercept, LL_Kp_b_1_slope, LL_Kp_b_1_intercept] = ...
    trilayer_LLs_lowest(gamma0, gamma1, gamma2, gamma5, delta, Delta2);

LL_K_b_1_analytical = LL_K_b_1_slope * B_fields + LL_K_b_1_intercept;
LL_Kp_b_1_analytical = LL_Kp_b_1_slope * B_fields + LL_Kp_b_1_intercept;

% % 找出数值计算中相应的Lowest LLs
start_index = floor(B_steps / 2);
end_index = B_steps;

slope_eps = 0.001;
gap_error = 1.5;
[LLb0_K_index, LLb0_Kp_index, LLm0_K_index, LLm0_Kp_index] = ...
    trilayer_LLs_find_LLLs(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, LL_K_b_0, LL_Kp_b_0, LL_K_m_0, LL_Kp_m_0, N_LL, start_index, end_index, slope_eps, gap_error);                                                                            

[LL_K_b_positive_indexs, LL_K_b_negative_indexs, LL_Kp_b_positive_indexs, LL_Kp_b_negative_indexs] = ...
    bilayer_LLs_find_indexs(LL_K_b(end,:), LL_Kp_b(end,:), LLb0_K_index, LLb0_Kp_index);

[LL_K_m_positive_indexs, LL_K_m_negative_indexs, LL_Kp_m_positive_indexs, LL_Kp_m_negative_indexs] = ...
    monolayer_LLs_find_indexs(LL_K_m(end,:), LL_Kp_m(end,:), LLm0_K_index, LLm0_Kp_index);

LLb1_K_index = LL_K_b_positive_indexs(1);
LLb1_Kp_index = LL_Kp_b_positive_indexs(1);
LLb2_K_index = LL_K_b_positive_indexs(2);
LLb2_Kp_index = LL_Kp_b_positive_indexs(2);

hold on
plot(B_fields, LL_K_b(:, LLb0_K_index), 'm*', 'LineWidth', 0.25) 
hold on
plot(B_fields, LL_Kp_b(:, LLb0_Kp_index), 'm*', 'LineWidth', 0.25)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

hold on
plot(B_fields, LL_K_b(:, LLb1_K_index), 'y*', 'LineWidth', 0.25)
hold on
plot(B_fields, LL_Kp_b(:, LLb1_Kp_index), 'y*', 'LineWidth', 0.25)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

hold on
plot(B_fields, LL_K_b(:, LLb2_K_index), 'c*', 'LineWidth', 0.25)
hold on
plot(B_fields, LL_Kp_b(:, LLb2_Kp_index), 'c*', 'LineWidth', 0.25)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

hold on
plot(B_fields, LL_K_m(:, LLm0_K_index), 'g*', 'LineWidth', 0.25)
hold on
plot(B_fields, LL_Kp_m(:, LLm0_Kp_index), 'g*', 'LineWidth', 0.25)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

%% plot DOS as a function of n & B
% %% plot DOS as a function of n & B
E_upper_bound = 300;
E_lower_bound = -300;
E_points = 6000; % 保证能量点之间的间隔远远小于Gamma
delta_ene = (E_upper_bound - E_lower_bound) / E_points;

Gamma = 1.5; % 由于disorder导致的LL展宽
% get DOS as a function of Ene & B
tic
[ene_vecs, DOS_EB_matrix_b] = helper_get_dos_asfo_EB(LL_K_b, LL_Kp_b, B_fields, E_upper_bound, E_lower_bound, E_points, B_steps, Gamma);
[~, DOS_EB_matrix_m] = helper_get_dos_asfo_EB(LL_K_m, LL_Kp_m, B_fields, E_upper_bound, E_lower_bound, E_points, B_steps, Gamma);
DOS_EB_matrix = DOS_EB_matrix_b + DOS_EB_matrix_m;
toc

figure
imagesc(B_fields, ene_vecs, DOS_EB_matrix); % B_fields是x轴, 与DOS_matrix的第二个指标对应；ene_vecs是y轴, 与DOS_matrix的第一个指标对应
colorbar

% % get carrier density as a function of Ene & B
tic
[carrier_density_matrix] = helper_get_n_asfo_EB(LL_K_b(:, LLb0_K_index), LL_K_b(:, LLb1_K_index), DOS_EB_matrix, ene_vecs, delta_ene, E_points, B_steps);
toc

figure
imagesc(B_fields, ene_vecs, carrier_density_matrix); % B_fields是x轴, 与DOS_matrix的第二个指标对应；ene_vecs是y轴, 与DOS_matrix的第一个指标对应
colorbar
hold on
for i = 1:dims_m
    plot(B_fields, LL_K_m(:,i), 'r', 'LineWidth', 0.5)
    plot(B_fields, LL_Kp_m(:,i), 'r--', 'LineWidth', 1.0)
end
for i = 1:dims_b
    plot(B_fields, LL_K_b(:,i), 'b', 'LineWidth', 0.5)
    plot(B_fields, LL_Kp_b(:,i), 'b--', 'LineWidth', 1.0)
end

%% get DOS as a function of n & B
% 首先确定出carrier_density的上限和下限
n_points = 1000;
carrier_density_ub = 3 * B_fields(end) / (25.66)^2;  % 最大磁场下六个朗道能级填满
carrier_density_lb = - 3 * B_fields(end) / (25.66)^2;
% carrier_density_ub = 0.05;
% carrier_density_lb = -0.05;
% carrier_density_vecs = linspace(carrier_density_lb, carrier_density_ub, n_points);

% 开始作插值
tic
[carrier_density_vecs, DOS_nB_matrix] = helper_get_LLs_nB(DOS_EB_matrix, carrier_density_matrix, carrier_density_ub, carrier_density_lb, n_points, B_steps);
toc

figure
im = imagesc(B_fields, carrier_density_vecs, DOS_nB_matrix); % B_fields是x轴, 与DOS_matrix的第二个指标对应；carrier_density_vecs是y轴, 与DOS_matrix的第一个指标对应
axis([2 14 -0.02 0.02])
colorbar
hold on

% 画几条填充线
fillings = [-4, -3, -2, -1, 0, 1, 2, 3, 4];
for index = 1:length(fillings)
    n_vecs = 1 / (pi * (25.66)^2) * fillings(index) * B_fields;
    plot(B_fields, n_vecs, 'm*', 'LineWidth', 0.1)
end
