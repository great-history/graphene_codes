gamma0 = 3100;
gamma1 = 390;
gamma2 = -28;
gamma3 = 315;
gamma4 = 41;
gamma5 = 50;

delta = 46;
Delta1 = 0;
Delta2 = 5.7;

% monolayer
U1 = Delta2 - gamma2 / 2;
U2 = Delta2 + delta - gamma5 / 2;

N_LL = 30;

B_start = 0.01;
B_end = 14;
B_steps = 1000;
B_fields = linspace(B_start, B_end, B_steps);

eps = 0.005 * gamma0;  % 取一个能够大致描述能隙大小的值
tic
[LL_K_m, LL_Kp_m] = monolayer_LLs_solver(gamma0, U1, U2, N_LL, B_fields, B_steps, eps);
toc

figure
axis([B_start B_end -300 300])
hold on

dims_m = 2 * N_LL + 3;

% plot monolayer-like LL
for i = 1:dims_m
    plot(B_fields, LL_K_m(:,i), 'r', 'LineWidth', 0.5)
    plot(B_fields, LL_Kp_m(:,i), 'r--', 'LineWidth', 1.0)
end

grid on

[LL_K_m_0, LL_Kp_m_0] = monolayer_LLs_lowest(gamma0, U1, U2);
% % 找出数值计算中相应的Lowest LLs
start_index = floor(B_steps / 2);
end_index = B_steps;
[LLm0_K_index, LLm0_Kp_index] = monolayer_LLs_find_LLLs(LL_K_m, LL_Kp_m, LL_K_m_0, LL_Kp_m_0, N_LL, start_index, end_index);

% hold on
% plot(B_fields, LL_K_m(:, LLm0_K_index), 'g*', 'LineWidth', 0.25)
% hold on
% plot(B_fields, LL_Kp_m(:, LLm0_Kp_index), 'g*', 'LineWidth', 0.25)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

%% plot DOS as a function of n & B
E_upper_bound = 300;
E_lower_bound = -300;
E_points = 6000; % 保证能量点之间的间隔远远小于Gamma
delta_ene = (E_upper_bound - E_lower_bound) / E_points;


Gamma = 2; % 由于disorder导致的LL展宽
% get DOS as a function of Ene & B
[ene_vecs, DOS_EB_matrix] = helper_get_dos_asfo_EB(LL_K_m, LL_Kp_m, B_fields, E_upper_bound, E_lower_bound, E_points, B_steps, Gamma);

figure
imagesc(B_fields, ene_vecs, DOS_EB_matrix); % B_fields是x轴, 与DOS_matrix的第二个指标对应；ene_vecs是y轴, 与DOS_matrix的第一个指标对应
colorbar

% % get carrier density as a function of Ene & B
[carrier_density_matrix] = helper_get_n_asfo_EB(LL_K_m(:, LLm0_K_index), LL_Kp_m(:, LLm0_Kp_index), DOS_EB_matrix, ene_vecs, delta_ene, E_points, B_steps);

% figure
% imagesc(B_fields, ene_vecs, carrier_density_matrix); % B_fields是x轴, 与DOS_matrix的第二个指标对应；ene_vecs是y轴, 与DOS_matrix的第一个指标对应
% colorbar
% hold on
% for i = 1:dims_m
%     plot(B_fields, LL_K_m(:,i), 'r', 'LineWidth', 0.5)
%     plot(B_fields, LL_Kp_m(:,i), 'r--', 'LineWidth', 1.0)
% end

%% get DOS as a function of n & B
% 首先确定出carrier_density的上限和下限
n_points = 1000;
carrier_density_ub = 3 * B_fields(end) / (25.66)^2;  % 最大磁场下六个朗道能级填满
carrier_density_lb = - 3 * B_fields(end) / (25.66)^2;
% carrier_density_ub = 0.05;
% carrier_density_lb = -0.05;
% carrier_density_vecs = linspace(carrier_density_lb, carrier_density_ub, n_points);

% 开始作插值
[carrier_density_vecs, DOS_nB_matrix] = helper_get_LLs_nB(DOS_EB_matrix, carrier_density_matrix, carrier_density_ub, carrier_density_lb, n_points, B_steps);

figure
im = imagesc(B_fields, carrier_density_vecs, DOS_nB_matrix); % B_fields是x轴, 与DOS_matrix的第二个指标对应；carrier_density_vecs是y轴, 与DOS_matrix的第一个指标对应
colorbar