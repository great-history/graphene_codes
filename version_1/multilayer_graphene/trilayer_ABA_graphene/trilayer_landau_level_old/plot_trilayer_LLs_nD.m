% Plot trilayer LL as a function of D
% parameters adopted from PRB 87,085424
gamma0 = 3100;
gamma1 = 390;
gamma2 = -28;
gamma3 = 315;  
gamma4 = 41;
gamma5 = 50;

Delta1_start = 0;
Delta1_end = 250;
Delta1_steps = 500;
Delta1s = linspace(Delta1_start, Delta1_end, Delta1_steps);

delta = 46;
Delta2 = 0;
N_LL = 30;
B_field = 10;
% B_field = 0.5;
eps = 0.005 * gamma0;  % 取一个能够大致描述能隙大小的值

tic
[LL_K, LL_Kp] = trilayer_LLs_asfo_Delta1(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_field, Delta1s, Delta1_steps, eps);
toc

figure
axis([Delta1_start Delta1_end -120 120])
hold on

% plot LL
% for i = 1:(6*N_LL+9)
%     plot(Delta1s, LL_K(:,i), 'b-o', 'LineWidth', 0.5, 'MarkerSize',2)
%     plot(Delta1s, LL_Kp(:,i), 'b--o', 'LineWidth', 0.5, 'MarkerSize',2)
% end

dims_m = 2 * N_LL + 3;
dims_b = 4 * N_LL + 6;

for i = 1:dims_m
%     plot(Delta1s, LL_K(:,i), 'r-o', 'LineWidth', 0.5, 'MarkerSize',5)
%     plot(Delta1s, LL_Kp(:,i), 'r--o', 'LineWidth', 0.5, 'MarkerSize',5)
    plot(Delta1s, LL_K(:,i), 'r-', 'LineWidth', 0.5, 'MarkerSize',5)
    plot(Delta1s, LL_Kp(:,i), 'r--', 'LineWidth', 0.5, 'MarkerSize',5)
end

for i = (dims_m + 1):(dims_m + dims_b)
%     plot(Delta1s, LL_K(:,i), 'b-o', 'LineWidth', 0.5, 'MarkerSize',3)
%     plot(Delta1s, LL_Kp(:,i), 'b--o', 'LineWidth', 0.5, 'MarkerSize',3)
    plot(Delta1s, LL_K(:,i), 'b-', 'LineWidth', 0.5, 'MarkerSize',3)
    plot(Delta1s, LL_Kp(:,i), 'b--', 'LineWidth', 0.5, 'MarkerSize',3)
end

grid on

%% convert E-D at fixed B to n-D at fixed B
%% step1:找到Delta1 = 0时的lowest Landau levels，我们从B=0一直演化到当前的B值
[LL_K_m_0, LL_Kp_m_0, LL_K_b_0, LL_Kp_b_0, LL_K_b_1_slope, LL_K_b_1_intercept, LL_Kp_b_1_slope, LL_Kp_b_1_intercept] = ...
    trilayer_LLs_lowest(gamma0, gamma1, gamma2, gamma5, delta, Delta2);

% [LLb0_K_index, LLb0_Kp_index, LLm0_K_index, LLm0_Kp_index, LLb1_K_index, LLb1_Kp_index, LLb2_K_index, LLb2_Kp_index] = ...
%     trilayer_LLs_find_LLLs(LL_K(1,(dims_m + 1):(dims_m + dims_b)), LL_Kp(1,(dims_m + 1):(dims_m + dims_b)), ...
%                            LL_K(1,1:dims_m), LL_Kp(1,1:dims_m), LL_K_b_0, LL_Kp_b_0, LL_K_m_0, LL_Kp_m_0, N_LL, 1, 1);

% 找到LL_K_m_0/LL_Kp_m_0/LL_K_b_0/LL_K_b_1/LL_Kp_b_0/LL_Kp_b_1相应的指标
% 寻找LLm0_K和LLm0_Kp
gap_K = 100; % meV
gap_Kp = 100; % meV
LLm0_K_index = 0;
LLm0_Kp_index = 0;
for i = 1:dims_m
    gap_K_now = sum(abs(LL_K(1,i) - LL_K_m_0)) / 1;
    gap_Kp_now = sum(abs(LL_Kp(1,i) - LL_Kp_m_0)) / 1;
    
    if gap_K_now < gap_K
        gap_K = gap_K_now;
        LLm0_K_index = i;
    end

    if gap_Kp_now < gap_Kp
        gap_Kp = gap_Kp_now;
        LLm0_Kp_index = i;
    end
end

% 寻找LLb0_K和LLb0_Kp
gap_K = 100; % meV
gap_Kp = 100; % meV
LLb0_K_index = 0;
LLb0_Kp_index = 0;
for i = (dims_m + 1):(dims_m + dims_b)
    gap_K_now = sum(abs(LL_K(1,i) - LL_K_b_0)) / 1;
    gap_Kp_now = sum(abs(LL_Kp(1,i) - LL_Kp_b_0)) / 1;
    
    if gap_K_now < gap_K
        gap_K = gap_K_now;
        LLb0_K_index = i;
    end

    if gap_Kp_now < gap_Kp
        gap_Kp = gap_Kp_now;
        LLb0_Kp_index = i;
    end
end


% 寻找LLb1_K和LLb1_Kp，用解析的会偏差较大
gap_K = 100; % meV
gap_Kp = 100; % meV
LLb1_K_index = 0;
LLb1_Kp_index = 0;
for i = (dims_m + 1):(dims_m + dims_b)
    if ~(i == LLb0_K_index)
        gap_K_now = sum(abs(LL_K(1,i) - LL_K(1, LLb0_K_index))) / 1;
        if gap_K_now < gap_K
            gap_K = gap_K_now;
            LLb1_K_index = i;
        end
    end

    if ~(i == LLb0_Kp_index)
        gap_Kp_now = sum(abs(LL_Kp(1,i) - LL_Kp(1, LLb0_Kp_index))) / 1;
        if gap_Kp_now < gap_Kp
            gap_Kp = gap_Kp_now;
            LLb1_Kp_index = i;
        end
    end
end

% % 找出LLb2_K 和 LLb2_K
gap_K = 100; % meV
gap_Kp = 100; % meV
LLb2_K_index = 0;
LLb2_Kp_index = 0;
for i = (dims_m + 1):(dims_m + dims_b)
    if (LL_K(1,i) - LL_K(1, LLb1_K_index)) > 0 && ~(i == LLb0_K_index) && ~(i == LLb1_K_index)  % 保证斜率为正
        gap_K_now = sum(abs(LL_K(1, i) - LL_K(1, LLb1_K_index))) / 1;
        if gap_K_now < gap_K
            gap_K = gap_K_now;
            LLb2_K_index = i;
        end
    end

    if (LL_Kp(1,i) - LL_Kp(1, LLb1_Kp_index)) > 0 && ~(i == LLb0_Kp_index) && ~(i == LLb1_Kp_index)
        gap_Kp_now = sum(abs(LL_Kp(1, i) - LL_Kp(1, LLb1_Kp_index))) / 1;
        if gap_Kp_now < gap_Kp
            gap_Kp = gap_Kp_now;
            LLb2_Kp_index = i;
        end
    end
end

% print LLL index
% LLm0_K_index
% LLm0_Kp_index
% LLb0_K_index
% LLb0_Kp_index
% LLb1_K_index
% LLb1_Kp_index
% LLb2_K_index
% LLb2_Kp_index

hold on
plot(Delta1s, LL_K(:, LLb0_K_index), 'm*', 'LineWidth', 0.025, 'MarkerSize',2) 
hold on
plot(Delta1s, LL_Kp(:, LLb0_Kp_index), 'm*', 'LineWidth', 0.025, 'MarkerSize',2)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

hold on
plot(Delta1s, LL_K(:, LLb1_K_index), 'y*', 'LineWidth', 0.025, 'MarkerSize',2)
hold on
plot(Delta1s, LL_Kp(:, LLb1_Kp_index), 'y*', 'LineWidth', 0.025, 'MarkerSize',2)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

hold on
plot(Delta1s, LL_K(:, LLb2_K_index), 'c*', 'LineWidth', 0.25, 'MarkerSize',2)
hold on
plot(Delta1s, LL_Kp(:, LLb2_Kp_index), 'c*', 'LineWidth', 0.25, 'MarkerSize',2)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

hold on
plot(Delta1s, LL_K(:, LLm0_K_index), 'g*', 'LineWidth', 0.25, 'MarkerSize',2)
hold on
plot(Delta1s, LL_Kp(:, LLm0_Kp_index), 'g*', 'LineWidth', 0.25, 'MarkerSize',2)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

%% step2:确定CNP
% % plot DOS as a function of n & D
E_upper_bound = max(max(max(LL_K)), max(max(LL_Kp))); % 为了保证不同Delta1下的CNP都对应n=0，必须把整个能量范围都考虑
E_lower_bound = min(min(min(LL_K)), min(min(LL_Kp)));
delta_ene = 0.1;
E_points = ceil((E_upper_bound - E_lower_bound) / delta_ene); % 保证能量点之间的间隔远远小于Gamma

Gamma = 1.5; % 由于disorder导致的LL展宽

% get DOS as a function of Ene & D
tic
[ene_vecs, DOS_ED_matrix] = helper_get_dos_asfo_ED(LL_K, LL_Kp, E_upper_bound, E_lower_bound, E_points, Delta1_steps, B_field, Gamma);
toc

figure
imagesc(Delta1s, ene_vecs, DOS_ED_matrix); % B_fields是x轴, 与DOS_matrix的第二个指标对应；ene_vecs是y轴, 与DOS_matrix的第一个指标对应
colorbar

% get carrier density as a function of Ene & B
tic
[carrier_density_matrix] = helper_get_n_asfo_ED(DOS_ED_matrix, delta_ene, E_points, Delta1_steps);
toc
% carrier density shift to the CNP
% 首先确定出在Delta1 = 0时CNP所对应的能量值，然后对能量E进行积分后得到对应的carrier density shift
E_CNP = (LL_K(1, LLb0_K_index) + LL_K(1, LLb1_K_index)) / 2;
for E_index = 1:E_points
    if ene_vecs(E_index) <= E_CNP && ene_vecs(E_index + 1) >= E_CNP
        E_CNP_index = E_index;
        break
    end
end

carrier_density_shift = sum(DOS_ED_matrix(1:E_CNP_index, 1)) * delta_ene;
carrier_density_matrix = carrier_density_matrix - carrier_density_shift * ones(size(carrier_density_matrix));

figure
imagesc(Delta1s, ene_vecs, carrier_density_matrix); % B_fields是x轴, 与DOS_matrix的第二个指标对应；ene_vecs是y轴, 与DOS_matrix的第一个指标对应
colorbar

% %% get DOS as a function of n & D
% 首先确定出carrier_density的上限和下限
n_points = 1000;
carrier_density_ub = 3 * B_field / (25.66)^2;  % 最大磁场下六个朗道能级填满
carrier_density_lb = - 3 * B_field / (25.66)^2;
% carrier_density_ub = 0.05;
% carrier_density_lb = -0.05;

% 开始作插值
tic
[carrier_density_vecs, DOS_nD_matrix] = helper_get_LLs_nD(DOS_ED_matrix, carrier_density_matrix, carrier_density_ub, carrier_density_lb, n_points, Delta1_steps);
toc

figure
im = imagesc(Delta1s, carrier_density_vecs, DOS_nD_matrix); % B_fields是x轴, 与DOS_matrix的第二个指标对应；carrier_density_vecs是y轴, 与DOS_matrix的第一个指标对应
% im = imagesc(carrier_density_vecs, Delta1s, DOS_nD_matrix');
axis([Delta1_start Delta1_end -0.02 0.02])
colorbar
hold on

% 画几条填充线
fillings = [-4, -3, -2, -1, 0, 1, 2, 3, 4];