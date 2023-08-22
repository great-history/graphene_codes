% add path
addpath("D:\matlab\graphene-package\trilayer-graphene\trilayer_landau_level")
%% Plot trilayer LL without D as a function of B
% 通过拟合得到的结果
gamma0 = 3100;    
gamma1 = 385.3;
gamma2 = -24.6;  
gamma3 = 247.8;    
gamma4 = 20.8;   
gamma5 = 51.5; 

delta = 40;
Delta1 = 0;
Delta2 = 1;

N_LL = 30;

B_start = 0.05;
B_end = 14;
B_steps = 1000;
B_fields = linspace(B_start, B_end, B_steps);

eps = 0.005 * gamma0;  % 取一个能够大致描述能隙大小的值
tic
[LL_K_m, LL_Kp_m, LL_K_b, LL_Kp_b] = trilayer_LLs_solver_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_fields, B_steps, eps);
toc

% figure
% % subplot(2,1,1)
% % axis([B_start B_end -100 100])
% axis([B_start B_end -120 120])
% hold on

dims_m = 2 * N_LL + 3;
dims_b = 4 * N_LL + 6;

% % plot bilayer-like LL
% for i = 1:dims_b
%     plot(B_fields, LL_K_b(:,i), 'b', 'LineWidth', 0.5)
%     plot(B_fields, LL_Kp_b(:,i), 'b--', 'LineWidth', 1.0)
% end
% 
% % plot monolayer-like LL
% for i = 1:dims_m
%     plot(B_fields, LL_K_m(:,i), 'r', 'LineWidth', 0.5)
%     plot(B_fields, LL_Kp_m(:,i), 'r--', 'LineWidth', 1.0)
% end
% grid on
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

figure
% subplot(2,1,1)
% axis([B_start B_end -100 100])
axis([B_start B_end -120 120])
hold on

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
grid off

hold on
plot(B_fields, LL_K_b(:, LLb0_K_index), 'm', 'LineWidth', 1.5) 
hold on
plot(B_fields, LL_Kp_b(:, LLb0_Kp_index), 'm', 'LineWidth', 1.5)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

hold on
plot(B_fields, LL_K_b(:, LL_K_b_positive_indexs(1)), 'm', 'LineWidth', 1.5)
hold on
plot(B_fields, LL_Kp_b(:, LL_Kp_b_positive_indexs(1)), 'm', 'LineWidth', 1.5)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

% hold on
% plot(B_fields, LL_K_b(:, LL_K_b_positive_indexs(2)), 'c', 'LineWidth', 0.25)
% hold on
% plot(B_fields, LL_Kp_b(:, LL_Kp_b_positive_indexs(2)), 'c', 'LineWidth', 0.25)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

hold on
plot(B_fields, LL_K_m(:, LLm0_K_index), 'g', 'LineWidth', 1.5)
hold on
plot(B_fields, LL_Kp_m(:, LLm0_Kp_index), 'g', 'LineWidth', 1.5)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差

% % 寻找交叉点
% LLm0与LLb1的交点，共有4个点
[B_cross, E_cross] = trilayer_LLs_find_cross_points(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, ...
    LL_K_b_positive_indexs(2), LL_Kp_b_positive_indexs(2), LLm0_K_index, LLm0_Kp_index, B_fields, B_steps);

% LLm0与LLb2的交点，共有4个点
[B_cross3, E_cross3] = trilayer_LLs_find_cross_points(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, ...
    LL_K_b_positive_indexs(3), LL_Kp_b_positive_indexs(3), LLm0_K_index, LLm0_Kp_index, B_fields, B_steps);

% LLm0与LLb3的交点，共有4个点
[B_cross4, E_cross4] = trilayer_LLs_find_cross_points(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, ...
    LL_K_b_positive_indexs(4), LL_Kp_b_positive_indexs(4), LLm0_K_index, LLm0_Kp_index, B_fields, B_steps);

% LLm0与LLb4的交点，共有4个点
[B_cross5, E_cross5] = trilayer_LLs_find_cross_points(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, ...
    LL_K_b_positive_indexs(5), LL_Kp_b_positive_indexs(5), LLm0_K_index, LLm0_Kp_index, B_fields, B_steps);

% sz = 5;
% scatter(B_cross, E_cross, sz, 'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor',[0 .7 .7], 'LineWidth',1.5)
% scatter(B_cross3, E_cross3, sz, 'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor',[0 .7 .7], 'LineWidth',1.5)
% scatter(B_cross4, E_cross4, sz, 'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor',[0 .7 .7], 'LineWidth',1.5)
% scatter(B_cross5, E_cross5, sz, 'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor',[0 .7 .7], 'LineWidth',1.5)

title('Landau Fan @ 1GPa','interpreter','latex', 'FontSize', 24);
ylabel('$E(meV)$','interpreter','latex', 'FontSize', 18);
xlabel('$B(T)$','interpreter','latex', 'FontSize', 18);
ylim([-10,20])

set(gcf,'position',[100,100, 600, 800]);
set(gca,'FontSize',16)  %是设置刻度字体大小
saveas(gcf,'LLL_1GPa.png');

% lowest bilayer LLs
% hold on
% plot(B_fields, LL_K_b_1_analytical, 'g*', 'LineWidth', 0.5)
% hold on
% plot(B_fields, LL_K_b_0 * ones(B_steps, 1), 'm*', 'LineWidth', 0.5)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差
% hold on
% plot(B_fields, LL_Kp_b_1_analytical, 'g*', 'LineWidth', 0.5)  % 对于bilayer，这些解析得到的低能LL一般不能与数值计算得到的结果完全吻合，会有1meV左右的偏差
% hold on
% plot(B_fields, LL_Kp_b_0 * ones(B_steps, 1), 'm*', 'LineWidth', 0.5)
% 
% % lowest monolayer LLs
% hold on
% plot(B_fields, LL_K_m_0 * ones(B_steps, 1), 'g*', 'LineWidth', 0.5)
% hold on
% plot(B_fields, LL_Kp_m_0 * ones(B_steps, 1), 'g*', 'LineWidth', 0.5)


%% 加入Zeeman splitting
% soc_ene = 9.27400949 / 1.6021766208 * 0.01;
% soc_enes = soc_ene * B_fields';
% 
% LL_K_b_spin_up = zeros(B_steps, dims_b);
% LL_K_b_spin_down = zeros(B_steps, dims_b);
% LL_Kp_b_spin_up = zeros(B_steps, dims_b);
% LL_Kp_b_spin_down = zeros(B_steps, dims_b);
% for i = 1:dims_b
%     LL_K_b_spin_up(:,i) = LL_K_b(:,i) - soc_enes;
%     LL_K_b_spin_down(:,i) = LL_K_b(:,i) + soc_enes;
%     LL_Kp_b_spin_up(:,i) = LL_Kp_b(:,i) - soc_enes;
%     LL_Kp_b_spin_down(:,i) = LL_Kp_b(:,i) + soc_enes;
% end
% 
% 
% LL_K_m_spin_up = zeros(B_steps, dims_m);
% LL_K_m_spin_down = zeros(B_steps, dims_m);
% LL_Kp_m_spin_up = zeros(B_steps, dims_m);
% LL_Kp_m_spin_down = zeros(B_steps, dims_m);
% 
% for i = 1:dims_m
%     LL_K_m_spin_up(:,i) = LL_K_m(:,i) - soc_enes;
%     LL_K_m_spin_down(:,i) = LL_K_m(:,i) + soc_enes;
%     LL_Kp_m_spin_up(:,i) = LL_Kp_m(:,i) - soc_enes;
%     LL_Kp_m_spin_down(:,i) = LL_Kp_m(:,i) + soc_enes;
% end
% 
% % plot bilayer-like LL
% for i = 1:dims_b
%     plot(B_fields, LL_K_b_spin_up(:,i), 'b', 'LineWidth', 0.5)
%     plot(B_fields, LL_Kp_b_spin_up(:,i), 'b--', 'LineWidth', 0.5)
%     plot(B_fields, LL_K_b_spin_down(:,i), 'b', 'LineWidth', 0.5)
%     plot(B_fields, LL_Kp_b_spin_down(:,i), 'b--', 'LineWidth', 0.5)
% end
% 
% % plot monolayer-like LL
% for i = 1:(2*N_LL+3)
%     plot(B_fields, LL_K_m_spin_up(:,i), 'r', 'LineWidth', 0.5)
%     plot(B_fields, LL_Kp_m_spin_up(:,i), 'r--', 'LineWidth', 0.5)
%     plot(B_fields, LL_K_m_spin_down(:,i), 'r', 'LineWidth', 0.5)
%     plot(B_fields, LL_Kp_m_spin_down(:,i), 'r--', 'LineWidth', 0.5)
% end

%% Plot trilayer LL with D as a function of B
% gamma0 = 3100;
% gamma1 = 390;
% gamma2 = -28;
% gamma3 = 315;  
% gamma4 = 41;
% gamma5 = 50;
% 
% delta = 46;
% Delta1_test = [0, 5, 10, 15, 20, 25, 30, 35];  % 0 / 25 / 50 / 100/ 150 / 180 / 250
% % Delta1_test = [0, 25, 50, 75, 100, 125, 150, 200];  % 0 / 25 / 50 / 100/ 150 / 180 / 250
% figure
% for test_idx = 1:length(Delta1_test)
%     Delta1 = Delta1_test(test_idx);
%     Delta1_start = 0;
%     Delta1_end = Delta1;
%     Delta1_steps = ceil(abs((Delta1_start - Delta1_end) / 0.1)) + 2;
%     Delta1s = linspace(Delta1_start, Delta1_end, Delta1_steps);
% 
%     Delta2 = 0;
%     N_LL = 30;
%     B_start = 0.05;
%     B_end = 6.05;
%     B_steps = 300;
% 
%     % B_start = 2.06;
%     % B_end = 3.06;
%     % B_steps = 10;
% 
%     B_fields = linspace(B_start, B_end, B_steps);
%     no_mix = true;
%     eps = 0.005 * gamma0;  % 取一个能够大致描述能隙大小的值
%     tic
%     [LL_K, LL_Kp] = trilayer_LLs_solver_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, N_LL, B_fields, B_steps, eps, no_mix, Delta1s, Delta1_steps);
%     toc
%     
%     subplot(2, 4, test_idx)
%     axis([B_start B_end -60 60])
%     hold on
% 
%     % plot monolayer-like LL
%     for i = 1:(2*N_LL+3)
%         plot(B_fields, LL_K(:,i), 'r', 'LineWidth', 0.5)
%         plot(B_fields, LL_Kp(:,i), 'r--', 'LineWidth', 1.5)
%     end
% 
%     % plot bilayer-like LL
%     for i = (2*N_LL+4):(6*N_LL+9)
%         plot(B_fields, LL_K(:,i), 'b', 'LineWidth', 0.5)
%         plot(B_fields, LL_Kp(:,i), 'b--', 'LineWidth', 1.5)
%     end
%     
%     title('D = ', num2str(Delta1,'%.3f'))
% end

%% Plot trilayer LL as a function of D
% % parameters adopted from PRB 87,085424
% gamma0 = 3100;
% gamma1 = 390;
% gamma2 = -28;
% gamma3 = 315;  
% gamma4 = 41;
% gamma5 = 50;
% 
% Delta1_start = 0;
% Delta1_end = 250;
% Delta1_steps = 500;
% Delta1s = linspace(Delta1_start, Delta1_end, Delta1_steps);
% 
% delta = 46;
% Delta2 = 0;
% N_LL = 30;
% B_field = 10;
% % B_field = 0.5;
% eps = 0.005 * gamma0;  % 取一个能够大致描述能隙大小的值
% 
% tic
% [LL_K, LL_Kp] = trilayer_LLs_asfo_Delta1(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_field, Delta1s, Delta1_steps, eps);
% toc
% 
% figure
% axis([Delta1_start Delta1_end -120 120])
% hold on
% 
% % plot LL
% % for i = 1:(6*N_LL+9)
% %     plot(Delta1s, LL_K(:,i), 'b-o', 'LineWidth', 0.5, 'MarkerSize',2)
% %     plot(Delta1s, LL_Kp(:,i), 'b--o', 'LineWidth', 0.5, 'MarkerSize',2)
% % end
% 
% for i = 1:(2*N_LL+3)
% %     plot(Delta1s, LL_K(:,i), 'r-o', 'LineWidth', 0.5, 'MarkerSize',5)
% %     plot(Delta1s, LL_Kp(:,i), 'r--o', 'LineWidth', 0.5, 'MarkerSize',5)
%     plot(Delta1s, LL_K(:,i), 'r-', 'LineWidth', 0.5, 'MarkerSize',5)
%     plot(Delta1s, LL_Kp(:,i), 'r--', 'LineWidth', 0.5, 'MarkerSize',5)
% end
% 
% for i = (2*N_LL+4):(6*N_LL+9)
% %     plot(Delta1s, LL_K(:,i), 'b-o', 'LineWidth', 0.5, 'MarkerSize',3)
% %     plot(Delta1s, LL_Kp(:,i), 'b--o', 'LineWidth', 0.5, 'MarkerSize',3)
%     plot(Delta1s, LL_K(:,i), 'b-', 'LineWidth', 0.5, 'MarkerSize',3)
%     plot(Delta1s, LL_Kp(:,i), 'b--', 'LineWidth', 0.5, 'MarkerSize',3)
% end