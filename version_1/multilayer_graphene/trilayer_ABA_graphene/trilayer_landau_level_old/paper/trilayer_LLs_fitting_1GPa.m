% add path
addpath("D:\matlab\graphene-package\trilayer-graphene\trilayer_landau_level")

%% 使用最小二乘法进行拟合
% 使用最小二乘法的优点是①可以知道标准差；②一般只有一个解，而ga则会存在local minima的问题会有好几个解
% 放在主函数中，主要是用来开启并行池，并将相应并行会使用的文件的路径加入
poolobj = gcp('nocreate');
if isempty(poolobj)
    disp('启动并行运算，核心数：8');
    % Perform a basic check by entering this code, where "local" is one kind of cluster profile.
    parpool('local', 8);
else
    disp(['并行运o算已启动，核心数：' num2str(poolobj.NumWorkers)]);
end

% addAttachedFiles(poolobj, {'D:\matlab\graphene-package\trilayer-graphene\trilayer_landau_level\trilayer_LLs_fitting_without_D_lsq.m'})

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','UseParallel', true);

%% 进行拟合
% 第一组
% x0 = [3100, 390, -30.3, 315, 40.6, 46.4, 34.2, 0];
% lb = [3100;390;-45;315;20;35;30;0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub = [3100;390;-25;315;50;65;45;2.0];

% % 第二组：固定所有变量除了gamma2和gamma5，因为它们最重要
% x0 = [3100, 390, -25, 315, 41, 46.4, 40, 1.0];
% lb = [3100;390;-80;315;41;20;40;1.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub = [3100;390;-10;315;41;100;40;1.0];
% % 得到的拟合结果为 3.1000    0.3900   -0.0247    0.3150    0.0410    0.0517    0.0400    0.0010

% 第三组：固定gamma2和gamma5，可以变动gamma1, gamma3, gamma4，因为它们次重要
% x0 = [3100, 385.3, -24.7, 315, 41, 51.7, 40, 1.0];
% lb = [3100; 385.3; -24.7; 205; 11; 51.7; 40; 1.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub = [3100; 385.3; -24.7; 345; 41; 51.7; 40; 1.0];
% % 得到的拟合结果为 3.1000    0.3853   -0.0247    0.2478    0.0208    0.0517  0.0400    0.0010, 发现gamma3/gamma4其实变化还蛮大的， 而gamma1几乎固定不变，说明加压对gamma1的影响不大

% 第四组：固定gamma1, gamma3, gamma4，再次可以变动gamma2, gamma5
x0 = [3100, 385.3, -24.7, 247.8, 20.8, 51.7, 40, 1.0];
lb = [3100; 385.3; -30; 247.8; 20.8; 41.7; 40; 1.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
ub = [3100; 385.3; -20; 247.8; 20.8; 61.7; 40; 1.0];
% % 得到的拟合结果为 3.1000    0.3853   -0.0246    0.2478    0.0208    0.0515    0.0400    0.0010

% 综上可以发现gamma2和gamma5的变化量其实相对较小，这原因是第一层和第三层石墨烯之间的距离大约是2d，gamma1的变化也不明显，而gamma3和gamma4的变化非常明显

xdata = [1,2,3,4,5,6];

B_cross2_max_exp = 5.85;
B_cross2_min_exp = 5.15;
B_cross3_max_exp = 3.32;
B_cross3_min_exp = 2.86;
B_cross4_max_exp = 2.32;
B_cross4_min_exp = 1.97;
ydata = [B_cross2_max_exp, B_cross2_min_exp, B_cross3_max_exp, B_cross3_min_exp, B_cross4_max_exp, B_cross4_min_exp];

[x2,resnorm2,residual2,exitflag2,output2] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x0, xdata, ydata, lb, ub, options);

%% 拟合函数
function ydata = trilayer_LLs_fitting_without_D_lsq_1GPa(gamma_params, xdata)
    % 用最小二乘求解非线性曲线拟合（数据拟合）问题
    
    % 该函数是用于“遗传算法”的函数
    gamma0 = gamma_params(1);
    gamma1 = gamma_params(2);
    gamma2 = gamma_params(3);
    gamma3 = gamma_params(4);  
    gamma4 = gamma_params(5);
    gamma5 = gamma_params(6);

    delta = gamma_params(7);
    Delta2 = gamma_params(8);

    N_LL = 20;

    B_start = 2;
    B_end = 26;
    B_steps = 100;
    B_fields = linspace(B_start, B_end, B_steps);

    eps = 0.005 * gamma0;  % 取一个能够大致描述能隙大小的值
    [LL_K_m, LL_Kp_m, LL_K_b, LL_Kp_b] = trilayer_LLs_solver_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_fields, B_steps, eps);
    
    % 理论上计算出来的最低几个朗道能级
    [LL_K_m_0, LL_Kp_m_0, LL_K_b_0, LL_Kp_b_0, ~, ~, ~, ~] = ...
        trilayer_LLs_lowest(gamma0, gamma1, gamma2, gamma5, delta, Delta2);
    
    % 找出数值计算结果中的最低几个朗道能级
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
        
    % % 寻找交叉点
    % LL_b_2 与 LL_m_0的交叉点
    [B_cross2, ~] = trilayer_LLs_find_cross_points(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, ...
        LL_K_b_positive_indexs(2), LL_Kp_b_positive_indexs(2), LLm0_K_index, LLm0_Kp_index, B_fields, B_steps);
    % LL_b_3 与 LL_m_0的交叉点
    [B_cross3, ~] = trilayer_LLs_find_cross_points(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, ...
        LL_K_b_positive_indexs(3), LL_Kp_b_positive_indexs(3), LLm0_K_index, LLm0_Kp_index, B_fields, B_steps);
    % LL_b_4 与 LL_m_0的交叉点
    [B_cross4, ~] = trilayer_LLs_find_cross_points(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, ...
        LL_K_b_positive_indexs(4), LL_Kp_b_positive_indexs(4), LLm0_K_index, LLm0_Kp_index, B_fields, B_steps);
    
    ydata = [max(B_cross2), min(B_cross2), max(B_cross3), min(B_cross3), max(B_cross4), min(B_cross4)];
end