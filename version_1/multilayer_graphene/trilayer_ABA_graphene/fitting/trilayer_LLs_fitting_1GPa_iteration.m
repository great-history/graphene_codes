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

x0 = [3100, 385.7, -30, 300.0,   43.9,   53.6,  40, 1.0];
B_cross2_max_exp = 5.85;
B_cross2_min_exp = 5.15;
B_cross3_max_exp = 3.32;
B_cross3_min_exp = 2.86;
B_cross4_max_exp = 2.32;
B_cross4_min_exp = 1.97;
ydata = [B_cross2_max_exp, B_cross2_min_exp, B_cross3_max_exp, B_cross3_min_exp, B_cross4_max_exp, B_cross4_min_exp];
error = 0.01;
turn_on = true;
B_error = 0.2;

% x1 = iteration_fitting_lsq(x0, ydata, error, turn_on);
% Pablo的结果：% gamma0 = 3100; gamma1 = 390; gamma2 = -28; gamma3 = 315; gamma4 = 41; gamma5 = 50; delta = 46; Delta1 = 0; Delta2 = 1.7;

% [x1, x1_down, x1_up] = get_error_bar(x0, ydata, error, turn_on, B_error);
% 得到的拟合结果为      x1 = 3.1000    0.3855   -0.0279    0.2478    0.0376    0.0555    0.0399    0.0010
% 得到的拟合结果为 x1_down = 3.1000    0.3857   -0.0271    0.3000    0.0439    0.0558    0.0400    0.0010
% 得到的拟合结果为   x1_up = 3.1000    0.3858   -0.0292    0.2000    0.0361    0.0541    0.0398    0.0011
% 可以发现gamma3 / gamma4 变化较大

turn_on = false;
[x1, x1_down, x1_up] = get_error_bar(x0, ydata, error, turn_on, B_error);

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


%% 迭代得到最佳解
% 因为要拟合的参数有gamma1,gamma2,gamma3,gamma4,gamm5五个参数，我会先按照重要性将它们分为[gamma2, gamma5], [gamma3, gamma4], [gamma1]这三类，然后依次拟合，每次拟合其中一类固定另外两类知道收敛为止
function x1 = iteration_fitting_lsq(x0, ydata, error, turn_on)       
    % errs = 0.0001;  %  我希望我的误差在0.0001左右
    % x0_init是一开始猜的一组解，这个函数是用来反复迭代得到最佳解
    xdata = [1,2,3,4,5,6];
    
    sub1 = [0, 0, 20, 0, 0, 20, 0, 0];
    sub2 = [0, 0, 0, 20, 20, 0, 0, 0];
    sub3 = [0, 20, 0, 0, 0, 0, 0, 0];
    sub4 = [0, 0, 0, 0, 0, 0, 5, 3];
    
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','UseParallel', true);
    
    if turn_on == true
        % 第一组：变动gamma2/gamma5
        lb = x0 - sub1;
        ub = x0 + sub1;
        [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x0, xdata, ydata, lb, ub, options);

        % 第二组：变动gamma3/gamma4
        lb = x1 - sub2;
        ub = x1 + sub2;
        [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x1, xdata, ydata, lb, ub, options);

        % 第三组：变动gamma1
        lb = x1 - sub3;
        ub = x1 + sub3;
        [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x1, xdata, ydata, lb, ub, options);

        % 第四组：变动delta/Delta2
        lb = x1 - sub4;
        ub = x1 + sub4;
        [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x1, xdata, ydata, lb, ub, options);

        while max(x1 - x0) >= error
            x0 = x1;
            % 第一组：变动gamma2/gamma5
            lb = x0 - sub1;
            ub = x0 + sub1;
            [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x0, xdata, ydata, lb, ub, options);

            % 第二组：变动gamma3/gamma4
            lb = x1 - sub2;
            ub = x1 + sub2;
            [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x1, xdata, ydata, lb, ub, options);

            % 第三组：变动gamma1
            lb = x1 - sub3;
            ub = x1 + sub3;
            [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x1, xdata, ydata, lb, ub, options);

            % 第四组：变动delta/Delta2
            lb = x1 - sub4;
            ub = x1 + sub4;
            [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x1, xdata, ydata, lb, ub, options);
        end
    else
        % 第一组：变动gamma2/gamma5
        lb = x0 - sub1;
        ub = x0 + sub1;
        [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x0, xdata, ydata, lb, ub, options);

        % 第二组：变动gamma3/gamma4
        lb = x1 - sub2;
        ub = x1 + sub2;
        [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x1, xdata, ydata, lb, ub, options);

        % 第三组：变动gamma1
        lb = x1 - sub3;
        ub = x1 + sub3;
        [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x1, xdata, ydata, lb, ub, options);

        while max(x1 - x0) >= error
            x0 = x1;
            % 第一组：变动gamma2/gamma5
            lb = x0 - sub1;
            ub = x0 + sub1;
            [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x0, xdata, ydata, lb, ub, options);

            % 第二组：变动gamma3/gamma4
            lb = x1 - sub2;
            ub = x1 + sub2;
            [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x1, xdata, ydata, lb, ub, options);

            % 第三组：变动gamma1
            lb = x1 - sub3;
            ub = x1 + sub3;
            [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_1GPa, x1, xdata, ydata, lb, ub, options);
        end
    end
end

%% 计算error_bar:(来源于LL展宽导致的误差)
function [x1, x1_down, x1_up] = get_error_bar(x0, ydata, error, turn_on, B_error)
    % B_error:LL broadening导致的交叉点处磁场的误差
    ydata_down = ydata - B_error;
    ydata_up = ydata + B_error;
    
    num_error = 0.01;
    turn_on = true;
    
    
    x1 = iteration_fitting_lsq(x0, ydata, num_error, turn_on);
    x1_down = iteration_fitting_lsq(x0, ydata_down, num_error, turn_on);
    x1_up = iteration_fitting_lsq(x0, ydata_up, num_error, turn_on);
end