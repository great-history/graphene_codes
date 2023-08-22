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
% x0 = [3100, 390, -22.3, 315, 45.6, 46.4, 34.2, 0];
% lb = [3100;390;-45;315;20;35;30;0];
% ub = [3100;390;-15;315;80;105;55;0];

% 第二组
% x0 = [3100, 390, -25, 315, 41, 46.4, 40, 1.0];
% lb = [3100;390;-80;315;41;20;40;1.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub = [3100;390;-10;315;41;100;40;1.0];

% 第三组
% x0 = [3100, 390, -25, 315, 41, 46.4, 40, 1.0];
% lb = [3100;390;-80;315;41;20;40;1.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub = [3100;390;-10;315;41;100;40;1.0];
% 得到的拟合结果为   3.1000    0.3900   -0.0188    0.3150    0.0410    0.0549    0.0400    0.0010
% 来自1GPa的拟合结果 3.1000    0.3900   -0.0247    0.3150    0.0410    0.0517    0.0400    0.0010

% 第四组：来自1GPa的数据，因为1GPa的crossing point比较多，比较好拟合
% 从1GPa的拟合中我们发现gamma1的变化几乎可以忽略, gamma2/gamma5变化相对较小，gamm3/gamma4变化相对较大
% x0 = [3100, 385.3, -24.6, 247.8, 50.8, 51.5, 40, 1.0];
% lb = [3100; 380.3; -28.6; 247.8; 50.8; 40;   40; 1.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub = [3100; 395.3; -15.6; 247.8; 50.8; 60;   40; 1.0];
% 得到的拟合结果为 3.1000    0.3803   -0.0198    0.2478    0.0508    0.0534    0.0400    0.0010

% 第五组：由于gamma2/gamma5变化相对较小，因此固定gamma2/gamma5，变动gamma3/gamma4
% x0 = [3100, 385.3, -19.8, 246.3,   50.8,   53.4,  40, 1.0];
% lb = [3100; 365.3; -19.8; 185;     50.8;   53.4;  40; 1.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub = [3100; 395.3; -19.8; 225;     65.8;   53.4;  40; 1.0];
% 得到的拟合结果为     3.1000    0.3876   -0.0198    0.2057    0.0597    0.0534    0.0400    0.0010
% 而从1GPa得到的结果为 3.1000    0.3853   -0.0246    0.2478    0.0208    0.0515    0.0400    0.0010

% 第五组：固定gamma1/gamma3/gamma4，变动gamma2/gamma5
% x0 = [3100, 387.6, -19.8, 205.7,   59.7,   53.4,  40, 1.0];
% lb = [3100; 387.6; -30.8; 205.7;   59.7;   43.4;  40; 1.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub = [3100; 387.6; -10.8; 205.7;   59.7;   63.4;  40; 1.0];
% 得到的拟合结果为 3.1000    0.3876   -0.0196    0.2057    0.0597    0.0536    0.0400    0.0010, 这说明gamma2/gamma5差不多就是这个拟合的值了

% 第六组：固定gamma2/gamma5，变动gamma1/gamma3/gamma4
x0 = [3100, 387.6, -19.8, 205.7,   59.7,   53.6,  40, 1.0];
lb = [3100; 387.6; -19.8; 165.7;   49.7;   53.6;  40; 1.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
ub = [3100; 387.6; -19.8; 245.7;   69.7;   53.6;  40; 1.0];
% 得到的拟合结果为     3.1000    0.3876   -0.0198    0.2130    0.0639    0.0536   0.0400    0.0010 说明gamma3/gamma4也拟合的差不多了
% 而从1GPa得到的结果为 3.1000    0.3853   -0.0246    0.2478    0.0208    0.0515    0.0400    0.0010

% 第七组：变动gamma1
x0 = [3100, 387.6, -19.8, 213.0,   63.9,   53.6,  40, 1.0];
lb = [3100; 377.6; -19.8; 213.0;   63.9;   53.6;  40; 1.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
ub = [3100; 397.6; -19.8; 213.0;   63.9;   53.6;  40; 1.0];
% 得到的拟合结果为     3.1000    0.3877   -0.0198    0.2130    0.0639    0.0536  0.0400    0.0010 说明gamma1拟合的也差不多了
% 而从1GPa得到的结果为 3.1000    0.3853   -0.0246    0.2478    0.0208    0.0515    0.0400    0.0010
% 通过对比可以发现gamma1发生了很小的变化，gamma2/gamma5变化相对较小，但是它们是真正决定了低能部分的结构因此即使是2—5meV左右的变化也会很大地影响LLL，
% 而gamma3和gamma4变化较大，在30-40meV的量级，但是它们主要影响的是能量较高的部分，对低能影响不大
% 因此我们可以把加压当作是可以调控gamma2/gamma5较小变动的一种手段，通过加压实验可以有效地改变低能的结构

% 最后一组：refine,主要是细调delta2, delta
x0 = [3100, 387.7, -19.8, 213.0,   63.9,   53.6,  40, 1.0];
lb = [3100, 387.7, -19.8, 213.0,   63.9,   53.6,  30, 0.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
ub = [3100, 387.7, -19.8, 213.0,   63.9,   53.6,  50, 2.0];
% 得到的拟合结果为     3.1000    0.3877   -0.0198    0.2130    0.0639    0.0536   0.0401    0.0010, 没有发生变化

xdata = [1,2,3,4];

B_cross2_max_exp = 5.13;
B_cross2_min_exp = 3.96;
B_cross3_max_exp = 2.84;
B_cross3_min_exp = 2.32;
ydata = [B_cross2_max_exp, B_cross2_min_exp, B_cross3_max_exp, B_cross3_min_exp];

[x2,resnorm2,residual2,exitflag2,output2] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_0GPa, x0, xdata, ydata, lb, ub, options);

x0 = [3100, 387.7, -19.8, 213.0,   63.9,   53.6,  40, 1.0];
error = 0.0001;
x1 = iteration_fitting_lsq(x0_init, error);


%% 拟合函数
function ydata = trilayer_LLs_fitting_without_D_lsq_0GPa(gamma_params, xdata)
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
    
    ydata = [max(B_cross2), min(B_cross2), max(B_cross3), min(B_cross3)];
end


%% 迭代得到最佳解
% 因为要拟合的参数有gamma1,gamma2,gamma3,gamma4,gamm5五个参数，我会先按照重要性将它们分为[gamma2, gamma5], [gamma3, gamma4], [gamma1]这三类，然后依次拟合，每次拟合其中一类固定另外两类知道收敛为止
function x1 = iteration_fitting_lsq(x0_init, error)
    % errs = 0.0001;  %  我希望我的误差在0.0001左右
    % x0_init是一开始猜的一组解，这个函数是用来反复迭代得到最佳解
    B_cross2_max_exp = 5.13;
    B_cross2_min_exp = 3.96;
    B_cross3_max_exp = 2.84;
    B_cross3_min_exp = 2.32;
    ydata = [B_cross2_max_exp, B_cross2_min_exp, B_cross3_max_exp, B_cross3_min_exp];
    xdata = [1,2,3,4];
    
    sub1 = [0, 0, 20, 0, 0, 20, 0, 0];
    sub2 = [0, 0, 0, 20, 20, 0, 0, 0];
    sub3 = [0, 20, 0, 0, 0, 0, 0, 0];
    
    
    x0 = x0_init;
    % 第一组：变动gamma2/gamma5
    lb = x0 - sub1;
    ub = x0 + sub1;
    [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_0GPa, x0, xdata, ydata, lb, ub, options);
    
    % 第二组：变动gamma3/gamma4
    lb = x1 - sub2;
    ub = x1 + sub2;
    [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_0GPa, x1, xdata, ydata, lb, ub, options);
    
    % 第三组：变动gamma1
    lb = x1 - sub3;
    ub = x1 + sub3;
    [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_0GPa, x1, xdata, ydata, lb, ub, options);
    
    while max(x1 - x0) >= error
        x0 = x1;
        % 第一组：变动gamma2/gamma5
        lb = x0 - sub1;
        ub = x0 + sub1;
        [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_0GPa, x0, xdata, ydata, lb, ub, options);

        % 第二组：变动gamma3/gamma4
        lb = x1 - sub2;
        ub = x1 + sub2;
        [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_0GPa, x1, xdata, ydata, lb, ub, options);

        % 第三组：变动gamma1
        lb = x1 - sub3;
        ub = x1 + sub3;
        [x1,~,~,~,~] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq_0GPa, x1, xdata, ydata, lb, ub, options);
    end
end