%% 脚本
% 该脚本用来对不同压强下的朗道扇形图进行拟合
% SWMC parameters from different references
% gamma0  gamma1  gamma2  gamma3  gamma4  gamma5  delta  Delta2
%   3.1    0.39   -0.028   0.315   0.041   0.05   0.046   0.0
%   3.1    0.38   -0.021   0.29    0.141   0.05   0.0355  0.0035   % Andrea Young PRL
%   3.16   0.39   -0.020   0.315   0.044   0.038  0.037   0.0      % Pablo Nature Physics
%   2.92   0.27   -0.022   0.15    0.10    0.0063 0.0362  0.0
%   3.0    0.40    0.0     0.3     0.15    0.0    0.018   0.0
%   3.1    0.38   -0.0257  0.288   0.050   0.056  0.0475  0.0     % 我们的0GPa的拟合结果
%% addpath
addpath("..\trilayer_ABA_hamiltonian\")
addpath("..\utilities\")
%% 使用最小二乘法进行拟合
% 使用最小二乘法的优点是①可以知道标准差；②一般只有一个解，而ga则会存在local minima的问题会有好几个解
% 放在主函数中，主要是用来开启并行池，并将相应并行会使用的文件的路径加入
poolobj = gcp('nocreate');
if isempty(poolobj)
    disp('启动并行运算，核心数：6');
    % Perform a basic check by entering this code, where "local" is one kind of cluster profile.
    parpool('local', 6);
else
    disp(['并行运算已启动，核心数：' num2str(poolobj.NumWorkers)]);
end

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','UseParallel', true); % 使用并行


%% 输入参数，主要就是hopping parameters
% LLL之间的交叉主要受一些对低能物理重要的参数，有gamma2 / gamma5 / delta / Delta2 等等
% gamma1 / gamma3 / gamma4的影响相对较小
% gamma0 取定为3.1(gamma0的具体取值是多少无所谓，因为相对能量才是有意义的)
% x0 = [3.1, 0.39,  -0.020,   0.3,  0.044, 0.05,  0.04, 0.0035];
% % bound for 1GPa
% lb_list = [3.1; 0.375; -0.040;  0.21;  0.040; 0.035; 0.03; 0.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub_list = [3.1; 0.395; -0.015;  0.39;  0.075; 0.07;  0.05; 0.005];

%     3.1  0.381  -0.032   0.291  0.0491 0.0561 0.0440  0.0     % our fitting results @ P = 1GPa
%     3.1  0.390  -0.0213  0.30   0.0442 0.052  0.0366  0.0028  % our fitting results @ P = 1GPa（似乎更可信）

% bound for 0GPa
% lb_list = [3.1; 0.380; -0.030;  0.15;  0.050; 0.055; 0.03; 0.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub_list = [3.1; 0.390; -0.015;  0.25;  0.085; 0.08;  0.06; 0.005];
% lb_list = [3.1; 0.380; -0.030;  0.20;  0.060; 0.055; 0.03; 0.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub_list = [3.1; 0.400; -0.010;  0.20;  0.085; 0.08;  0.06; 0.007];
% 固定gamma4参数，其影响应该是这里最小的
lb_list = [3.1; 0.370; -0.035;  0.25;  0.045; 0.04;  0.03; 0.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
ub_list = [3.1; 0.400; -0.010;  0.45;  0.085; 0.08;  0.06; 0.005];

%% 输出参数，主要就是交叉点对应的磁场强度
% % crossing points @ P = 1GPa
% % 对应LL_m_0与LL_b_K_2和LL_b_Kp_2之间的交叉点位置
% B_cross2_max_exp = 5.85;
% B_cross2_min_exp = 5.15;
% % 对应LL_m_0与LL_b_K_3和LL_b_Kp_3之间的交叉点位置
% B_cross3_max_exp = 3.32;
% B_cross3_min_exp = 2.86;
% % 对应LL_m_0与LL_b_K_4和LL_b_Kp_4之间的交叉点位置
% B_cross4_max_exp = 2.32;
% B_cross4_min_exp = 1.97;

% crossing points @ P = 0GPa
% 对应LL_m_0与LL_b_K_2和LL_b_Kp_2之间的交叉点位置
B_cross2_max_exp = 4.80;
B_cross2_min_exp = 4.10;
% 对应LL_m_0与LL_b_K_3和LL_b_Kp_3之间的交叉点位置
B_cross3_max_exp = 2.70;
B_cross3_min_exp = 2.30;
% 对应LL_m_0与LL_b_K_4和LL_b_Kp_4之间的交叉点位置
B_cross4_max_exp = 1.85;
B_cross4_min_exp = 1.75;

% input_index_list存放的是每个交叉点对应的指标，分别是1，2，3，4，5，6
input_index_list = [1,2,3,4,5,6];
% output_list存放的是每个交叉点对应的B field 
output_list = [B_cross2_max_exp, B_cross2_min_exp, B_cross3_max_exp, B_cross3_min_exp, B_cross4_max_exp, B_cross4_min_exp];

%% 随机生成初始值x0，然后保存到变量input_value_array，得到的结果保存到output_value_list中
% x0 = [3.1, 0.385, -0.025,  0.20,  0.064, 0.06,  0.04, 0.0035];
num_test = 3;
input_value_array = zeros(num_test, 8);
output_value_array = zeros(num_test, 8);
for jj = 1:num_test % 这里就不要用parfor了，会报错
    input_value_list = get_random_input_list(ub_list, lb_list, 8);
    input_value_array(jj, :) = input_value_list;
    tic
    [fit_hopping_params, resnorm2, residual2, exitflag2, output2] = ...
                    lsqcurvefit(@trilayer_ABA_LLs_fitting_without_D_lsq, input_value_list, input_index_list, output_list, lb_list, ub_list, options);
    toc
    output_value_array(jj, :) = fit_hopping_params;
end

%% 保存数据
filename = '_0GPa_without_Delta1.mat';
save_path = ['D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\test_data\', datestr(datetime, 'yy-mm-dd-HH-MM-SS'), filename];
save(save_path, 'input_value_array', 'output_value_array');