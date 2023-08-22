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
addpath(".\trilayer_ABA_hamiltonian\")
addpath(".\utilities\")

%% 开启并行
poolobj = gcp('nocreate');
if isempty(poolobj)
    disp('启动并行运算，核心数：6');
    % Perform a basic check by entering this code, where "local" is one kind of cluster profile.
    parpool('local', 6);
else
    disp(['并行运算已启动，核心数：' num2str(poolobj.NumWorkers)]);
end

%% fmincon
ene_width_disorder = 0.0008;

options = [];
% A = [0, 0, -0.5, 0, 0, 0.5, -1, 0; 0, 0, 0.5, 0, 0, -0.5, 1, 0];
% b = [2 * ene_width_disorder; 2 * ene_width_disorder];
% Aeq = [];
% beq = [];
% A = [];
% b = [];
% Aeq = [0, 0, -0.5, 0, 0, 0.5, -1, 0];
% beq = 0.002; % 单位是eV  在B=4.5T时，Zeeman Energy约为0.5meV, 
% beq = 0.002; % 单位是eV  在B=4.5T时，Zeeman Energy约为0.5meV, 
Aeq = [];
beq = [];
A = [];
b = [];

% lb_list = [3.1; 0.370; -0.035;  0.25;  0.045; 0.04;  0.03; 0.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub_list = [3.1; 0.400; -0.010;  0.45;  0.085; 0.08;  0.06; 0.005];
lb_list = [3.1; 0.360; -0.035;  0.25;  0.046; 0.06;  0.04; 0.0]; % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
ub_list = [3.1; 0.400; -0.010;  0.55;  0.086; 0.10;  0.07; 0.005];

options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'MaxIterations', 1500, 'UseParallel', true); % 'interior-point'

num_test = 2;
input_value_array = zeros(num_test, 8);
output_value_array = zeros(num_test, 8);
for ii = 1:num_test
    initial_value_list = get_random_input_list(ub_list, lb_list, 8);
    input_value_array(ii, :) = initial_value_list; 
    tic
    % 对于fmincon没有的项一定要用空格站位
    [output_value_list, fval, exitflag, output, lambda, grad, hessian] = fmincon(@trilayer_ABA_LLs_1GPa_fitting_without_D_fmincon, initial_value_list, A, b, Aeq, beq, lb_list, ub_list, @nonlcon, options);
    toc
    output_value_array(ii, :) = output_value_list;
end

%% 保存数据
filename = '_1GPa_without_Delta1.mat';
save_path = ['D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\test_data\', datestr(datetime, 'yy-mm-dd-HH-MM-SS'), filename];
save(save_path, 'input_value_array', 'output_value_array');

function [c,ceq] = nonlcon(hopping_params)
    c = [];
    ceq = [];
end