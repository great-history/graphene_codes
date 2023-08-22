%% 使用遗传算法进行拟合
% % 放在主函数中，主要是用来开启并行池，并将相应并行会使用的文件的路径加入
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     disp('启动并行运算，核心数：8');
%     % Perform a basic check by entering this code, where "local" is one kind of cluster profile.
%     parpool('local', 8);
% else
%     disp(['并行运o算已启动，核心数：' num2str(poolobj.NumWorkers)]);
% end
% 
% addAttachedFiles(poolobj, {'D:\matlab\graphene-package\trilayer-graphene\trilayer_landau_level\trilayer_LLs_fitting_without_D_ga'})
% 
% A=[];
% b=[];
% Aeq=[];
% beq=[];
% % Pablo parameters: gamma0 = 3100; gamma1 = 390; gamma2 = -28; gamma3 = 315; gamma4 = 41; gamma5 = 50; delta = 46;
% 
% lb = [3100;390;-36;315;40;40;30;0];  % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% ub = [3100;390;-14;315;50;52;40;1.5];
% % options = optimoptions('solvername','UseParallel',true); % The argument 'solvername' represents a nonlinear solver that supports parallel evaluation.
% options = optimoptions('ga','UseParallel', true);
% tic
% [opt_x,fval] = ga(@trilayer_LLs_fitting_without_D_ga,8,A,b,Aeq,beq,lb,ub,[],options);
% toc
% 
% % delete(poolobj)  % 用完之后，关闭并行池
% 
% % 3.1000    0.3900   -0.0206    0.3150    0.0437    0.0495    0.0382    0.0002
% % 3.1000    0.3900   -0.0223    0.3150    0.0500    0.0453    0.0305    0.0012
% % 3.1000    0.3900   -0.0183    0.3150    0.0456    0.0434    0.0342    0.0009
% % 3.1000    0.3900   -0.0224    0.3150    0.0458    0.0509    0.0334    0.0012


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

% gamma0 = 3100;
% gamma1 = 390;
% gamma2 = -28;
% gamma3 = 315;
% gamma4 = 41;
% gamma5 = 50;
% 
% delta = 46;
% Delta1 = 0;
% Delta2 = 1.7;
% x0 = [3100, 390, -28, 315, 41, 50, 46, 1];  % gamma0 / gamma1 / gamma2 / gamma3 / gamma4 / gamma5 / delta / Delta2
% x0 = [3100, 390, -20.6, 315, 43.7, 49.5, 38.2, 0.2];
% x0 = [3100, 390, -18.3, 315, 45.6, 43.4, 34.2, 0.9];
x0 = [3100, 390, -22.3, 315, 45.6, 46.4, 34.2, 0];
% 
lb = [3100;390;-45;315;20;35;30;0];
ub = [3100;390;-15;315;80;105;55;0];
% ub = [3100;390;-14;315;50;55;50;1.5];
% lb = [3100;390;-35;315;40;47;30;0];
% ub = [3100;390;-25;315;100;80;100;0];

xdata = [1,2,3,4];

B_cross2_max_exp = 5.13;
B_cross2_min_exp = 3.96;
B_cross3_max_exp = 2.84;
B_cross3_min_exp = 2.32;
ydata = [B_cross2_max_exp, B_cross2_min_exp, B_cross3_max_exp, B_cross3_min_exp];

[x2,resnorm2,residual2,exitflag2,output2] = lsqcurvefit(@trilayer_LLs_fitting_without_D_lsq, x0, xdata, ydata, lb, ub, options);