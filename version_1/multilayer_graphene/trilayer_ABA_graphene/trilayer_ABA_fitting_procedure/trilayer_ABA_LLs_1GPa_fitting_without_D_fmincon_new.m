function sse_crossing_points = trilayer_ABA_LLs_1GPa_fitting_without_D_fmincon_new(hopping_params)
    %% 使用fmincon进行参数的优化
    % hopping_params存放各种跃迁参数  //  input_index_list存放的是每个交叉点对应的指标，分别是1，2，3，4，5，6
    % 用最小二乘求解非线性曲线拟合（数据拟合）问题
    %% 准备参数
    % 基本参数
    h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
    d_interlayer = 0.335; % 单位是nm
    a_intralayer = 0.246; % 单位是nm
    
    % 导入所有的hopping parameters
    gamma0 = hopping_params(1);
    gamma1 = hopping_params(2);
    gamma2 = hopping_params(3);
    gamma3 = hopping_params(4);  
    gamma4 = hopping_params(5);
    gamma5 = hopping_params(6);
    
    delta = hopping_params(7);
    Delta2 = hopping_params(8);
    
    % delta_ene
    delta_ene = 0.002;
    hopping_params1 = hopping_params;
    hopping_params1(7) = gamma5 / 2  - gamma2 / 2 + delta_ene;
    sse_crossing_points_1 = trilayer_ABA_LLs_1GPa_fitting_without_D_fmincon(hopping_params1);
    hopping_params2 = hopping_params;
    hopping_params2(7) = gamma5 / 2  - gamma2 / 2 - delta_ene;
    sse_crossing_points_2 = trilayer_ABA_LLs_1GPa_fitting_without_D_fmincon(hopping_params2);
    
    sse_crossing_points = min(sse_crossing_points_1, sse_crossing_points_2);
end