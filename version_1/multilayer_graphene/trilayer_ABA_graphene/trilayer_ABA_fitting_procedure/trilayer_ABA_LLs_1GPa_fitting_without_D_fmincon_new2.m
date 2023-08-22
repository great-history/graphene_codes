function sse_crossing_points = trilayer_ABA_LLs_1GPa_fitting_without_D_fmincon_new2(hopping_params)
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
    
    delta_ene = hopping_params(7);
    delta = gamma5 / 2 - gamma2 / 2 + delta_ene;
    Delta2 = hopping_params(8);
    Delta1 = 0.0; 
    
    % 将gamma_i转换为v_i，同时某些gamma_i需要扩大sqrt(2)倍!!!!
    gamma1 = gamma1 * sqrt(2);
    v0 = sqrt(3) / 2 * gamma0 * (a_intralayer * 10^(-9)) / h_bar;
    v3 = sqrt(3) / 2 * sqrt(2) * gamma3 * (a_intralayer * 10^(-9)) / h_bar; % 不要忘了乘以sqrt(2) !!!
    v4 = sqrt(3) / 2 * sqrt(2) * gamma4 * (a_intralayer * 10^(-9)) / h_bar; % 不要忘了乘以sqrt(2) !!!
    
    % 磁场
    B_start = 1.0;
    B_end = 9;
    B_steps = 401;
    B_fields_list = linspace(B_start, B_end, B_steps);
    
    % 朗道能级参数
    LL_index_cutoff = 30;
    dims_m = 2 * LL_index_cutoff + 1;
    dims_b = 4 * LL_index_cutoff;
    dims = dims_b + dims_m;
    % energy window
    ene_ub = 0.1; % 100meV
    ene_lb = - 0.1; % -100meV
    
    %% 计算LL fan diagram
    [eigvals_LL_K, eigvals_LL_Kp, eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, eig_info_HK_m_select_cell, eig_info_HKp_m_select_cell] = ...
                trilayer_ABA_LL_solver_EB_without_Delta1_with_ene_win(v0, v3, v4, gamma1, delta, gamma2, gamma5, Delta2, B_fields_list, B_steps, LL_index_cutoff, dims_m, dims_b, ene_ub, ene_lb);
    
    %% 对selected_LL eig info进行重新排列（从大磁场开始）
    % ene of LLL
    [LL_K_m_0, LL_Kp_m_0, LL_K_b_0, LL_Kp_b_0] = trilayer_LLL_ene(gamma2, gamma5, delta, Delta2);
    % LL_K_b_1 = LL_K_b_0;
    % LL_Kp_b_1 = LL_Kp_b_0;

    % 如果不考虑gamma3的影响，bilayer-like branch的0th LL和monolayer-like branch的0th LL都是有解析形式的并且是完全极化的
    weight = 0.8; % 如果不考虑gamma3的影响，那么weight应该是1，因此我们只需要关心gamma3的大小即可，即weight实际上应该是gamma3的函数, weight(gamma3 = 0) = 1
    ene_eps = 0.005; % 5meV的误差
    
    % 先确定单层对应的LL_K_m0, LL_Kp_m0, 
    [eig_info_HK_m_select_cell, eig_info_HKp_m_select_cell] = sort_LLs_m(eig_info_HK_m_select_cell, eig_info_HKp_m_select_cell, ...
                                                                         LL_K_m_0, LL_Kp_m_0, LL_index_cutoff, B_steps, ene_eps, weight);

    % 寻找LL_K_b0和LL_K_b1的指标
    weight_K_b0 = weight;
    weight_Kp_b0 = weight;
    % flag_slope_K_b0 = false;
    % flag_slope_Kp_b0 = false;
    [eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, eigval_LL_K_b0_list, eigval_LL_Kp_b0_list] = ...
                        sort_LLs_b(eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, LL_K_b_0, LL_Kp_b_0, LL_index_cutoff, B_steps, ene_eps, weight_K_b0, weight_Kp_b0);
                    
    %% 确定了所有LL的index之后可以确定交叉点(crossing point)了
    % 确定 LL_K_b2, LL_K_b3, LL_K_b4, LL_Kp_b2, LL_Kp_b3, LL_Kp_b4 与 LL_K_m0, LL_Kp_m0 的交叉点的位置
    % 其中 LL_K_m0, LL_Kp_m0 是严格的直线，因此在考虑disorder之后就是一个矩形
    
    ene_width_disorder = 0.0008; % 假设有1meV的disorder带来的能量展开
    % 要求LL_Kp_m0与LL_K_m0之间的能量间隔应小于ene_width_disorder
    [LL_K_m0_crossing_points_list, LL_Kp_m0_crossing_points_list, poly_LL_K_m0, poly_LL_Kp_m0, eigval_LL_K_b_list, eigval_LL_Kp_b_list] = ...
            trilayer_ABA_LLs_find_cross_points(eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, LL_K_m_0, LL_Kp_m_0, ene_width_disorder, B_steps, B_fields_list);
    LL_m0_crossing_points_list = zeros(1, 6);
    % output_list = [B_cross2_max_exp, B_cross2_min_exp, B_cross3_max_exp, B_cross3_min_exp, B_cross4_max_exp, B_cross4_min_exp];
    for ii = 2:4
        cp_K_left = LL_K_m0_crossing_points_list(1, ii);
        cp_K_right = LL_K_m0_crossing_points_list(2, ii);
        cp_Kp_left = LL_Kp_m0_crossing_points_list(1, ii);
        cp_Kp_right = LL_Kp_m0_crossing_points_list(2, ii);

        LL_m0_crossing_points_list(2*(ii - 2) + 1) = min(cp_K_left, cp_Kp_left);
        LL_m0_crossing_points_list(2*(ii - 2) + 2) = max(cp_K_right, cp_Kp_right);
    end

    % crossing points @ P = 1GPa
    % 对应LL_m_0与LL_b_K_2和LL_b_Kp_2之间的交叉点位置
    B_cross2_max_exp = 5.70;
    B_cross2_min_exp = 5.30;
    % 对应LL_m_0与LL_b_K_3和LL_b_Kp_3之间的交叉点位置
    B_cross3_max_exp = 3.10;
    B_cross3_min_exp = 2.90;
    % 对应LL_m_0与LL_b_K_4和LL_b_Kp_4之间的交叉点位置
    B_cross4_max_exp = 2.30;
    B_cross4_min_exp = 2.15;
    B_crossing_points_list = [B_cross2_max_exp, B_cross2_min_exp, B_cross3_max_exp, B_cross3_min_exp, B_cross4_max_exp, B_cross4_min_exp];

    sse_crossing_points = sum((LL_m0_crossing_points_list - B_crossing_points_list).^2);
end