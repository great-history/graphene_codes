function error = trilayer_LLs_fitting_without_D_ga(gamma_params)
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
    B_end = 14;
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
    [B_cross2, E_cross] = trilayer_LLs_find_cross_points(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, ...
        LL_K_b_positive_indexs(2), LL_Kp_b_positive_indexs(2), LLm0_K_index, LLm0_Kp_index, B_fields, B_steps);
    % LL_b_3 与 LL_m_0的交叉点
    [B_cross3, ~] = trilayer_LLs_find_cross_points(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, ...
        LL_K_b_positive_indexs(3), LL_Kp_b_positive_indexs(3), LLm0_K_index, LLm0_Kp_index, B_fields, B_steps);
    
    % 与实验值进行作差
    % 0Gpa : 5.13/3.96/2.84/2.32T         % 1Gpa
    B_cross2_max_exp = 5.13;
    B_cross2_min_exp = 3.96;
    B_cross3_max_exp = 2.84;
    B_cross3_min_exp = 2.32;
    
    error1 = abs(B_cross2_max_exp - max(B_cross2));
    error2 = abs(B_cross2_min_exp - min(B_cross2));
    error3 = abs(B_cross3_max_exp - max(B_cross3));
    error4 = abs(B_cross3_min_exp - min(B_cross3));
    
    error = error1 + error2 + error3 + error4;
end