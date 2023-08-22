function [LL_K_m_0, LL_Kp_m_0, LL_K_b_0, LL_Kp_b_0, LL_K_b_1_slope, LL_K_b_1_intercept, LL_Kp_b_1_slope, LL_Kp_b_1_intercept] = trilayer_LLs_lowest(gamma0, gamma1, gamma2, gamma5, delta, Delta2)
    % 最低的几个朗道能级是可以被确定下来的(在忽略gamma3的假设下)
    
    % monolayer
    U5 = Delta2 - gamma2 / 2;
    U6 = Delta2 + delta - gamma5 / 2;

    % 确定monolayer LL0+ ( + : K valley)
    LL_K_m_0 = U6;
    
    % 确定monolayer LL0- ( - : K' valley)
    LL_Kp_m_0 = U5;
    
    % bilayer
    U1 = Delta2 + gamma2 / 2;
    U2 = Delta2 + delta + gamma5 / 2;
    U3 = -2*Delta2 + delta;
    U4 = -2*Delta2;
    
    % 确定bilayer LL0+ ( + : K valley)
    LL_K_b_0 = U4;
    
    % 确定bilayer LL1+ ( + : K valley)
    LL_K_b_1_slope = (3 / 4) * (0.246 / 25.66)^2 * (gamma0 / gamma1)^2 * U2;
    LL_K_b_1_intercept = U4;
    
    % 确定monolayer LL0- ( - : K' valley)
    LL_Kp_b_0 = U1;
    
    % 确定bilayer LL1- ( - : K' valley)
    LL_Kp_b_1_slope = (3 / 4) * (0.246 / 25.66)^2 * (gamma0 / gamma1)^2 * U3;
    LL_Kp_b_1_intercept = U1;
end