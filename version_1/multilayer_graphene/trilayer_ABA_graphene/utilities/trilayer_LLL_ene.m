function [LL_K_m_0, LL_Kp_m_0, LL_K_b_0, LL_Kp_b_0] = trilayer_LLL_ene(gamma2, gamma5, delta, Delta2)
    % 最低的几个朗道能级是可以被确定下来的(在忽略gamma3的假设下)
    
    % monolayer
    % 确定monolayer LL0+ ( + : K valley)
    LL_K_m_0 = Delta2 + delta - gamma5 / 2;
    
    % 确定monolayer LL0- ( - : K' valley)
    LL_Kp_m_0 = Delta2 - gamma2 / 2;
    
    % bilayer
    % 确定bilayer LL0+ ( + : K valley)
    LL_K_b_0 = - 2 * Delta2;
    
    % 确定bilayer LL0- ( - : Kp valley)
    LL_Kp_b_0 = gamma2 / 2 + Delta2;
end