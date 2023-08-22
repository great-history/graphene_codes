function [LL_K_m_0, LL_Kp_m_0] = monolayer_LLs_lowest(gamma0, U1, U2)
    % 最低的几个朗道能级是可以被确定下来的(在忽略gamma3的假设下)

    % 确定monolayer LL0+ ( + : K valley)
    LL_K_m_0 = U2;
    
    % 确定monolayer LL0- ( - : K' valley)
    LL_Kp_m_0 = U1;
end