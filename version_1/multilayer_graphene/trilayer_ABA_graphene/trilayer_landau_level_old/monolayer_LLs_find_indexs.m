function [LL_K_m_positive_indexs, LL_K_m_negative_indexs, LL_Kp_m_positive_indexs, LL_Kp_m_negative_indexs] = ...
    monolayer_LLs_find_indexs(LL_K_m, LL_Kp_m, LLm0_K_index, LLm0_Kp_index)
    dims_m = length(LL_K_m);
    
    [~,I] = sort(LL_K_m);
    ii = find(I == LLm0_K_index);
    LL_K_m_positive_indexs = I(ii+1:dims_m);
    LL_K_m_negative_indexs = I(1:ii);

    [~,I] = sort(LL_Kp_m);
    ii = find(I == LLm0_Kp_index);
    LL_Kp_m_positive_indexs = I(ii+1:dims_m);
    LL_Kp_m_negative_indexs = I(1:ii);
end