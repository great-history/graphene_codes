function [LL_K_b_positive_indexs, LL_K_b_negative_indexs, LL_Kp_b_positive_indexs, LL_Kp_b_negative_indexs, LL_K_m_positive_indexs, LL_K_m_negative_indexs, LL_Kp_m_positive_indexs, LL_Kp_m_negative_indexs] = ...
    trilayer_LLs_find_indexs(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, LLb0_K_index, LLb0_Kp_index, LLm0_K_index, LLm0_Kp_index)
    [~,I] = sort(LL_K_b(end,:));
    ii = find(I == LLb0_K_index);  % 说明比LLb0_K低的LL有(ii-1)个，而比LLb1_K高的LL有(dims_b - ii - 1)个
    LL_K_b_positive_indexs = I(ii+1:dims_b);
    LL_K_b_negative_indexs = I(1:ii);

    [~,I] = sort(LL_Kp_b(end,:));
    ii = find(I == LLb0_Kp_index);
    LL_Kp_b_positive_indexs = I(ii+1:dims_b);
    LL_Kp_b_negative_indexs = I(1:ii);

    [~,I] = sort(LL_K_m(end,:));
    ii = find(I == LLm0_K_index);
    LL_K_m_positive_indexs = I(ii+1:dims_m);
    LL_K_m_negative_indexs = I(1:ii);

    [~,I] = sort(LL_Kp_m(end,:));
    ii = find(I == LLm0_Kp_index);
    LL_Kp_m_positive_indexs = I(ii+1:dims_m);
    LL_Kp_m_negative_indexs = I(1:ii);
end