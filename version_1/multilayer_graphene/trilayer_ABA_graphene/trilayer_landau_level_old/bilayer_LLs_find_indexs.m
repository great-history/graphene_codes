function [LL_K_b_positive_indexs, LL_K_b_negative_indexs, LL_Kp_b_positive_indexs, LL_Kp_b_negative_indexs] = ...
    bilayer_LLs_find_indexs(LL_K_b, LL_Kp_b, LLb0_K_index, LLb0_Kp_index)
    dims_b = length(LL_K_b);

    [~,I] = sort(LL_K_b);
    ii = find(I == LLb0_K_index);  % 说明比LLb0_K低的LL有(ii-1)个，而比LLb1_K高的LL有(dims_b - ii - 1)个
    LL_K_b_positive_indexs = I(ii+1:dims_b);
    LL_K_b_negative_indexs = I(1:ii);

    [~,I] = sort(LL_Kp_b);
    ii = find(I == LLb0_Kp_index);
    LL_Kp_b_positive_indexs = I(ii+1:dims_b);
    LL_Kp_b_negative_indexs = I(1:ii);
end