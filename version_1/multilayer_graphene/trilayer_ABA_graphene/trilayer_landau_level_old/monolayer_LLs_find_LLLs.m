function [LLm0_K_index, LLm0_Kp_index] = monolayer_LLs_find_LLLs(LL_K_m, LL_Kp_m, LL_K_m_0, LL_Kp_m_0, N_LL, start_index, end_index)
    % % 找出数值计算中相应的Lowest LLs
    dims_m = 2 * N_LL + 3;
    % 寻找LLm0_K和LLm0_Kp
    gap_K = 100; % meV
    gap_Kp = 100; % meV
    LLm0_K_index = 0;
    LLm0_Kp_index = 0;
    for i = 1:dims_m
        gap_K_now = sum(abs(LL_K_m(start_index:end_index,i) - LL_K_m_0)) / (end_index - start_index);
        gap_Kp_now = sum(abs(LL_Kp_m(start_index:end_index,i) - LL_Kp_m_0)) / (end_index - start_index);
        if gap_K_now < gap_K
            gap_K = gap_K_now;
            LLm0_K_index = i;
        end

        if gap_Kp_now < gap_Kp
            gap_Kp = gap_Kp_now;
            LLm0_Kp_index = i;
        end
    end
end