function [LLb0_K_index, LLb0_Kp_index, LLm0_K_index, LLm0_Kp_index] = ...
    trilayer_LLs_find_LLLs(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, LL_K_b_0, LL_Kp_b_0, LL_K_m_0, LL_Kp_m_0, N_LL, start_index, end_index, slope_eps, gap_error)
    % 确定出在Delta1 = 0时最低的几个朗道能级：monolayer的0 和 bilayer的0
    
    % % 找出数值计算中相应的Lowest LLs
    dims_m = 2 * N_LL + 3;
    dims_b = 4 * N_LL + 6;
    
    % 寻找LLb0_K和LLb0_Kp
%     gap_K = 100; % meV
%     gap_Kp = 100; % meV
    LLb0_K_index = 0;
    LLb0_Kp_index = 0;
    for i = 1:dims_b
        slope_K_now = abs(LL_K_b(start_index, i) - LL_K_b(end_index, i)) / (end_index - start_index + 1);
%         gap_K_now = sum(abs(LL_K_b(start_index:end_index,i) - LL_K_b_0)) / (end_index - start_index + 1);
        slope_Kp_now = abs(LL_Kp_b(start_index, i) - LL_Kp_b(end_index, i)) / (end_index - start_index + 1);
%         gap_Kp_now = sum(abs(LL_Kp_b(start_index:end_index,i) - LL_Kp_b_0)) / (end_index - start_index + 1);

        if slope_K_now < slope_eps
            gap_K_now = sum(abs(LL_K_b(start_index:end_index,i) - LL_K_b_0)) / (end_index - start_index + 1);
            if gap_K_now < gap_error
                LLb0_K_index = i;
            end
        end
        
        if slope_Kp_now < slope_eps
            gap_Kp_now = sum(abs(LL_Kp_b(start_index:end_index,i) - LL_Kp_b_0)) / (end_index - start_index + 1);
            
            if gap_Kp_now < gap_error
                LLb0_Kp_index = i;
            end
        end
        
        if ~(LLb0_Kp_index == 0) && ~(LLb0_K_index == 0)
            break
        end
%         if gap_K_now < gap_K 
%             gap_K = gap_K_now;
%             LLb0_K_index = i;
%         end
% 
%         if gap_Kp_now < gap_Kp
%             gap_Kp = gap_Kp_now;
%             LLb0_Kp_index = i;
%         end
    end

    % 寻找LLm0_K和LLm0_Kp
    gap_K = 100; % meV
    gap_Kp = 100; % meV
    LLm0_K_index = 0;
    LLm0_Kp_index = 0;
    for i = 1:dims_m
        gap_K_now = sum(abs(LL_K_m(start_index:end_index,i) - LL_K_m_0)) / (end_index - start_index + 1);
        gap_Kp_now = sum(abs(LL_Kp_m(start_index:end_index,i) - LL_Kp_m_0)) / (end_index - start_index + 1);
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