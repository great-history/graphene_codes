function [B_cross, E_cross] = trilayer_LLs_find_cross_points(LL_K_b, LL_Kp_b, LL_K_m, LL_Kp_m, LLb2_K_index, LLb2_Kp_index, LLm0_K_index, LLm0_Kp_index, B_fields, B_steps)
    % % 寻找LLb2_K, LLb2_Kp与LLm0_K, LLm0_Kp的交叉点位置,共有四个
    B_cross = zeros(4,1);
    E_cross = zeros(4,1);
    
    [B_cross(1), cross_index] = helper_find_cross_points(LL_K_b, LL_K_m, LLb2_K_index, LLm0_K_index, B_fields, B_steps);
    if cross_index == -1  % 交叉点落在磁场范围之外
        E_cross(1) = LL_K_m(end, LLm0_K_index);
    else
        E_cross(1) = LL_K_m(cross_index, LLm0_K_index);
    end
    
    [B_cross(2), cross_index] = helper_find_cross_points(LL_Kp_b, LL_K_m, LLb2_Kp_index, LLm0_K_index, B_fields, B_steps);
    if cross_index == -1  % 交叉点落在磁场范围之外
        E_cross(2) = LL_K_m(end, LLm0_K_index);
    else
        E_cross(2) = LL_K_m(cross_index, LLm0_K_index);
    end
    
    [B_cross(3), cross_index] = helper_find_cross_points(LL_K_b, LL_Kp_m, LLb2_K_index, LLm0_Kp_index, B_fields, B_steps);
    if cross_index == -1  % 交叉点落在磁场范围之外
        E_cross(3) = LL_Kp_m(end, LLm0_Kp_index);
    else
        E_cross(3) = LL_Kp_m(cross_index, LLm0_Kp_index);
    end
    
    [B_cross(4), cross_index] = helper_find_cross_points(LL_Kp_b, LL_Kp_m, LLb2_Kp_index, LLm0_Kp_index, B_fields, B_steps);
    if cross_index == -1  % 交叉点落在磁场范围之外
        E_cross(4) = LL_Kp_m(end, LLm0_Kp_index);
    else
        E_cross(4) = LL_Kp_m(cross_index, LLm0_Kp_index);
    end
end