function [LL_K_m0_crossing_points_list, LL_Kp_m0_crossing_points_list, poly_LL_K_m0, poly_LL_Kp_m0, eigval_LL_K_b_list, eigval_LL_Kp_b_list] = ...
                trilayer_ABA_LLs_find_cross_points(eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, LL_K_m_0, LL_Kp_m_0, ene_width_disorder, B_steps, B_fields_list)
    B_start = B_fields_list(1);
    B_end = B_fields_list(end);
    % 确定LL_K_m0与(LL_K_b, LL_Kp_b)的交点， 确定LL_Kp_m0与(LL_K_b, LL_Kp_b)的交点
    % 首先把LL_K_b2, LL_K_b3, LL_K_b4, LL_Kp_b2, LL_Kp_b3, LL_Kp_b4 都找出来
    eigval_LL_K_b_list = zeros(4, B_steps); % 存放b1/b2/b3/b4
    eigval_LL_Kp_b_list = zeros(4, B_steps);
    
    % 找出 LL_b_K_1 // LL_b_K_2 // LL_b_K_3 // LL_b_K_4 和 LL_b_Kp_1 // LL_b_Kp_2 // LL_b_Kp_3 // LL_b_Kp_4
    for B_index = 1:B_steps
        ind1_K = find(eig_info_HK_b_select_cell{B_index, 4} == 1);
        for ii = 1:4 % 确定LL_K_b_1 // LL_K_b_2 // LL_K_b_3 // LL_K_b_4
            eigval_LL_K_b_list(ii, B_index) = eig_info_HK_b_select_cell{B_index, 2}(ind1_K + ii - 1);
        end

        ind1_Kp = find(eig_info_HKp_b_select_cell{B_index, 4} == 1);
        for ii = 1:4 % 确定LL_Kp_b_1 // LL_Kp_b_2 // LL_Kp_b_3 // LL_Kp_b_4
            eigval_LL_Kp_b_list(ii, B_index) = eig_info_HKp_b_select_cell{B_index, 2}(ind1_Kp + ii - 1);
        end
    end

    % 找出LL_K_m_0与各LL_K_b_i和LL_Kp_b_i交叉点
    LL_K_m0_crossing_points_list = zeros(2, 4); % 存放关于LL_K_m0交叉点的信息,第一行放left point,第二行放right point
    poly_LL_K_m0 = polyshape([B_start B_start B_end B_end], ...
                             [LL_K_m_0 + ene_width_disorder LL_K_m_0 - ene_width_disorder LL_K_m_0 - ene_width_disorder LL_K_m_0 + ene_width_disorder]);
    % figure
    % hold on
    % plot(poly_LL_K_m0)
    for ii = 1:4
        lineseg_K = [B_fields_list', eigval_LL_K_b_list(ii, :)'];
        [in_K, out_K] = intersect(poly_LL_K_m0, lineseg_K);
        if isempty(in_K)
            crossing_field_left = B_start;
            crossing_field_right = B_start;
        else
            if in_K(1,2) < LL_K_m_0
                crossing_field_left = in_K(1,1);
                crossing_field_right = in_K(end,1);
            else
                crossing_field_left = in_K(end,1);
                crossing_field_right = in_K(1,1);
            end
        end

        lineseg_Kp = [B_fields_list', eigval_LL_Kp_b_list(ii, :)'];
        [in_Kp, out_Kp] = intersect(poly_LL_K_m0, lineseg_Kp);
        if isempty(in_Kp)
            crossing_field_left = crossing_field_left + B_start;
            crossing_field_right = crossing_field_right + B_start;
        else
            if in_Kp(1,2) < LL_K_m_0
                crossing_field_left = crossing_field_left + in_Kp(1,1);
                crossing_field_right = crossing_field_right + in_Kp(end,1);
            else
                crossing_field_left = crossing_field_left + in_Kp(end,1);
                crossing_field_right = crossing_field_right + in_Kp(1,1);
            end
        end
        crossing_field_left = crossing_field_left / 2;
        crossing_field_right = crossing_field_right / 2;
        LL_K_m0_crossing_points_list(1, ii) = crossing_field_left;
        LL_K_m0_crossing_points_list(2, ii) = crossing_field_right;

        % plot(in_K(:,1),in_K(:,2),'b',out_K(:,1),out_K(:,2),'r')
        % plot(in_Kp(:,1),in_Kp(:,2),'b',out_Kp(:,1),out_Kp(:,2),'r--')
    end
    
    % 找出LL_Kp_m_0与各LL_K_b_i和LL_Kp_b_i交叉点
    LL_Kp_m0_crossing_points_list = zeros(2, 4); % 存放关于LL_Kp_m0交叉点的信息,第一行放left point,第二行放right point                         
    poly_LL_Kp_m0 = polyshape([B_start, B_start, B_end, B_end], ...
                              [LL_Kp_m_0 + ene_width_disorder, LL_Kp_m_0 - ene_width_disorder, LL_Kp_m_0 - ene_width_disorder, LL_Kp_m_0 + ene_width_disorder]);
    % hold on
    % plot(poly_LL_Kp_m0)
    for ii = 1:4
        lineseg_K = [B_fields_list', eigval_LL_K_b_list(ii, :)'];
        [in_K, out_K] = intersect(poly_LL_Kp_m0, lineseg_K);
        if isempty(in_K)
            crossing_field_left = B_start;
            crossing_field_right = B_start;
        else
            if in_K(1,2) < LL_Kp_m_0
                crossing_field_left = in_K(1,1);
                crossing_field_right = in_K(end,1);
            else
                crossing_field_left = in_K(end,1);
                crossing_field_right = in_K(1,1);
            end
        end

        lineseg_Kp = [B_fields_list', eigval_LL_Kp_b_list(ii, :)'];
        [in_Kp, out_Kp] = intersect(poly_LL_Kp_m0, lineseg_Kp);
        if isempty(in_Kp)
            crossing_field_left = crossing_field_left + B_start;
            crossing_field_right = crossing_field_right + B_start;
        else
            if in_Kp(1,2) < LL_Kp_m_0
                crossing_field_left = crossing_field_left + in_Kp(1,1);
                crossing_field_right = crossing_field_right + in_Kp(end,1);
            else
                crossing_field_left = crossing_field_left + in_Kp(end,1);
                crossing_field_right = crossing_field_right + in_Kp(1,1);
            end
        end
        crossing_field_left = crossing_field_left / 2;
        crossing_field_right = crossing_field_right / 2;
        LL_Kp_m0_crossing_points_list(1, ii) = crossing_field_left;
        LL_Kp_m0_crossing_points_list(2, ii) = crossing_field_right;

        % plot(in_K(:,1),in_K(:,2),'b',out_K(:,1),out_K(:,2),'r')
        % plot(in_Kp(:,1),in_Kp(:,2),'b',out_Kp(:,1),out_Kp(:,2),'r--')
    end
    
end