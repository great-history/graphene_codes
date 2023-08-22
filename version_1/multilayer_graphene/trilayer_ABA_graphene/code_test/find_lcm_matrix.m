%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 找到二维矩阵中的所有极小值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [peak_value_list, peak_location_list] = find_lcm_matrix(mat)
    if length(size(mat)) > 2
        peak_value_list = [];
        peak_location_list = [];
        disp("维度问题")
        return
    end
    
    if length(size(mat)) < 1
        peak_value_list = [];
        peak_location_list = [];
        disp("维度问题")
        return
    end
    
    % 搜索极值点，边界不考虑
    peak_value_list = [];
    peak_location_list = [];
    for i = 2:(size(mat, 1) - 1)
        temp_vec = squeeze(mat(i, :));
        [temp_pk_list,temp_lc_list] = findpeaks(temp_vec);
        for j = 1:length(temp_lc_list)
            temp_lc = temp_lc_list(j);
            temp_pk = temp_pk_list(j);
            flag_true = (temp_pk > mat(i + 1, temp_lc)) & (temp_pk > mat(i - 1, temp_lc));
            if flag_true % 找到一个local maximum
                peak_value_list = [peak_value_list; temp_pk];
                peak_location_list = [peak_location_list; [i, temp_lc]];
            end
        end
    end
end