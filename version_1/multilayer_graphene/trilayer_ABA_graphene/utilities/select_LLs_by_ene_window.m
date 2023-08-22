function [eigvec_select_array, eigval_select_list] = select_LLs_by_ene_window(eigvec_array, eigval_list, ene_ub, ene_lb, dim)
    % 找出位于energy window [ene_lb, ene_ub]范围内的所有本征态和相应的本征能量
    eigvec_select_array = [];
    eigval_select_list = [];
    for ii = 1:dim
        if eigval_list(ii) >= ene_lb && eigval_list(ii) <= ene_ub
            eigval_select_list = [eigval_select_list, eigval_list(ii)];
            eigvec_select_array = [eigvec_select_array, eigvec_array(:, ii)];
        end
    end
end