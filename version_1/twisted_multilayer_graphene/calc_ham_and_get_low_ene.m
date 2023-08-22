function [eig_vals, eig_vals_low_ene_temp, eig_vecs_low_ene_temp, num_low_band] = calc_ham_and_get_low_ene(H_ham, eig_vals, low_ene_bound, dims, i)
    %% 对(akx, aky)的哈密顿量进行对角化
    % valley K
    [eig_vecs_temp, eigvals] = eig(H_ham);
    eig_vals_temp = diag(eigvals);
    eig_vals(i, :) = eig_vals_temp;
    % 按照绝对值重新进行排序, 可以减少计算次数
    [~, new_order] = sort(abs(eig_vals_temp));
    eig_vals_temp = eig_vals_temp(new_order);
    eig_vecs_temp = eig_vecs_temp(:, new_order);
    
    %% 得到低能能带数num_K_low_band
    num_low_band = 0;
    for jj = 1:dims
        % if abs(eig_vals_K_temp(jj)) < low_ene_bound + 0.01  % 适用于low_ene_bound会自动调整的情况
        if abs(eig_vals_temp(jj)) < low_ene_bound  % 适用于low_ene_bound固定的情况
            num_low_band = num_low_band + 1;
        else
            break
        end
    end
    
    if ~(num_low_band == 0)
        eig_vals_low_ene_temp = eig_vals_temp(1:num_low_band);
        eig_vecs_low_ene_temp = eig_vecs_temp(:, 1:num_low_band);
        
        [~, new_order] = sort(eig_vals_low_ene_temp);
        eig_vals_low_ene_temp = eig_vals_low_ene_temp(new_order); % 不按绝对值从小到大排
        eig_vecs_low_ene_temp = eig_vecs_low_ene_temp(:, new_order);
    else
        eig_vals_low_ene_temp = [];
        eig_vecs_low_ene_temp = [];
    end
end