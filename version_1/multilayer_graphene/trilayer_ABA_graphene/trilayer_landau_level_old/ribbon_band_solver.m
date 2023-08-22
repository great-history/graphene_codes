function [eig_enes, guide_centers, vel_array, vel_dirs] = ribbon_band_solver(H0, H1, H_1, y_coords, aks, eps, order_flag)
    dims = size(H0,1);
    k_points = length(aks);
    eig_enes = zeros(dims, k_points);

    guide_centers = zeros(dims, k_points);
    vel_dirs = zeros(dims, k_points);
    vel_array = zeros(dims, k_points);

    tic
    
    idx = 1;
    Im = 1j;
    Hk = H0 + H1 * exp(Im * aks(idx)) + H_1 * exp(-Im * aks(idx));
    [eigvecs_last, eigvals_last] = eig(Hk);
    eigvals_last = diag(eigvals_last);
    [~, degeneracy, eig_num, first_idxs, last_idxs] = helper_get_degeneracy(eigvals_last, dims);

    Vel_k = Im * ( H1 * exp(Im * aks(idx)) - H_1 * exp(-Im * aks(idx)) );
    [eigvecs_last, vel_dirs(:, idx), vel_array(:, idx)] = get_bloch_velocity(Vel_k, eigvecs_last, degeneracy, eig_num, first_idxs, last_idxs, dims);

    eig_enes(:,idx) = eigvals_last;

    for num = 1:dims
        guide_centers(num, idx) = dot(y_coords, abs(eigvecs_last(:, num)).^2);
    end
    
    if order_flag == false
        for idx = 2:k_points
            Hk = H0 + H1 * exp(Im * aks(idx)) + H_1 * exp(-Im * aks(idx));
            [eigvecs_now, eigvals_now] = eig(Hk);
            eigvals_now = diag(eigvals_now);
            [~, degeneracy, eig_num, first_idxs, last_idxs] = helper_get_degeneracy(eigvals_now, dims); % 得到简并度

            % bloch velocity 
            Vel_k = Im * (H1 * exp(Im * aks(idx)) - H_1 * exp(-Im * aks(idx)));
            [eigvecs_now, vel_dirs(:, idx), vel_array(:, idx)] = get_bloch_velocity(Vel_k, eigvecs_now, degeneracy, eig_num, first_idxs, last_idxs, dims);

            % 计算guide center
            for num = 1:dims
                guide_centers(num, idx) = dot(y_coords, abs(eigvecs_now(:, num)).^2);
            end

            eig_enes(:,idx) = eigvals_now;
        end
    else
        for idx = 2:k_points
            Hk = H0 + H1 * exp(Im * aks(idx)) + H_1 * exp(-Im * aks(idx));
            [eigvecs_now, eigvals_now] = eig(Hk);
            eigvals_now = diag(eigvals_now);
            [~, degeneracy, eig_num, first_idxs, last_idxs] = helper_get_degeneracy(eigvals_now, dims); % 得到简并度

            % bloch velocity 
            Vel_k = Im * (H1 * exp(Im * aks(idx)) - H_1 * exp(-Im * aks(idx)));
            [eigvecs_now, vel_dirs(:, idx), vel_array(:, idx)] = get_bloch_velocity(Vel_k, eigvecs_now, degeneracy, eig_num, first_idxs, last_idxs, dims);

            delta_k = aks(idx) - aks(idx - 1);
            eigvals_predict = eigvals_now - delta_k * vel_array(:, idx);
            new_order = helper_get_new_order(eigvecs_last, eigvals_last, eigvecs_now, eigvals_predict, degeneracy, first_idxs, last_idxs, dims, eig_num, eps);
                
            % 交换位置
            eigvecs_now(:,:) = eigvecs_now(:,new_order);
            eigvals_now(:) = eigvals_now(new_order);
            vel_dirs(:, idx) = vel_dirs(new_order, idx);
            vel_array(:, idx) = vel_array(new_order, idx);

            % 计算guide center
            for num = 1:dims
                guide_centers(num, idx) = dot(y_coords, abs(eigvecs_now(:, num)).^2);
            end

            eig_enes(:,idx) = eigvals_now;

            eigvecs_last = eigvecs_now;
            eigvals_last = eigvals_now;
        end
    end
    
    toc
end