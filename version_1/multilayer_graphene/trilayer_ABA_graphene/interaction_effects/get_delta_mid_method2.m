function Delta_mid = get_delta_mid_method2(density_matrix_down_temp, density_matrix_up_temp, transform_mat_cell, LL_index_max)
    % % 第二种方式 Delta_mid = 0.0; for a = 1:3
    Delta_mid = 0.0;
    for n = 0:LL_index_max
        W_left = conj(transform_mat_cell{n + 1});
        W_right = transpose(transform_mat_cell{n + 1});
        
        density_mat_down = W_left * density_matrix_down_temp * W_right;
        density_mat_up = W_left * density_matrix_up_temp * W_right;
        
        Delta_mid = Delta_mid + density_mat_down(3,3) + density_mat_down(4,4);
        Delta_mid = Delta_mid + density_mat_up(3,3) + density_mat_up(4,4);
    end
end

% 
% function Delta_mid = get_delta_mid_method2(density_matrix_down_temp, density_matrix_up_temp, T2_eigvecs_cell, LL_index_max)
%     % % 第二种方式 Delta_mid = 0.0; for a = 1:3
%     Delta_mid = 0.0;
%     for a = 1:3
%         for b = 1:3
%             
%             density = 0.0;
%             for n = 0:LL_index_max
%                density = density + conj(T2_eigvecs_cell{n + 1}(5, a)) * T2_eigvecs_cell{n + 1}(5, b);
%                density = density + conj(T2_eigvecs_cell{n + 1}(6, a)) * T2_eigvecs_cell{n + 1}(6, b);
%             end
%             Delta_mid = Delta_mid + (density_matrix_down_temp(a, b) + density_matrix_up_temp(a, b)) * density;
%             
%         end
%     end
% end