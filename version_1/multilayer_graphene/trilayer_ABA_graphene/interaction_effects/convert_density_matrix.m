function [density_matrix_up_temp_LL, density_matrix_down_temp_LL] = convert_density_matrix(density_matrix_up_temp, density_matrix_down_temp, ...
                                                                                           transform_mat_cell, num_LL, LL_index_max)
    % 将密度矩阵从eigenstate basis转换到LL basis
    density_matrix_down_temp_LL = zeros(6*num_LL); % 只有一个谷,并且是spin down，所以是6*num_LL 【|A1,n>，|B1,n>，|A2,n>，|B2,n>，|A3,n>，|B3,n>】
    density_matrix_up_temp_LL = zeros(6*num_LL); % 只有一个谷,并且是spin down，所以是6*num_LL 【|A1,n>，|B1,n>，|A2,n>，|B2,n>，|A3,n>，|B3,n>
    for n = 0:LL_index_max
        index_left = 6 * n;
        transform_mat_left = conj(transform_mat_cell{n + 1});

        for m = 0:LL_index_max
            index_right = 6 * m;
            transform_mat_right = transpose(transform_mat_cell{m + 1});
            
            for alpha = 1:6
                for beta = 1:6
                    val_up = 0.0;
                    val_down = 0.0;
                    for a = 1:3
                        for b = 1:3
                            val_up = val_up + transform_mat_left(alpha, a) * density_matrix_up_temp(a, b) * transform_mat_right(b, beta);
                            val_down = val_down + transform_mat_left(alpha, a) * density_matrix_down_temp(a, b) * transform_mat_right(b, beta);
                        end
                    end
                    
                    
                    density_matrix_down_temp_LL(index_left + alpha, index_right + beta) = val_down;
                    density_matrix_up_temp_LL(index_left + alpha, index_right + beta) = val_up;
                end
            end
%             density_matrix_down_temp_LL((index_left + 1):(index_left + 6), (index_right + 1):(index_right + 6)) = ...
%                                                     transform_mat_left * density_matrix_down_temp * transform_mat_right;
%             density_matrix_up_temp_LL((index_left + 1):(index_left + 6), (index_right + 1):(index_right + 6)) = ...
%                                                     transform_mat_left * density_matrix_up_temp * transform_mat_right;
        end
    end
    % helper_check_hermite(density_matrix_down_temp_LL, 1e-16)  1
    % helper_check_hermite(density_matrix_down_temp_LL, 1e-20)  0
    
end