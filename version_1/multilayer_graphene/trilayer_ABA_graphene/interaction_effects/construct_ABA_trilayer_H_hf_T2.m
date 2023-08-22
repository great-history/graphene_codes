function [H_hf_up, H_hf_down] = construct_ABA_trilayer_H_hf_T2(density_matrix_up_temp, density_matrix_down_temp, S_tensor, R_mat, ...
                                                               T2_eigvals_list, E_zeeman, E_H, Delta_mid, E_exchange)
    % S_tensor是用来计算fock term，而R_mat是用来计算hartree term
    % energy scale : T2_eigvals_list // E_zeeman // E_H * Delta_mid / 2 // E_exchange
    
    %% 构造在T2 basis下的Hartree-Fock Hamiltonian
    % eigenvalues
    H_hf_down = diag(T2_eigvals_list);
    H_hf_up = diag(T2_eigvals_list);
    
    % Zeeman term
    H_hf_down = H_hf_down - E_zeeman * diag([1,1,1]);
    H_hf_up = H_hf_up + E_zeeman * diag([1,1,1]);
    
    % Hartree-term
    H_hf_down = H_hf_down + E_H / 2 * Delta_mid * R_mat;
    H_hf_up = H_hf_up + E_H / 2 * Delta_mid * R_mat;
    
    % fock term
    fock_mat_down = zeros(3);
    fock_mat_up = zeros(3);
    for a = 1:3
        for b = 1:3
            for c = 1:3
                for d = 1:3
                    fock_mat_down(a,b) = fock_mat_down(a,b) + density_matrix_down_temp(c,d) * S_tensor(a,d,c,b);
                    fock_mat_up(a,b) = fock_mat_up(a,b) + density_matrix_up_temp(c,d) * S_tensor(a,d,c,b);
                end
            end
        end
    end
    fock_mat_down = - E_exchange * fock_mat_down;
    fock_mat_up = - E_exchange * fock_mat_up;
    
    H_hf_down = H_hf_down + fock_mat_down;
    H_hf_up = H_hf_up + fock_mat_up;
end

% function [H_hf_up, H_hf_down] = construct_ABA_trilayer_H_hf_T2(density_matrix_up_temp_LL, density_matrix_down_temp_LL, exchange_integrals_list, ...
%                                                                dist_index_array, indice_cell, transform_mat_cell, ...
%                                                                T2_eigvals_list, E_zeeman, E_H, Delta_mid, E_exchange, LL_index_max)
%     % energy scale : T2_eigvals_list // E_zeeman // 
%     %% 构造在T2 basis下的Hartree-Fock Hamiltonian
%     % eigenvalues
%     H_hf_down = diag(T2_eigvals_list);
%     H_hf_up = diag(T2_eigvals_list);
%     
%     % Zeeman term
%     H_hf_down = H_hf_down - E_zeeman * diag([1,1,1]);
%     H_hf_up = H_hf_up + E_zeeman * diag([1,1,1]);
% 
%     % Hartree term
%     Ham_hartree = get_hartree_mat(Delta_mid, E_H, transform_mat_cell, LL_index_max);
%     H_hf_down = H_hf_down + Ham_hartree;
%     H_hf_up = H_hf_up + Ham_hartree;
%     
%     %% 如果没有加入Fock term, 密度矩阵是可以很快收敛的，说明问题出在Fock term
%     % Fock term of H_hf_down & H_hf_up
%     % 开始计算矩阵元<a|U_ex|b>
%     for n1 = 0:LL_index_max
%         for m2 = 0:LL_index_max
%             [X_mat_up, X_mat_down] = get_fock_mat(density_matrix_up_temp_LL, density_matrix_down_temp_LL, exchange_integrals_list, ...
%                                                   transform_mat_cell, indice_cell, n1, m2, dist_index_array, LL_index_max);
%             H_hf_up = H_hf_up - E_exchange * X_mat_up;
%             H_hf_down = H_hf_down - E_exchange * X_mat_down;
%         end
%     end
%     
%     %     for n1 = 0:LL_index_max
%     %         for m2 = 0:LL_index_max
%     %             transform_mat_left = (transform_mat_cell{n1 + 1})';
%     %             transform_mat_right = transform_mat_cell{m2 + 1};
%     %
%     %             X_mat_up = zeros(6);
%     %             X_mat_down = zeros(6);
%     %
%     %             for m1 = 0:LL_index_max
%     %                 n2 = n1 + m1 - m2;
%     %                 if n2 < 0 || n2 > LL_index_max
%     %                     continue
%     %                 else
%     %                     a = abs(n1 - n2);
%     %                     b_prime = min(n1, n2);
%     %                     c_prime = min(m1, m2);
%     %                     b = max(b_prime, c_prime);
%     %                     c = min(b_prime, c_prime);
%     %                     ex_index = indice_cell{a + 1}{b + 1}{c + 1};
%     %
%     %                     index_left = 6 * m1;
%     %                     index_right = 6 * n2;
%     %
%     %                     for alpha = 1:6
%     %                         for beta = 1:6
%     %                             dist_index = dist_index_array(alpha, beta);
%     %                             X_mat_down(alpha, beta) = X_mat_down(alpha, beta) + ...
%     %                                                       exchange_integrals_list(ex_index, dist_index) * density_matrix_down_temp_LL(index_left + beta, index_right + alpha);
%     %                             X_mat_up(alpha, beta) = X_mat_up(alpha, beta) + ...
%     %                                                     exchange_integrals_list(ex_index, dist_index) * density_matrix_up_temp_LL(index_left + beta, index_right + alpha);
%     %                         end
%     %                     end
%     %                 end
%     %             end
%     %
%     %             X_mat_up = - E_exchange * X_mat_up;
%     %             X_mat_down = - E_exchange * X_mat_down;
%     %             UX_mat_up = transform_mat_left * X_mat_up * transform_mat_right;
%     %             UX_mat_down = transform_mat_left * X_mat_down * transform_mat_right;
%     %
%     %             H_hf_up = H_hf_up + UX_mat_up;
%     %             H_hf_down = H_hf_down + UX_mat_down;
%     %             % H_hf_up = H_hf_up + transform_mat_left * X_mat_up * transform_mat_right;
%     %             % H_hf_down = H_hf_down + transform_mat_left * X_mat_down * transform_mat_right;
%     %         end
%     %     end
%     
%     % H_hf_down
%     % helper_check_hermite(H_hf_down, 1e-16)
% end