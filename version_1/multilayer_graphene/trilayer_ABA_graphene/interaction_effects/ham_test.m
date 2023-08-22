H_hf_down = diag(T2_eigvals_list);
H_hf_down = H_hf_down - E_zeeman * diag([1,1,1]);

layer_potential_mat = diag([-1,-1,1,1,-1,-1]);
mat_temp = zeros(3);
for n = 0:LL_index_max
    transform_mat_left = (transform_mat_cell{n + 1})';
    transform_mat_right = transform_mat_cell{n + 1};
    mat_temp = mat_temp + (transform_mat_left * layer_potential_mat * transform_mat_right);
end

mat_temp = E_H / 2 * real(Delta_mid) * mat_temp;
H_hf_down = H_hf_down + mat_temp;
H_hf_up = H_hf_up + mat_temp;

for n1 = 0:LL_index_max
    for m2 = 0:LL_index_max
        transform_mat_left = (transform_mat_cell{n1 + 1})';
        transform_mat_right = transform_mat_cell{m2 + 1};

        X_mat_up = 0.025 * ones(6, 6);
        X_mat_down = 0.025 * ones(6, 6);
        
%         for m1 = 0:LL_index_max
%             n2 = n1 + m1 - m2;
%             if n2 < 0 || n2 > LL_index_max
%                 continue
%             else
%                 a = abs(n1 - n2);
%                 b_prime = min(n1, n2);
%                 c_prime = min(m1, m2);
%                 b = max(b_prime, c_prime);
%                 c = min(b_prime, c_prime);
%                 ex_index = indice_cell{a + 1}{b + 1}{c + 1};
%             end
% 
%             index_left = 6 * m1;
%             index_right = 6 * n2;
% 
%             for alpha = 1:6
%                 for beta = 1:6
%                     dist_index = dist_index_array(alpha, beta);
%                     X_mat_down(alpha, beta) = X_mat_down(alpha, beta) ...
%                         + exchange_integrals_list(ex_index, dist_index) * density_matrix_down_temp_LL(index_left + beta, index_right + alpha);
%                     X_mat_up(alpha, beta) = X_mat_up(alpha, beta) ...
%                         + exchange_integrals_list(ex_index, dist_index) * density_matrix_up_temp_LL(index_left + beta, index_right + alpha);
%                 end
%             end
% 
%         end
        
        UX_mat_up = - E_exchange * (transform_mat_left * X_mat_up * transform_mat_right);
        UX_mat_down = - E_exchange * (transform_mat_left * X_mat_down * transform_mat_right);

        % H_hf_up = H_hf_up + transform_mat_left * X_mat_up * transform_mat_right;
        % H_hf_down = H_hf_down + transform_mat_left * X_mat_down * transform_mat_right;
    end
end