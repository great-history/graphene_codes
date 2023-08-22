function [X_mat_up, X_mat_down] = get_fock_mat(density_matrix_up_temp_LL, density_matrix_down_temp_LL, exchange_integrals_list, ...
                                               transform_mat_cell, indice_cell, n1, m2, dist_index_array, LL_index_max)
    transform_mat_left = (transform_mat_cell{n1 + 1})';
    transform_mat_right = transform_mat_cell{m2 + 1};

    for m1 = 0:LL_index_max
        n2 = n1 + m1 - m2;
        if n2 < 0 || n2 > LL_index_max
            continue
        else
            [X_mat_up, X_mat_down] = get_fock_term_LL_basis(density_matrix_up_temp_LL, density_matrix_down_temp_LL, exchange_integrals_list, ...
                                                            transform_mat_left, transform_mat_right, indice_cell, n1, n2, m1, m2, dist_index_array);
            
        end
    end

end

function [X_mat_up, X_mat_down] = get_fock_term_LL_basis(density_matrix_up_temp_LL, density_matrix_down_temp_LL, exchange_integrals_list, ...
                                                         transform_mat_left, transform_mat_right, indice_cell, n1, n2, m1, m2, dist_index_array)
    a = abs(n1 - n2);
    b_prime = min(n1, n2);
    c_prime = min(m1, m2);
    b = max(b_prime, c_prime);
    c = min(b_prime, c_prime);
    ex_index = indice_cell{a + 1}{b + 1}{c + 1};
    
    X_mat_up = zeros(6);                                                 
    X_mat_down = zeros(6); 
    index_left = 6 * m1;
    index_right = 6 * n2;
    for alpha = 1:6
        for beta = 1:6
            dist_index = dist_index_array(alpha, beta);
            X_mat_up(alpha, beta) = exchange_integrals_list(ex_index, dist_index) * density_matrix_up_temp_LL(index_left + beta, index_right + alpha);
            X_mat_down(alpha, beta) = exchange_integrals_list(ex_index, dist_index) * density_matrix_down_temp_LL(index_left + beta, index_right + alpha);
        end
    end
    
    X_mat_up = transform_mat_left * X_mat_up * transform_mat_right;
    X_mat_down = transform_mat_left * X_mat_down * transform_mat_right;
    
end