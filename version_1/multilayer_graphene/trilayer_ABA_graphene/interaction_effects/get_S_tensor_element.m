function S_tensor_element = get_S_tensor_element(exchange_integrals_list, transform_mat_cell, dist_index_array, a, b, c, d, indice_cell, LL_index_max, num_sublattice)
    % a : 产生        b : 湮灭        c : 产生        d : 湮灭
    S_tensor_element = 0.0;
    for n1 = 0:LL_index_max
        for m2 = 0:LL_index_max
            Wa = (transform_mat_cell{n1 + 1})';
            Wb = (transform_mat_cell{m2 + 1});
            for m1 = 0:LL_index_max
                for n2 = 0:LL_index_max
                    if ~(n1 - n2 + m1 - m2 == 0)
                        continue
                    end
                    
                    Wc = (transform_mat_cell{m1 + 1})';
                    Wd = (transform_mat_cell{n2 + 1});
                    
                    xx = abs(n1 - n2);
                    yy_prime = min(n1, n2);
                    zz_prime = min(m1, m2);
                    yy = max(yy_prime, zz_prime);
                    zz = min(yy_prime, zz_prime);
                    
                    ex_index = indice_cell{xx + 1}{yy + 1}{zz + 1};
                    
                    S_value = 0.0;
                    for alpha = 1:num_sublattice % 对于trilayer ABA而言，num_sublattice = 6，对于bilayer而言，num_sublattice = 4，对于monolayer而言，num_sublattice = 2
                        for beta = 1:num_sublattice
                            dist_index = dist_index_array(alpha, beta);
                            X_nm = exchange_integrals_list(ex_index, dist_index);
                            S_value = S_value + X_nm * Wa(a, alpha) * Wb(beta, b) * Wc(c, beta) * Wd(alpha, d);
                        end
                    end
                    
                    S_tensor_element = S_tensor_element + S_value;
                end
            end
        end
    end
end