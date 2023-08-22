function S_tensor = get_S_tensor_T2(exchange_integrals_list, transform_mat_cell, dist_index_array, indice_cell, LL_index_max)
    S_tensor = zeros(3,3,3,3); % (a,a',b,b')
    for a = 1:3
        for b = 1:3
            for c = 1:3
                for d = 1:3
                    S_tensor_element = get_S_tensor_element(exchange_integrals_list, transform_mat_cell, dist_index_array, a, b, c, d, indice_cell, LL_index_max);
                    S_tensor(a,d,c,b) = S_tensor_element;
                end
            end
        end
    end
end

function S_tensor_element = get_S_tensor_element(exchange_integrals_list, transform_mat_cell, dist_index_array, a, b, c, d, indice_cell, LL_index_max)
    % a : 产生 b : 湮灭 c : 产生  d : 湮灭
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
                    for alpha = 1:6 %% 这里有问题，应该是6，因为trilayer ABA 的num_sublattice是6
                        for beta = 1:6
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


%% test
% for a = 1:3
%     for d = 1:3
%         for c = 1:3
%             for b = 1:3
%                 if abs(S_tensor(a, d, c, b) - conj(S_tensor(b, c, d, a))) >= 1e-20
%                     disp("不厄密")
%                     S_tensor(a, d, c, b)
%                     conj(S_tensor(b, c, d, a))
%                     S_tensor(a, d, c, b) - conj(S_tensor(b, c, d, a))
%                 end
%             end
%         end
%     end
% end