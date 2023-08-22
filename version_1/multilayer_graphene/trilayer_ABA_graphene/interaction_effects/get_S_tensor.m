function S_tensor = get_S_tensor(exchange_integrals_list, transform_mat_cell, dist_index_array, indice_cell, LL_index_max, dim_subspace)
    S_tensor = zeros(dim_subspace, dim_subspace, dim_subspace, dim_subspace); % (a,a',b,b')
    for a = 1:dim_subspace
        for b = 1:dim_subspace
            for c = 1:dim_subspace
                for d = 1:dim_subspace
                    
                    S_tensor_element = get_S_tensor_element(exchange_integrals_list, transform_mat_cell, dist_index_array, a, b, c, d, indice_cell, LL_index_max, 6);
                    S_tensor(a,d,c,b) = S_tensor_element;
                    
                end
            end
        end
    end
    
end