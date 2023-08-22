function [eigvec_A_component, eigvec_B_component] = get_LL_components_each_sublattice_two_bands(eigvec, LL_index_cutoff, valley)
    % 得到sublattice上朗道能级对应的分量
    % two bands model : (B2, A1) @ K valley & (A1, B2) @ Kp valley
    eigvec_A_component = zeros(LL_index_cutoff + 1, 1);
    eigvec_B_component = zeros(LL_index_cutoff + 1, 1);
    
    if valley == 1 % for K valley
        eigvec_B_component(1) = eigvec(1); % |B2, 0>的分量
        eigvec_B_component(2) = eigvec(2); % |B2, 1>的分量
        for n = 2:LL_index_cutoff
            eigvec_B_component(n + 1) = eigvec(2*(n - 1) + 1); % |B2, n>的分量
            eigvec_A_component(n - 1) = eigvec(2*n); % |A1, n - 2>的分量
        end
   
    else % for Kp valley
        eigvec_A_component(1) = eigvec(1); % |A1, 0>的分量
        eigvec_A_component(2) = eigvec(2); % |A1, 1>的分量
        for n = 2:LL_index_cutoff
            eigvec_A_component(n + 1) = eigvec(2*(n - 1) + 1); % |A1, n>的分量
            eigvec_B_component(n - 1) = eigvec(2*n); % |B2, n - 2>的分量
        end
        
    end
end