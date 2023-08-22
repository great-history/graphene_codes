function eigvec_m_phi_k_component_array = get_LL_m_components_each_phi_k(eigvec_m, LL_index_cutoff, valley)
    % 因为我们在ABA trilayer graphene中用的基组为【phi1, phi2, phi3, phi4, phi5, phi6】
    % six bands model : phi_k @ K valley & phi_k @ Kp valley
    % monolayer-like brach : |phi1> = 1 / sqrt(2) * 【|A1> - |A3>】 ; |phi2> = 1 / sqrt(2) * 【|B1> - |B3>】
    % 得到 phi_k 上朗道能级对应的分量
    
    % 对于单层，维数为dims_m = 2 * LL_index_cutoff + 1 ; 对于双层，维数为dims_b = 4 * LL_index_cutoff
    % valley : + 1 for K valley ; -1 for Kp valley
    
    % in K valley :
    % 对于monolayer graphene, the Landau levels in the order : [|0,phi2> // |1,phi2>, |0,phi1> // |2,phi2>, |1,phi1> // ... ...]
    % in Kp valley : 
    % 对于monolayer graphene, the Landau levels in the order : [|0,phi1> // |1,phi1>, |0,phi2> // |2,phi1>, |1,phi2> // ... ...]
    
    eigvec_m_phi_k_component_array = zeros(LL_index_cutoff + 1, 2); % 第一个维度是第i个矢量，第二个维度存放的是phi_k
    
    if valley == 1 % for K valley
        % monolayer-like branch
        eigvec_m_phi_k_component_array(1, 2) = eigvec_m(1); % |0,phi2>上的分量
        for n = 1:LL_index_cutoff
            eigvec_m_phi_k_component_array(n + 1, 2) = eigvec_m(2*n); % |n, phi2>上的分量
            eigvec_m_phi_k_component_array(n, 1) = eigvec_m(2*n + 1); % |n-1, phi1>上的分量
        end
        
    else % for Kp valley
        % monolayer-like branch
        eigvec_m_phi_k_component_array(1, 1) = eigvec_m(1); % |0,phi1>上的分量
        for n = 1:LL_index_cutoff
            eigvec_m_phi_k_component_array(n + 1, 1) = eigvec_m(2*n); % |n, phi1>上的分量
            eigvec_m_phi_k_component_array(n, 2) = eigvec_m(2*n + 1); % |n-1, phi2>上的分量
        end
        
    end
end