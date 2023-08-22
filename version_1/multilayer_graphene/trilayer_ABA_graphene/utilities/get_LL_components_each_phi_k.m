function eigvec_phi_k_component_array = get_LL_components_each_phi_k(eigvec, LL_index_cutoff, valley, dims_m, dims_b)
    % 因为我们在ABA trilayer graphene中用的基组为【phi1, phi2, phi3, phi4, phi5, phi6】
    % six bands model : phi_k @ K valley & phi_k @ Kp valley
    % monolayer-like brach : |phi1> = 1 / sqrt(2) * 【|A1> - |A3>】 ; |phi2> = 1 / sqrt(2) * 【|B1> - |B3>】
    % bilayer-like brach : |phi3> = 1 / sqrt(2) * 【|A1> + |A3>】 ; |phi4> = 1 / sqrt(2) * 【|B1> + |B3>】 ; |phi5> = |A2> ; |phi6> = |B2>
    % 得到 phi_k 上朗道能级对应的分量
    
    % 对于单层，维数为dims_m = 2 * LL_index_cutoff + 1 ; 对于双层，维数为dims_b = 4 * LL_index_cutoff
    % valley : + 1 for K valley ; -1 for Kp valley
    
    % in K valley :
    % 对于monolayer graphene, the Landau levels in the order : [|0,phi2> // |1,phi2>, |0,phi1> // |2,phi2>, |1,phi1> // ... ...]
    % 对于bilayer graphene, the Landau levels in the order : [|0,phi6>, |1,phi6>, |0,phi4>, |0,phi5> // |2,phi6>, |0,phi3>, |1,phi4>, |1,phi5> // |3,phi6>, |1,phi3>, |2,phi4>, |2,phi5> // ... ...]
    
    % in Kp valley : 
    % 对于monolayer graphene, the Landau levels in the order : [|0,phi1> // |1,phi1>, |0,phi2> // |2,phi1>, |1,phi2> // ... ...]
    % 对于bilayer graphene, the Landau levels in the order : [|0,phi3>, |1,phi3>, |0,phi5>, |0,phi4> // |2,phi3>, |0,phi6>, |1,phi5>, |1,phi4> // |3,phi3>, |1,phi6>, |2,phi5>, |2,phi4> // ... ...]
    
    eigvec_phi_k_component_array = zeros(LL_index_cutoff + 1, 6); % 第一个维度是第i个矢量，第二个维度存放的是每个phi_k下第n个LL上的权重，第三个维度对应phi_k
    
    if valley == 1 % for K valley
        % monolayer-like branch
        eigvec_phi_k_component_array(1, 2) = eigvec(1); % |0,phi2>上的分量
        for n = 1:LL_index_cutoff
            eigvec_phi_k_component_array(n + 1, 2) = eigvec(2*n); % |n, phi2>上的分量
            eigvec_phi_k_component_array(n, 1) = eigvec(2*n + 1); % |n-1, phi1>上的分量
        end
        
        % bilayer-like branch
        eigvec_phi_k_component_array(1, 6) = eigvec(dims_m + 1); % |0,phi6>上的分量
        eigvec_phi_k_component_array(2, 6) = eigvec(dims_m + 2); % |1,phi6>上的分量
        eigvec_phi_k_component_array(1, 4) = eigvec(dims_m + 3); % |0,phi4>上的分量
        eigvec_phi_k_component_array(1, 5) = eigvec(dims_m + 4); % |0,phi5>上的分量
        for n = 2:LL_index_cutoff
            eigvec_phi_k_component_array(n + 1, 6) = eigvec(dims_m + 4 * n - 3); % |n, phi6>上的分量
            eigvec_phi_k_component_array(n - 1, 3) = eigvec(dims_m + 4 * n - 2); % |n-2, phi3>上的分量
            eigvec_phi_k_component_array(n, 4) = eigvec(dims_m + 4 * n - 1); % |n-1, phi4>上的分量
            eigvec_phi_k_component_array(n, 5) = eigvec(dims_m + 4 * n); % |n-1, phi5>上的分量
        end
        
    else % for Kp valley
        % monolayer-like branch
        eigvec_phi_k_component_array(1, 1) = eigvec(1); % |0,phi1>上的分量
        for n = 1:LL_index_cutoff
            eigvec_phi_k_component_array(n + 1, 1) = eigvec(2*n); % |n, phi1>上的分量
            eigvec_phi_k_component_array(n, 2) = eigvec(2*n + 1); % |n-1, phi2>上的分量
        end
        
        % bilayer-like branch
        eigvec_phi_k_component_array(1, 3) = eigvec(dims_m + 1); % |0,phi3>上的分量
        eigvec_phi_k_component_array(2, 3) = eigvec(dims_m + 2); % |1,phi3>上的分量
        eigvec_phi_k_component_array(1, 5) = eigvec(dims_m + 3); % |0,phi5>上的分量
        eigvec_phi_k_component_array(1, 4) = eigvec(dims_m + 4); % |0,phi4>上的分量
        for n = 2:LL_index_cutoff
            eigvec_phi_k_component_array(n + 1, 3) = eigvec(dims_m + 4 * n - 3); % |n, phi3>上的分量
            eigvec_phi_k_component_array(n - 1, 6) = eigvec(dims_m + 4 * n - 2); % |n-2, phi6>上的分量
            eigvec_phi_k_component_array(n, 5) = eigvec(dims_m + 4 * n - 1); % |n-1, phi5>上的分量
            eigvec_phi_k_component_array(n, 4) = eigvec(dims_m + 4 * n); % |n-1, phi4>上的分量
        end
        
    end
end