function Ham_LL_valley = construct_monolayer_LL_for_ABA(x0, delta, gamma2, gamma5, Delta2, valley, LL_index_cutoff, dims)
    %% 函数说明
    % LL_index_cutoff = 10 说明最高LL的指标为10，即只考虑LL_index = 0,1,2,3,4,5,6,7,8,,9,10，维数为dims = 2 * LL_index_cutoff + 1
    % 对于trilayer ABA而言，通过基矢的变换，可以得到:
    % monolayer-like brach : |phi1> = 1 / sqrt(2) * 【|A1> - |A3>】 ; |phi2> = 1 / sqrt(2) * 【|B1> - |B3>】
    % bilayer-like brach : |phi3> = 1 / sqrt(2) * 【|A1> + |A3>】 ; |phi4> = 1 / sqrt(2) * 【|B1> + |B3>】 ; |phi5> = |A2> ; |phi6> = |B2>
    
    % valley : + 1 for K valley ; -1 for Kp valley
    % in K valley, the Landau levels in the order : [|0,phi2> // |1,phi2>, |0,phi1> // |2,phi2>, |1,phi1> // ... ...]
    % in Kp valley, the Landau levels in the order : [|0,phi1> // |1,phi1>, |0,phi2> // |2,phi1>, |1,phi2> // ... ...]
    
    %% 构造哈密顿量
    Ham_LL_valley = zeros(dims);
    
    if valley == 1
        Ham_LL_valley(1,1) = delta - gamma5 / 2 + Delta2;
        for n = 1:LL_index_cutoff
            H_n = get_diagonal_H_n(x0, delta, gamma2, gamma5, Delta2, valley, n);
            Ham_LL_valley(2*n:(2*n + 1), 2*n:(2*n + 1)) = H_n;
        end
    else
        Ham_LL_valley(1,1) = - gamma2 / 2 + Delta2;
        for n = 1:LL_index_cutoff
            H_n = get_diagonal_H_n(x0, delta, gamma2, gamma5, Delta2, valley, n);
            Ham_LL_valley(2*n:(2*n + 1), 2*n:(2*n + 1)) = H_n;
        end
    end
end


function H_n = get_diagonal_H_n(x0, delta, gamma2, gamma5, Delta2, valley, n)
    % diagonal elements
    if valley == 1
        H_n = diag([delta - gamma5 / 2 + Delta2, - gamma2 / 2 + Delta2]);
        H_n(1, 2) = - 1j * x0 * sqrt(2 * n);
        H_n(2, 1) = conj(H_n(1, 2));
    else
        H_n = diag([- gamma2 / 2 + Delta2, delta - gamma5 / 2 + Delta2]);
        H_n(1, 2) = 1j * x0 * sqrt(2 * n);
        H_n(2, 1) = conj(H_n(1, 2));
    end
end