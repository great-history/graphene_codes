function Ham_LL_valley = construct_bilayer_LL_four_bands_for_ABA(x0, x3, x4, gamma1, delta, gamma2, gamma5, Delta2, valley, LL_index_cutoff, dims)
    %% 该函数是用来计算trilayer ABA graphene对应的bilayer-like brach，而不应该被用来计算Bilayer graphene（x4可能会不对）
    % 与bilayer_graphene那里不完全一样的是x4的处理，这里会比那里多一个负号()
    % LL_index_cutoff = 10 说明最高LL的指标为10，即只考虑LL_index = 0,1,2,3,4,5,6,7,8,,9,10，维数为dims = 2 * LL_index_cutoff
    % 对于trilayer ABA而言，通过基矢的变换，可以得到:
    % monolayer-like brach : |phi1> = 1 / sqrt(2) * 【|A1> - |A3>】 ; |phi2> = 1 / sqrt(2) * 【|B1> - |B3>】
    % bilayer-like brach : |phi3> = 1 / sqrt(2) * 【|A1> + |A3>】 ; |phi4> = 1 / sqrt(2) * 【|B1> + |B3>】 ; |phi5> = |A2> ; |phi6> = |B2>
    % 
    % valley : + 1 for K valley ; -1 for Kp valley
    % in K valley, the Landau levels in the order : [|0,phi6>, |1,phi6>, |0,phi4>, |0,phi5> // |2,phi6>, |0,phi3>, |1,phi4>, |1,phi5> // |3,phi6>, |1,phi3>, |2,phi4>, |2,phi5> // ... ...]
    % in Kp valley, the Landau levels in the order : [|0,phi3>, |1,phi3>, |0,phi5>, |0,phi4> // |2,phi3>, |0,phi6>, |1,phi5>, |1,phi4> // |3,phi3>, |1,phi6>, |2,phi5>, |2,phi4> // ... ...]
    
    %% 注意点：
    % 传入gamma1之前应将gamma1扩大sqrt(2)倍
    % 传入x3之前应将gamma3扩大sqrt(2)倍
    % 传入x4之前应将gamma4扩大sqrt(2)倍
    
    %%
    % 构建在谷valley处的LL能级哈密顿量
    Ham_LL_valley = zeros(dims);
    
    % 对 n = 1 block进行单独处理
    n = 1; 
    index_start = 4 * (n - 1) + 1;
    index_end = 4 * n;
    
    % H_n
    % H_n = get_diagonal_H_n(x0, x4, gamma1, delta, u, valley, n);
    H_n = get_diagonal_H_n(x0, x4, gamma1, delta, gamma2, gamma5, Delta2, valley, n);
    Ham_LL_valley(index_start:index_end, index_start:index_end) = H_n;
    
    % W1
    W_n = get_off_diagonal_W_n(x3, valley, n);
    Ham_LL_valley(index_start:index_end, (index_start + 8):(index_end + 8)) = W_n;
    Ham_LL_valley((index_start + 8):(index_end + 8), index_start:index_end) = W_n';

    % W2
    W_n = get_off_diagonal_W_n(x3, valley, n + 1);
    Ham_LL_valley(index_start:index_end, (index_start + 12):(index_end + 12)) = W_n;
    Ham_LL_valley((index_start + 12):(index_end + 12), index_start:index_end) = W_n';
    
    % 对 n = 2 ： (LL_index_cutoff - 3) block进行处理
    for n = 2:(LL_index_cutoff - 3)
        index_start = 4 * (n - 1) + 1;
        index_end = 4 * n;
        
        % H_n
        H_n = get_diagonal_H_n(x0, x4, gamma1, delta, gamma2, gamma5, Delta2, valley, n);
        Ham_LL_valley(index_start:index_end, index_start:index_end) = H_n;

        % Wn (n >= 3)
        W_n = get_off_diagonal_W_n(x3, valley, n + 1);
        Ham_LL_valley(index_start:index_end, (index_start + 12):(index_end + 12)) = W_n;
        Ham_LL_valley((index_start + 12):(index_end + 12), index_start:index_end) = W_n';
    end
    
    % 对 n = (LL_index_cutoff - 2) : LL_index_cutoff block 进行单独处理
    for n = (LL_index_cutoff - 2):LL_index_cutoff
        index_start = 4 * (n - 1) + 1;
        index_end = 4 * n;

        % H_n
        H_n = get_diagonal_H_n(x0, x4, gamma1, delta, gamma2, gamma5, Delta2, valley, n);
        Ham_LL_valley(index_start:index_end, index_start:index_end) = H_n;
        Ham_LL_valley(index_start:index_end, index_start:index_end) = H_n;
    end
    
end

function H_n = get_diagonal_H_n(x0, x4, gamma1, delta, gamma2, gamma5, Delta2, valley, n)
    % 与bilayer_graphene那里不完全一样的是x4的处理，这里会比那里多一个负号()
    if n == 1
        % diagonal elements
        if valley == 1
            H_n = diag([-2*Delta2, -2*Delta2, delta + gamma5 / 2 + Delta2, delta - 2 * Delta2]); 
        else
            H_n = diag([gamma2 / 2 + Delta2, gamma2 / 2 + Delta2, delta - 2 * Delta2, delta + gamma5 / 2 + Delta2]);
        end
        
        % off-diagonal elements
        H_n(2, 4) = - valley * 1j * x0 * sqrt(2);
        H_n(4, 2) = valley * 1j * x0 * sqrt(2);
        
        H_n(2, 3) = valley * 1j * x4 * sqrt(2);
        H_n(3, 2) = - valley * 1j * x4 * sqrt(2);
        
        H_n(3, 4) = gamma1;
        H_n(4, 3) = gamma1;
    else
        % diagonal elements
        if valley == 1
            H_n = diag([-2*Delta2, gamma2 / 2 + Delta2, delta + gamma5 / 2 + Delta2, delta - 2 * Delta2]);
        else
            H_n = diag([gamma2 / 2 + Delta2, -2*Delta2, delta - 2 * Delta2, delta + gamma5 / 2 + Delta2]);
        end
        
        % off-diagonal elements
        H_n(1, 4) = - valley * 1j * x0 * sqrt(2 * n);
        H_n(4, 1) = valley * 1j * x0 * sqrt(2 * n);
        
        H_n(1, 3) = valley * 1j * x4 * sqrt(2 * n);
        H_n(3, 1) = - valley * 1j * x4 * sqrt(2 * n);
        
        H_n(2, 3) = valley * 1j * x0 * sqrt(2 * (n - 1));
        H_n(3, 2) = - valley * 1j * x0 * sqrt(2 * (n - 1));
        
        H_n(2, 4) = - valley * 1j * x4 * sqrt(2 * (n - 1));
        H_n(4, 2) = valley * 1j * x4 * sqrt(2 * (n - 1));
        
        H_n(3, 4) = gamma1;
        H_n(4, 3) = gamma1;
    end
end

function W_n = get_off_diagonal_W_n(x3, valley, n)
    W_n = zeros(4);
    if n == 2
        W_n(2, 2) = valley * 1j * x3 * 2;
    else
        W_n(1, 2) = valley * 1j * x3 * sqrt(2 * n);
    end
end