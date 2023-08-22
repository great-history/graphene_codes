function Ham_LL_valley = construct_bilayer_LL_four_bands(x0, x3, x4, gamma1, delta, u, valley, LL_index_cutoff, dims)
    % LL_index_cutoff = 10 说明最高LL的指标为10，即只考虑LL_index = 0,1,2,3,4,5,6,7,8,,9,10，维数为dims = 2 * LL_index_cutoff
    % valley : + 1 for K valley ; -1 for Kp valley
    % in K valley, the Landau levels in the order : [|0,B2>, |1,B2>, |0,B1>, |0,A2> // |2,B2>, |0,A1>, |1,B1>, |1,A2> // |3,B2>, |1,A1>, |2,B1>, |2,A2> // ... ...]
    % in Kp valley, the Landau levels in the order : [|0,A1>, |1,A1>, |0,A2>, |0,B1> // |2,A1>, |0,B2>, |1,A2>, |1,B1> // |3,A1>, |1,B2>, |2,A2>, |2,B1> // ... ...]
    
    % 构建在谷valley处的LL能级哈密顿量
    Ham_LL_valley = zeros(dims);
    
    % 对 n = 1 block进行单独处理
    n = 1; 
    index_start = 4 * (n - 1) + 1;
    index_end = 4 * n;
    
    % H_n
    H_n = get_diagonal_H_n(x0, x4, gamma1, delta, u, valley, n);
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
        H_n = get_diagonal_H_n(x0, x4, gamma1, delta, u, valley, n);
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
        H_n = get_diagonal_H_n(x0, x4, gamma1, delta, u, valley, n);
        Ham_LL_valley(index_start:index_end, index_start:index_end) = H_n;
        Ham_LL_valley(index_start:index_end, index_start:index_end) = H_n;
    end
    
end

function H_n = get_diagonal_H_n(x0, x4, gamma1, delta, u, valley, n)
    if n == 1
        % diagonal elements
        H_n = valley * (u / 2) * diag([-1, -1, 1, -1]);
        H_n(3, 3) = H_n(3, 3) + delta;
        H_n(4, 4) = H_n(4, 4) + delta;
        
        % off-diagonal elements
        H_n(2, 4) = - valley * 1j * x0 * sqrt(2);
        H_n(4, 2) = valley * 1j * x0 * sqrt(2);
        
        H_n(2, 3) = - valley * 1j * x4 * sqrt(2);
        H_n(3, 2) = valley * 1j * x4 * sqrt(2);
        
        H_n(3, 4) = gamma1;
        H_n(4, 3) = gamma1;
    else
        % diagonal elements
        H_n = valley * u / 2 * diag([-1, 1, 1, -1]);
        H_n(3, 3) = H_n(3, 3) + delta;
        H_n(4, 4) = H_n(4, 4) + delta;
        
        % off-diagonal elements
        H_n(1, 4) = - valley * 1j * x0 * sqrt(2 * n);
        H_n(4, 1) = valley * 1j * x0 * sqrt(2 * n);
        
        H_n(1, 3) = - valley * 1j * x4 * sqrt(2 * n);
        H_n(3, 1) = valley * 1j * x4 * sqrt(2 * n);
        
        H_n(2, 3) = valley * 1j * x0 * sqrt(2 * (n - 1));
        H_n(3, 2) = - valley * 1j * x0 * sqrt(2 * (n - 1));
        
        H_n(2, 4) = valley * 1j * x4 * sqrt(2 * (n - 1));
        H_n(4, 2) = - valley * 1j * x4 * sqrt(2 * (n - 1));
        
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