function [Ham_m_LL_valley, Ham_b_LL_valley, D_matrix] = construct_trilayer_ABA_LL_six_bands(x0, x3, x4, gamma1, delta, gamma2, gamma5, Delta1, Delta2, ...
                                                                                            valley, LL_index_cutoff, dims_m, dims_b)
    %% 该函数只适用于不加电位移场时的计算，加入电位移场后最好是用six_band model
    % 由于该函数不考虑电位移场的影响，因此也就不需要让单层与双层的维数一致了
    
    %% 函数说明
    % LL_index_cutoff = 10 说明最高LL的指标为10，即只考虑LL_index = 0,1,2,3,4,5,6,7,8,,9,10，即每个分量上的LL_index最大不能超过LL_index_cutoff
    % 对于trilayer ABA而言，通过基矢的变换，可以得到:
    % monolayer-like brach : |phi1> = 1 / sqrt(2) * 【|A1> - |A3>】 ; |phi2> = 1 / sqrt(2) * 【|B1> - |B3>】
    % bilayer-like brach : |phi3> = 1 / sqrt(2) * 【|A1> + |A3>】 ; |phi4> = 1 / sqrt(2) * 【|B1> + |B3>】 ; |phi5> = |A2> ; |phi6> = |B2>
    
    % 对于单层，维数为dims = 2 * LL_index_cutoff + 1 ; 对于双层，维数为dims = 4 * LL_index_cutoff

    % valley : + 1 for K valley ; -1 for Kp valley
    % 对于monolayer graphene
    % in K valley, the Landau levels in the order : [|0,phi2> // |1,phi2>, |0,phi1> // |2,phi2>, |1,phi1> // ... ...]
    % in Kp valley, the Landau levels in the order : [|0,phi1> // |1,phi1>, |0,phi2> // |2,phi1>, |1,phi2> // ... ...]
    % 对于bilayer graphene
    % in K valley, the Landau levels in the order : [|0,phi6>, |1,phi6>, |0,phi4>, |0,phi5> // |2,phi6>, |0,phi3>, |1,phi4>, |1,phi5> // |3,phi6>, |1,phi3>, |2,phi4>, |2,phi5> // ... ...]
    % in Kp valley, the Landau levels in the order : [|0,phi3>, |1,phi3>, |0,phi5>, |0,phi4> // |2,phi3>, |0,phi6>, |1,phi5>, |1,phi4> // |3,phi3>, |1,phi6>, |2,phi5>, |2,phi4> // ... ...]
    
    %% 注意点：
    % 传入gamma1之前应将gamma1扩大sqrt(2)倍
    % 传入x3之前应将gamma3扩大sqrt(2)倍
    % 传入x4之前应将gamma4扩大sqrt(2)倍
    
    % dims_m = 2 * LL_index_cutoff + 1;
    % dims_b = 4 * LL_index_cutoff;
    
    Ham_m_LL_valley = construct_monolayer_LL_for_ABA(x0, delta, gamma2, gamma5, Delta2, valley, LL_index_cutoff, dims_m);
    Ham_b_LL_valley = construct_bilayer_LL_four_bands_for_ABA(x0, x3, x4, gamma1, delta, gamma2, gamma5, Delta2, valley, LL_index_cutoff, dims_b);
    
    % Ham_LL_valley = zeros(dims);
    % Ham_LL_valley(1:dims_m, 1:dims_m) = Ham_m_LL_valley;
    % Ham_LL_valley((dims_m + 1):end, (dims_m + 1):end) = Ham_b_LL_valley;
    
    %% 加入电位移场的影响 Delta1
    if ~(Delta1 == 0)
        D_matrix = add_D_field(Delta1, valley, LL_index_cutoff, dims_m, dims_b);
        % Ham_LL_valley((dims_m + 1):end, 1:dims_m) = D_matrix;
        % Ham_LL_valley(1:dims_m, (dims_m + 1):end) = D_matrix'; 
    else
        D_matrix = zeros(dims_b, dims_m);
    end
    
end

function D_matrix = add_D_field(Delta1, valley, LL_index_cutoff, dims_m, dims_b)
    % 加入电位移场之后， phi1与phi3之间有一一对应，phi2与phi4之间有一一对应
    D_matrix = zeros(dims_b, dims_m);
    if valley == 1
        % phi2与phi4之间
        D_matrix(3,1) = Delta1;
        for n = 1:(LL_index_cutoff - 1)
            D_matrix(4 * n + 3, 2 * n) = Delta1;
        end
        
        % phi1与phi3之间
        for n = 1:(LL_index_cutoff - 1)
            D_matrix(4 * n + 2, 2 * n + 1) = Delta1;
        end
    else
        % phi1与phi3之间
        D_matrix(1,1) = Delta1;
        D_matrix(2,2) = Delta1;
        for n = 1:(LL_index_cutoff - 1)
            D_matrix(4 * n + 1, 2 * (n + 1)) = Delta1;
        end
        
        % phi2与phi4之间
        for n = 1:LL_index_cutoff
            D_matrix(4 * n, 2 * n + 1) = Delta1;
        end
    end
end