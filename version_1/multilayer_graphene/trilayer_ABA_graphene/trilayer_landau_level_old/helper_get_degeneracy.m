function [eigvals_unique, degeneracy, eig_num, first_idxs, last_idxs] = helper_get_degeneracy(eigvals, dims)
    % 假定eigval_H_m_K_diag_last是从小到大排列的，所以返回的简并度可以很方便地转换到指标，
    % 比如[1,1,2,3,1] → [1,2,4,7,8]就是本征值的第一位置
    degeneracy = zeros(dims,1);
    eigvals_unique = zeros(dims,1);
    
    eig_val_now = eigvals(1);
    eig_num = 1;
    eigvals_unique(eig_num) = eig_val_now;
    count = 1;
    
    eps = 1e-8;
    for i = 2:dims
        eig_val_last = eig_val_now;
        eig_val_now = eigvals(i);
        if abs(eig_val_now - eig_val_last) <= eps
            count = count + 1;
        else
            degeneracy(eig_num) = count;
            eig_num = eig_num + 1;
            eigvals_unique(eig_num) = eig_val_now;
            count = 1;
        end
    end
    degeneracy(eig_num) = count;
    
    % pull out elements
    eigvals_unique = eigvals_unique(1:eig_num);
    degeneracy = degeneracy(1:eig_num);
    
    first_idxs = zeros(eig_num,1);
    last_idxs = zeros(eig_num,1);
    
    last = 0;
    for num = 1:eig_num
        first = last + 1;
        last = last + degeneracy(num); 
        first_idxs(num) = first;
        last_idxs(num) = last;
    end
end