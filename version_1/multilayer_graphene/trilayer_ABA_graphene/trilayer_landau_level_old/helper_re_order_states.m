function [eigvecs_now, eigvals_now, eig_num] = ...
            helper_re_order_states(eigvecs_last, eigvals_last, eigvecs_now, eigvals_now, dims, eps)
    % 要保证eigvals_last和eigvals_last都是非简并的才行！！！
    
    % 首先检查eigval_H_m_K_last和eigval_H_m_Kp_last是否都是非简并的
    [eigvals_unique, degeneracy, eig_num, first_idxs, last_idxs] = helper_get_degeneracy(eigvals_now, dims);
    
    % 分流：eigvecs_now中简并度为1的本征态与eigvecs_last中简并度为1的本征态有一一对应关系，即内积有且仅有一个为1，其余全为0
    % 而eigvecs_now中多重简并态与eigvecs_last中的多个本征态有交叠，即内积可能有多个不为0，但都小于1
    
    new_order = helper_get_new_order(eigvecs_last, eigvals_last, eigvecs_now, eigvals_now, degeneracy, first_idxs, last_idxs, dims, eig_num, eps);
    
     % 交换位置
    eigvecs_now(:,:) = eigvecs_now(:,new_order);
    eigvals_now(:) = eigvals_now(new_order);
end