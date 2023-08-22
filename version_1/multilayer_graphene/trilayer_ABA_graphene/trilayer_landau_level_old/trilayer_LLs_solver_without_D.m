function [LL_K_m, LL_Kp_m, LL_K_b, LL_Kp_b] = trilayer_LLs_solver_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_fields, B_steps, eps)
    dims_m = 2 * N_LL + 3;
    dims_b = 4 * N_LL + 6;
    LL_K_m = zeros(B_steps, dims_m);
    LL_Kp_m = zeros(B_steps, dims_m);
    LL_K_b = zeros(B_steps, dims_b);
    LL_Kp_b = zeros(B_steps, dims_b);
    
    % 先将第一个磁场先算出来，作为下一步的一个参考标准
    B_field = B_fields(1);
    [HK_b, HK_m, HKp_b, HKp_m] = construct_HK_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, B_field, N_LL);
    [eigvec_HK_m_last, eigval_HK_m] = eig(HK_m);
    [eigvec_HKp_m_last, eigval_HKp_m] = eig(HKp_m);
    [eigvec_HK_b_last, eigval_HK_b] = eig(HK_b);
    [eigvec_HKp_b_last, eigval_HKp_b] = eig(HKp_b);

    eigval_HK_m_diag_last = diag(eigval_HK_m);
    eigval_HKp_m_diag_last = diag(eigval_HKp_m);
    eigval_HK_b_diag_last = diag(eigval_HK_b);
    eigval_HKp_b_diag_last = diag(eigval_HKp_b);

    LL_K_m(1,:) = eigval_HK_m_diag_last;
    LL_Kp_m(1,:) = eigval_HKp_m_diag_last;
    LL_K_b(1,:) = eigval_HK_b_diag_last;
    LL_Kp_b(1,:) = eigval_HKp_b_diag_last;
    
    for B_index = 2:B_steps
        B_field = B_fields(B_index);
        % construct Hamiltonian
        [HK_b, HK_m, HKp_b, HKp_m] = construct_HK_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, B_field, N_LL);

        % call the eig sovler
        [eigvec_HK_m_now, eigval_HK_m] = eig(HK_m);
        [eigvec_HKp_m_now, eigval_HKp_m] = eig(HKp_m);

        [eigvec_HK_b_now, eigval_HK_b] = eig(HK_b);
        [eigvec_HKp_b_now, eigval_HKp_b] = eig(HKp_b);

        eigval_HK_m_diag_now = diag(eigval_HK_m);
        eigval_HKp_m_diag_now = diag(eigval_HKp_m);
        eigval_HK_b_diag_now = diag(eigval_HK_b);
        eigval_HKp_b_diag_now = diag(eigval_HKp_b);


        % re-order of bilayer
%         [eigvec_HK_b_now, eigval_HK_b_diag_now, eigvec_HKp_b_now, eigval_HKp_b_diag_now] = re_order_bilayer(eigvec_HK_b, eigvec_HK_b_now, eigval_HK_b_diag_now, ...
%                                                                                                             eigvec_HKp_b, eigvec_HKp_b_now, eigval_HKp_b_diag_now, N_LL, eps);

%         [eigvec_HK_m_now, eigval_HK_m_diag_now, eigvec_HKp_m_now, eigval_HKp_m_diag_now] = re_order_monolayer(eigvec_HK_m, eigvec_HK_m_now, eigval_HK_m_diag_now, ...
%                                                                                                               eigvec_HKp_m, eigvec_HKp_m_now, eigval_HKp_m_diag_now, N_LL, eps);
        
        [eigvec_HK_m_now, eigval_HK_m_diag_now, eig_num_K_m] = helper_re_order_states(eigvec_HK_m_last, eigval_HK_m_diag_last, ...
                                                                                      eigvec_HK_m_now, eigval_HK_m_diag_now, dims_m, eps);
        [eigvec_HKp_m_now, eigval_HKp_m_diag_now, eig_num_Kp_m] = helper_re_order_states(eigvec_HKp_m_last, eigval_HKp_m_diag_last, ...
                                                                                         eigvec_HKp_m_now, eigval_HKp_m_diag_now, dims_m, eps);
                                                                                     
        [eigvec_HK_b_now, eigval_HK_b_diag_now, eig_num_K_b] = helper_re_order_states(eigvec_HK_b_last, eigval_HK_b_diag_last, ...
                                                                                      eigvec_HK_b_now, eigval_HK_b_diag_now, dims_b, eps);
        [eigvec_HKp_b_now, eigval_HKp_b_diag_now, eig_num_Kp_b] = helper_re_order_states(eigvec_HKp_b_last, eigval_HKp_b_diag_last, ...
                                                                                         eigvec_HKp_b_now, eigval_HKp_b_diag_now, dims_b, eps);
        
        % push into the LLs
        LL_K_m(B_index,:) = eigval_HK_m_diag_now;
        LL_Kp_m(B_index,:) = eigval_HKp_m_diag_now;
        LL_K_b(B_index,:) = eigval_HK_b_diag_now;
        LL_Kp_b(B_index,:) = eigval_HKp_b_diag_now;
        
        % 如果全是非简并，那么我们需要将这一次的结果，即eigvec_H_m_K_now/eigval_H_m_K_diag_now和eigvec_H_m_Kp_now/eigval_H_m_Kp_diag_now用到下一次
        % 如果存在简并的情形，那么我们需要将上一次的结果，即eigvec_H_m_K_last/eigval_H_m_K_diag_last和eigvec_H_m_Kp_last/eigval_H_m_Kp_diag_last用到下一次
        if (eig_num_K_m == dims_m)
           eigvec_HK_m_last = eigvec_HK_m_now;
           eigval_HK_m_diag_last = eigval_HK_m_diag_now;
        end

        if (eig_num_Kp_m == dims_m)
           eigvec_HKp_m_last = eigvec_HKp_m_now;
           eigval_HKp_m_diag_last = eigval_HKp_m_diag_now;
        end
        
        if (eig_num_K_b == dims_b)
           eigvec_HK_b_last = eigvec_HK_b_now;
           eigval_HK_b_diag_last = eigval_HK_b_diag_now;
        end

        if (eig_num_Kp_b == dims_b)
           eigvec_HKp_b_last = eigvec_HKp_b_now;
           eigval_HKp_b_diag_last = eigval_HKp_b_diag_now;
        end
        % eigvec_HK_m = eigvec_HK_m_now(:,:);
        % eigvec_HKp_m = eigvec_HKp_m_now(:,:);
        % eigvec_HK_b = eigvec_HK_b_now(:,:);
        % eigvec_HKp_b = eigvec_HKp_b_now(:,:);
    end
    
end

%% 下面这两个re_order函数可能算起来会稍快，但我还是觉得用helper_re_order_states更好
% function [eigvec_HK_b_now, eigval_HK_b_diag_now, eigvec_HKp_b_now, eigval_HKp_b_diag_now] = re_order_bilayer(eigvec_HK_b_last, eigvec_HK_b_now, eigval_HK_b_diag_now, ...
%                                                                                                              eigvec_HKp_b_last, eigvec_HKp_b_now, eigval_HKp_b_diag_now, ...
%                                                                                                              N_LL, eps)
%     thres = 1 - eps;
%     for i = 1:(4*N_LL+6)
%         vec = eigvec_HK_b_last(:,i);
%         for j = i:(4*N_LL+6)
%             inner_product = abs(dot(vec, eigvec_HK_b_now(:,j)));
%             if inner_product > thres
%                 vec = eigvec_HK_b_now(:,j);
%                 eigvec_HK_b_now(:,j) = eigvec_HK_b_now(:,i);
%                 eigvec_HK_b_now(:,i) = vec;
%                 
%                 val = eigval_HK_b_diag_now(j);
%                 eigval_HK_b_diag_now(j) = eigval_HK_b_diag_now(i);
%                 eigval_HK_b_diag_now(i) = val;
%                 
%                 break
%             end
%         end
%         
%         vec = eigvec_HKp_b_last(:,i);
%         for k = i:(4*N_LL+6)
%             inner_product = abs(dot(vec, eigvec_HKp_b_now(:,k)));
%             if inner_product > thres
%                 vec = eigvec_HKp_b_now(:,k);
%                 eigvec_HKp_b_now(:,k) = eigvec_HKp_b_now(:,i);
%                 eigvec_HKp_b_now(:,k) = vec;
%                 
%                 val = eigval_HKp_b_diag_now(k);
%                 eigval_HKp_b_diag_now(k) = eigval_HKp_b_diag_now(i);
%                 eigval_HKp_b_diag_now(i) = val;
%                 
%                 break
%             end
%         end
%         
%     end
% end
% 
% function [eigvec_HK_m_now, eigval_HK_m_diag_now, eig_num_K, eigvec_HKp_m_now, eigval_HKp_m_diag_now, eig_num_Kp] = ...
%             re_order_monolayer(eigvec_HK_m_last, eigval_HK_m_diag_last, eigvec_HK_m_now, eigval_HK_m_diag_now, ...
%                                eigvec_HKp_m_last, eigval_HKp_m_diag_last, eigvec_HKp_m_now, eigval_HKp_m_diag_now, ...
%                                N_LL, eps)
%     % 要保证eigval_H_m_K_last和eigval_H_m_Kp_last都是非简并的才行！！！
%     dims = 2*N_LL + 3;
%     % 首先检查eigval_H_m_K_last和eigval_H_m_Kp_last是否都是非简并的
%     [eigvals_K_unique, degeneracy_K, eig_num_K, first_idxs_K, last_idxs_K] = get_degeneracy(eigval_HK_m_diag_now, dims);
%     [eigvals_Kp_unique, degeneracy_Kp, eig_num_Kp, first_idxs_Kp, last_idxs_Kp] = get_degeneracy(eigval_HKp_m_diag_now, dims);
%     
%     % 分流：eigvec_HK_m_now中简并度为1的本征态与eigvec_HK_m_last中简并度为1的本征态有一一对应关系，即内积有且仅有一个为1，其余全为0
%     % 而eigvec_HK_m_now中多重简并态与eigvec_HK_m_last中的多个本征态有交叠，即内积可能有多个不为0，但都小于1
%     
%     new_order_K = get_new_order(eigvec_HK_m_last, eigval_HK_m_diag_last, eigvec_HK_m_now, eigvals_K_unique, degeneracy_K, first_idxs_K, last_idxs_K, dims, eig_num_K, eps);
%     new_order_Kp = get_new_order(eigvec_HKp_m_last, eigval_HKp_m_diag_last, eigvec_HKp_m_now, eigvals_Kp_unique, degeneracy_Kp, first_idxs_Kp, last_idxs_Kp, dims, eig_num_Kp, eps);
%     
%      % 交换位置
%     eigvec_HK_m_now(:,:) = eigvec_HK_m_now(:,new_order_K);
%     eigval_HK_m_diag_now(:) = eigval_HK_m_diag_now(new_order_K);
%     
%     eigvec_HKp_m_now(:,:) = eigvec_HKp_m_now(:,new_order_Kp);
%     eigval_HKp_m_diag_now(:) = eigval_HKp_m_diag_now(new_order_Kp);
% end


% function [eigvec_HK_m_now, eigval_HK_m_diag_now, eigvec_HKp_m_now, eigval_HKp_m_diag_now] = re_order_monolayer_old(eigvec_HK_m_last, eigvec_HK_m_now, eigval_HK_m_diag_now, ...
%                                                                                                                eigvec_HKp_m_last, eigvec_HKp_m_now, eigval_HKp_m_diag_now, ...
%                                                                                                                N_LL, eps)
%     thres = 1 - eps;                                                                            
%     for i = 1:(2*N_LL+3)
%         for j = i:(2*N_LL+3)
%             inner_product = abs(dot(eigvec_HK_m_last(:,i), eigvec_HK_m_now(:,j)));
%             if inner_product > thres
%                 vec = eigvec_HK_m_now(:,j);
%                 eigvec_HK_m_now(:,j) = eigvec_HK_m_now(:,i);
%                 eigvec_HK_m_now(:,i) = vec;
%                 
%                 val = eigval_HK_m_diag_now(j);
%                 eigval_HK_m_diag_now(j) = eigval_HK_m_diag_now(i);
%                 eigval_HK_m_diag_now(i) = val;
%                 break
%             end
%         end
%         
%         for k = i:(2*N_LL+3)
%             inner_product = abs(dot(eigvec_HKp_m_last(:,i), eigvec_HKp_m_now(:,k)));
%             if inner_product > thres
%                 vec = eigvec_HKp_m_now(:,k);
%                 eigvec_HKp_m_now(:,k) = eigvec_HKp_m_now(:,i);
%                 eigvec_HKp_m_now(:,i) = vec;
%                 
%                 val = eigval_HKp_m_diag_now(k);
%                 eigval_HKp_m_diag_now(k) = eigval_HKp_m_diag_now(i);
%                 eigval_HKp_m_diag_now(i) = val;
%                 
%                 break
%             end
%         end
%         
%     end
% end

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
% % example
% gamma0 = 3100;
% gamma1 = 355;
% gamma2 = -20.7904;
% gamma3 = 346.7153;  
% gamma4 = 91.2578;
% gamma5 = 35.7223;
% 
% delta = 1.5+(gamma5-gamma2)/2;
% Delta2 = 1.8;
% 
% N_LL = 30;
% 
% B_start = 0.05;
% B_end = 6;
% B_steps = 400;
% B_fields = linspace(B_start, B_end, B_steps);
% 
% eps = 0.1
% [LL_K_m, LL_Kp_m, LL_K_b, LL_Kp_b] =
% trilayer_LLs_solver_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_fields, B_steps, eps);