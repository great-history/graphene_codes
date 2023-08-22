function [LL_K, LL_Kp] = trilayer_LLs_solver_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, N_LL, B_fields, B_steps, eps, no_mix, varargin)
    if no_mix == true
        var1 = varargin(1);
        Delta1s = var1{1};
        var2 = varargin(2);
        Delta1_steps = var2{1};
        
        [LL_K, LL_Kp] = trilayer_LLs_solver_with_D_no_mix(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_fields, B_steps, eps, Delta1s, Delta1_steps);
    else
        [LL_K, LL_Kp] = trilayer_LLs_solver_with_D_mix(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, N_LL, B_fields, B_steps, eps);
    end
end

function [LL_K, LL_Kp] = trilayer_LLs_solver_with_D_no_mix(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_fields, B_steps, eps, Delta1s, Delta1_steps)
    % 这个函数是先调用trilayer_LLs_asfo_Delta1,将monolayer的继承者与bilayer的继承者区分开来(故取名为mix)，这样画图的时候就知道了
    
    [LL_K_test, LL_Kp_test] = trilayer_LLs_asfo_Delta1(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_fields(1), Delta1s, Delta1_steps, eps);
    
    dims = 6 * N_LL + 9;
    LL_K = zeros(B_steps, 6*N_LL+9);
    LL_Kp = zeros(B_steps, 6*N_LL+9);
    
    % % 下面是开始计算固定D，随着B的演化
    % 先将第一个磁场先算出来，作为下一步的一个参考标准
    Delta1 = Delta1s(end);
    B_field = B_fields(1);
    [HK, HKp] = construct_HK_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, B_field, N_LL);
    
    [eigvec_HK_last, eigval_HK] = eig(HK);
    [eigvec_HKp_last, eigval_HKp] = eig(HKp);

    eigval_HK_diag_last = diag(eigval_HK);
    eigval_HKp_diag_last = diag(eigval_HKp);
    
    % 先做一下排序:以LL_K_test / LL_Kp_test为参考
    [eigvec_HK_last, eigval_HK_diag_last] = re_order_mix(LL_K_test(end,:), eigvec_HK_last, eigval_HK_diag_last, dims);
    [eigvec_HKp_last, eigval_HKp_diag_last] = re_order_mix(LL_Kp_test(end,:), eigvec_HKp_last, eigval_HKp_diag_last, dims);
    clear LL_K_test LL_Kp_test
    
    LL_K(1,:) = eigval_HK_diag_last;
    LL_Kp(1,:) = eigval_HKp_diag_last;
    
    for B_index = 2:B_steps
        B_field = B_fields(B_index);
        % construct Hamiltonian
        [HK, HKp] = construct_HK_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, B_field, N_LL);

        % call the eig sovler
        [eigvec_HK_now, eigval_HK] = eig(HK);
        eigval_HK_diag_now = diag(eigval_HK);
        
        [eigvec_HKp_now, eigval_HKp] = eig(HKp);
        eigval_HKp_diag_now = diag(eigval_HKp);
        
        [eigvec_HK_now, eigval_HK_diag_now, eig_num_K] = helper_re_order_states(eigvec_HK_last, eigval_HK_diag_last, ...
                                                                                eigvec_HK_now, eigval_HK_diag_now, dims, eps);
        [eigvec_HKp_now, eigval_HKp_diag_now, eig_num_Kp] = helper_re_order_states(eigvec_HKp_last, eigval_HKp_diag_last, ...
                                                                                   eigvec_HKp_now, eigval_HKp_diag_now, dims, eps);
        
        % push into the LLs
        LL_K(B_index,:) = eigval_HK_diag_now;
        LL_Kp(B_index,:) = eigval_HKp_diag_now;
        
        % 如果全是非简并，那么我们需要将这一次的结果，即eigvec_H_m_K_now/eigval_H_m_K_diag_now和eigvec_H_m_Kp_now/eigval_H_m_Kp_diag_now用到下一次
        % 如果存在简并的情形，那么我们需要将上一次的结果，即eigvec_H_m_K_last/eigval_H_m_K_diag_last和eigvec_H_m_Kp_last/eigval_H_m_Kp_diag_last用到下一次
        if (eig_num_K == dims)
           eigvec_HK_last = eigvec_HK_now;
           eigval_HK_diag_last = eigval_HK_diag_now;
        end

        if (eig_num_Kp == dims)
           eigvec_HKp_last = eigvec_HKp_now;
           eigval_HKp_diag_last = eigval_HKp_diag_now;
        end
    end    
end

function [eigvec_last, eigval_last] = re_order_mix(eigval_test, eigvec_last, eigval_last, dims)
    % 先做一下排序:以LL_K_test / LL_Kp_test为参考
    new_order = zeros(dims,1);
    for i = 1:dims
        val = eigval_test(end, i);
        for j = 1:dims
            if val == eigval_last(j)
                new_order(i) = j;
            end
        end
    end
    
    % 交换位置
    eigvec_last(:,:) = eigvec_last(:,new_order);
    eigval_last(:) = eigval_last(new_order);
end

function [LL_K, LL_Kp] = trilayer_LLs_solver_with_D_mix(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, N_LL, B_fields, B_steps, eps)
    dims = 6 * N_LL + 9;
    LL_K = zeros(B_steps, 6*N_LL+9);
    LL_Kp = zeros(B_steps, 6*N_LL+9);
    
    % 先将第一个磁场先算出来，作为下一步的一个参考标准
    B_field = B_fields(1);
    [HK, HKp] = construct_HK_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, B_field, N_LL);
    
    [eigvec_HK_last, eigval_HK] = eig(HK);
    [eigvec_HKp_last, eigval_HKp] = eig(HKp);

    eigval_HK_diag_last = diag(eigval_HK);
    eigval_HKp_diag_last = diag(eigval_HKp);

    LL_K(1,:) = eigval_HK_diag_last;
    LL_Kp(1,:) = eigval_HKp_diag_last;
    
    for B_index = 2:B_steps
        B_field = B_fields(B_index);
        % construct Hamiltonian
        [HK, HKp] = construct_HK_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, B_field, N_LL);

        % call the eig sovler
        [eigvec_HK_now, eigval_HK] = eig(HK);
        eigval_HK_diag_now = diag(eigval_HK);
        
        [eigvec_HKp_now, eigval_HKp] = eig(HKp);
        eigval_HKp_diag_now = diag(eigval_HKp);
%         [eigvec_HK_now, eigval_HK_diag_now, eigvec_HKp_now, eigval_HKp_diag_now] = re_order_trilayer(eigvec_HK, eigval_HK_diag, eigvec_HK_now, eigval_HK_diag_now, ...
%                                                                                                          eigvec_HKp, eigval_HKp_diag, eigvec_HKp_now, eigval_HKp_diag_now, ...
%                                                                                                          N_LL, err);
        
        [eigvec_HK_now, eigval_HK_diag_now, eig_num_K] = helper_re_order_states(eigvec_HK_last, eigval_HK_diag_last, ...
                                                                                eigvec_HK_now, eigval_HK_diag_now, dims, eps);
        [eigvec_HKp_now, eigval_HKp_diag_now, eig_num_Kp] = helper_re_order_states(eigvec_HKp_last, eigval_HKp_diag_last, ...
                                                                                   eigvec_HKp_now, eigval_HKp_diag_now, dims, eps);
        
        

        % push into the LLs
        LL_K(B_index,:) = eigval_HK_diag_now;
        LL_Kp(B_index,:) = eigval_HKp_diag_now;
        
        % 如果全是非简并，那么我们需要将这一次的结果，即eigvec_H_m_K_now/eigval_H_m_K_diag_now和eigvec_H_m_Kp_now/eigval_H_m_Kp_diag_now用到下一次
        % 如果存在简并的情形，那么我们需要将上一次的结果，即eigvec_H_m_K_last/eigval_H_m_K_diag_last和eigvec_H_m_Kp_last/eigval_H_m_Kp_diag_last用到下一次
        if (eig_num_K == dims)
           eigvec_HK_last = eigvec_HK_now;
           eigval_HK_diag_last = eigval_HK_diag_now;
        end

        if (eig_num_Kp == dims)
           eigvec_HKp_last = eigvec_HKp_now;
           eigval_HKp_diag_last = eigval_HKp_diag_now;
        end
        
%         eigvec_HK = eigvec_HK_now(:,:);
%         eigval_HK_diag = eigval_HK_diag_now(:);
%         eigvec_HKp = eigvec_HKp_now(:,:);
%         eigval_HKp_diag = eigval_HKp_diag_now(:);
    end    
end

%% 下面是旧版本的re_order函数，可能会稍快写
% function [eigvec_HK_now, eigval_HK_diag_now, eigvec_HKp_now, eigval_HKp_diag_now] = re_order_trilayer_old(eigvec_HK_last, eigvec_HK_now, eigval_HK_diag_now, ...
%                                                                                                      eigvec_HKp_last, eigvec_HKp_now, eigval_HKp_diag_now, ...
%                                                                                                      N_LL, eps)
%     thres = 1 - eps;
%     for i = 1:(6*N_LL+9)
%         vec = eigvec_HK_last(:,i);
%         val = eigval_HK_diag_now(i);
%         for j = i:(6*N_LL+9)
%             if abs(val - eigval_HK_diag_now(j)) < 0.1 
%                 inner_product = abs(dot(vec, eigvec_HK_now(:,j)));
%                 if inner_product > thres
% %                     if ~(j == i)
% %                         j
% %                         inner_product
% %                     end
%                     vec = eigvec_HK_now(:,j);
%                     eigvec_HK_now(:,j) = eigvec_HK_now(:,i);
%                     eigvec_HK_now(:,i) = vec;
% 
%                     eigval_HK_diag_now(i) = eigval_HK_diag_now(j);
%                     eigval_HK_diag_now(j) = val;
% 
%                     break
%                 end
%             end
%         end
%         
%         vec = eigvec_HKp_last(:,i);
%         val = eigval_HKp_diag_now(i);
%         
%         for k = i:(6*N_LL+9)
%             if abs(val - eigval_HKp_diag_now(k)) < 0.1
%                 inner_product = abs(dot(vec, eigvec_HKp_now(:,k)));
%                 if inner_product > thres
% %                     if ~(k == i)
% %                         k
% %                         inner_product
% %                     end
%                     vec = eigvec_HKp_now(:,k);
%                     eigvec_HKp_now(:,k) = eigvec_HKp_now(:,i);
%                     eigvec_HKp_now(:,i) = vec;
% 
%                     eigval_HKp_diag_now(i) = eigval_HKp_diag_now(k);
%                     eigval_HKp_diag_now(k) = val;
% 
%                     break
%                 end
%             end
%             
%         end
%         
%     end
% end
% 
% 
% function [eigvec_HK_now, eigval_HK_diag_now, eigvec_HKp_now, eigval_HKp_diag_now] = re_order_trilayer(eigvec_HK_last, eigval_HK_diag_last, eigvec_HK_now, eigval_HK_diag_now, ...
%                                                                                                           eigvec_HKp_last, eigval_HKp_diag_last, eigvec_HKp_now, eigval_HKp_diag_now, ...
%                                                                                                           N_LL, err)
%     dims = 6 * N_LL + 9;
%     new_order_K = zeros(dims,1);
%     new_order_Kp = zeros(dims,1);
%     
%     for idx = 1:dims
%         vec_K = eigvec_HK_last(:,idx);
%         val_K = eigval_HK_diag_last(idx);
%         
%         vec_Kp = eigvec_HKp_last(:,idx);
%         val_Kp = eigval_HKp_diag_last(idx);
%         
%         % 地毯式搜索
%         max_id_K = 0;
%         target_K = 1;
%         max_id_Kp = 0;
%         target_Kp = 1;
%         for jj = 1:dims
%             if abs(eigval_HK_diag_now(jj) - val_K) <= err
%                 inner_dot = abs(dot(vec_K, eigvec_HK_now(:,jj)));
%                 if inner_dot > max_id_K
%                     max_id_K = inner_dot;
%                     target_K = jj;
%                 end
%             end
%             
%             if abs(eigval_HKp_diag_now(jj) - val_Kp) <= err
%                 inner_dot = abs(dot(vec_Kp, eigvec_HKp_now(:,jj)));
%                 if inner_dot > max_id_Kp
%                     max_id_Kp = inner_dot;
%                     target_Kp = jj;
%                 end
%             end
%         end
%         new_order_K(idx) = target_K;
%         new_order_Kp(idx) = target_Kp;
%     end
%     
%      % 交换位置
%     eigvec_HK_now(:,:) = eigvec_HK_now(:,new_order_K);
%     eigval_HK_diag_now(:) = eigval_HK_diag_now(new_order_K);
%     
%     eigvec_HKp_now(:,:) = eigvec_HKp_now(:,new_order_Kp);
%     eigval_HKp_diag_now(:) = eigval_HKp_diag_now(new_order_Kp);
% end

%% 下面是测试代码
% example:
% gamma0 = 3100;
% gamma1 = 355;
% gamma2 = -20.7904;
% gamma3 = 346.7153;  
% gamma4 = 91.2578;
% gamma5 = 35.7223;
% D_field = 1;
% delta = 1.5+(gamma5-gamma2)/2;
% Delta1 = D_field * 0.1 * 1000;
% Delta2 = 1.8;
% N_LL = 30;
% B_start = 0.05;
% B_end = 6;
% B_steps = 400;
% B_fields = linspace(B_start, B_end, B_steps);
% eps = 0.1

% [LL_K, LL_Kp] = trilayer_LLs_solver_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, N_LL, B_fields, B_steps, eps);