function [LL_K, LL_Kp] = trilayer_LLs_asfo_Delta1(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_field, Delta1s, Delta1_steps, eps)
    dims = 6 * N_LL + 9;
    LL_K = zeros(Delta1_steps, dims);
    LL_Kp = zeros(Delta1_steps, dims);
    
    % 先将第一个D场先算出来，作为下一步的一个参考标准
    Delta1 = Delta1s(1);
    if Delta1 == 0  % 此时单层与双层解耦
        [HK_b, HK_m, HKp_b, HKp_m] = construct_HK_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, B_field, N_LL);
        
        [eigvec_HK_m, eigval_HK_m] = eig(HK_m);
        [eigvec_HKp_m, eigval_HKp_m] = eig(HKp_m);
        [eigvec_HK_b, eigval_HK_b] = eig(HK_b);
        [eigvec_HKp_b, eigval_HKp_b] = eig(HKp_b);

        eigval_HK_m_diag = diag(eigval_HK_m);
        eigval_HKp_m_diag = diag(eigval_HKp_m);
        eigval_HK_b_diag = diag(eigval_HK_b);
        eigval_HKp_b_diag = diag(eigval_HKp_b);
        
        eigvec_HK_last = [eigvec_HK_m, zeros(2 * N_LL + 3, 4 * N_LL + 6); zeros(4 * N_LL + 6, 2 * N_LL + 3), eigvec_HK_b];
        eigvec_HKp_last = [eigvec_HKp_m, zeros(2 * N_LL + 3, 4 * N_LL + 6); zeros(4 * N_LL + 6, 2 * N_LL + 3), eigvec_HKp_b];
        eigval_HK_diag_last = [eigval_HK_m_diag; eigval_HK_b_diag];
        eigval_HKp_diag_last = [eigval_HKp_m_diag; eigval_HKp_b_diag];
    else
        [HK, HKp] = construct_HK_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, B_field, N_LL);
        
        [eigvec_HK_last, eigval_HK] = eig(HK);
        [eigvec_HKp_last, eigval_HKp] = eig(HKp);
        
        eigval_HK_diag_last = diag(eigval_HK);
        eigval_HKp_diag_last = diag(eigval_HKp);
    end

    LL_K(1,:) = eigval_HK_diag_last;
    LL_Kp(1,:) = eigval_HKp_diag_last;
    
    for D_index = 2:Delta1_steps
        Delta1 = Delta1s(D_index);
        % construct Hamiltonian
        [HK, HKp] = construct_HK_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, B_field, N_LL);

        % call the eig sovler
        [eigvec_HK_now, eigval_HK] = eig(HK);
        [eigvec_HKp_now, eigval_HKp] = eig(HKp);

        eigval_HK_diag_now = diag(eigval_HK);
        eigval_HKp_diag_now = diag(eigval_HKp);

        % % re-order of trlayer
        % [eigvec_HK_now, eigval_HK_diag_now] = re_order_trilayer(eigvec_HK, eigval_HK_diag, eigvec_HK_now, eigval_HK_diag_now, N_LL, err);
        % [eigvec_HKp_now, eigval_HKp_diag_now] = re_order_trilayer(eigvec_HKp, eigval_HKp_diag, eigvec_HKp_now, eigval_HKp_diag_now, N_LL, err);
        
        [eigvec_HK_now, eigval_HK_diag_now, eig_num_K] = helper_re_order_states(eigvec_HK_last, eigval_HK_diag_last, ...
                                                                                eigvec_HK_now, eigval_HK_diag_now, dims, eps);
        [eigvec_HKp_now, eigval_HKp_diag_now, eig_num_Kp] = helper_re_order_states(eigvec_HKp_last, eigval_HKp_diag_last, ...
                                                                                   eigvec_HKp_now, eigval_HKp_diag_now, dims, eps);
        
        % push into the LLs
        LL_K(D_index,:) = eigval_HK_diag_now;
        LL_Kp(D_index,:) = eigval_HKp_diag_now;
        
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

function [eigvec_HK_now, eigval_HK_diag_now] = re_order_trilayer(eigvec_HK_last, eigval_HK_diag_last, eigvec_HK_now, eigval_HK_diag_now, N_LL, err)
%     thres = 1 - eps;
    dims = 6 * N_LL + 9;
    new_order = zeros(dims,1);
    
    for idx = 1:(6*N_LL+9)
        vec = eigvec_HK_last(:,idx);
        val = eigval_HK_diag_last(idx);
        
        % 地毯式搜索
        max_id = 0;
        target = 1;
        for jj = 1:(6*N_LL+9)
            if abs(eigval_HK_diag_now(jj) - val) <= err
                inner_dot = abs(dot(vec, eigvec_HK_now(:,jj)));
                if inner_dot > max_id
                    max_id = inner_dot;
                    target = jj;
                end
            end
        end
        new_order(idx) = target;
    end
    
     % 交换位置
    eigvec_HK_now(:,:) = eigvec_HK_now(:,new_order);
    eigval_HK_diag_now(:) = eigval_HK_diag_now(new_order);
end


% function [eigvec_HK_now, eigval_HK_diag_now] = re_order_trilayer_new(eigvec_HK_last, eigval_HK_diag_last, eigvec_HK_now, eigval_HK_diag_now, N_LL, eps, err)
%     dims = 6 * N_LL + 9;
%     new_order = zeros(dims,1);
%     inner_dots_K = zeros(6 * N_LL +9, 1);
%     
%     for idx = 1:(6*N_LL+9)
%         vec = eigvec_HK_last(:,idx);
%         val = eigval_HK_diag_last(idx);
%         upb = val + err;
%         dwb = val - err;
%         % 先用二分法找出搜索范围，利用eigval_HK_diag_now往往是从小到大排列的性质
%         if eigval_HK_diag_now(1) > upb
%             target = 1;
%         elseif eigval_HK_diag_now(end) < dwb
%             target = dims;
%         elseif eigval_HK_diag_now(1) > dwb
%             left = 1;
%             inner_dot = abs(dot(vec, eigvec_HK_now(:,left)));
%             inner_dots_K(left) = inner_dot;
%             
%             max_id = inner_dot;
%             target = left;
%             
%             right = 2;
%             while eigval_HK_diag_now(right) <= upb
%                 inner_dot = abs(dot(vec, eigvec_HK_now(:,right)));
%                 inner_dots_K(right) = inner_dot;
%                 if inner_dot > max_id
%                     max_id = inner_dot;
%                     target = right;
%                 end
%                 
%                 right = right + 1;
%             end
%             
%         elseif eigval_HK_diag_now(end) < upb
%             right = dims;
%             inner_dot = abs(dot(vec, eigvec_HK_now(:,right)));
%             inner_dots_K(right) = inner_dot;
%             
%             
%         else
%             left = 1;
%             right = 6*N_LL + 9;
%         end
%         
%         % push target into new_order
%         new_order(idx) = target;
%     end
% 
%     % 交换位置
%     eigvec_HK_now(:,:) = eigvec_HK_now(:,new_order);
%     eigval_HK_diag_now(:) = eigval_HK_diag_now(new_order);
% end

% example:
% gamma0 = 3100;
% gamma1 = 355;
% gamma2 = -20.7904;
% gamma3 = 346.7153;  
% gamma4 = 91.2578;
% gamma5 = 35.7223;
% 
% Delta1_start = 0;
% Delta1_end = 250;
% Delta1_steps = 1000;
% Delta1s = linspace(Delta1_start, Delta1_end, Delta1_steps);
% 
% delta = 1.5+(gamma5-gamma2)/2;
% Delta2 = 1.8;
% N_LL = 30;
% B_field = 5;
% eps = 0.1;
% 
% [LL_K, LL_Kp] = trilayer_LLs_asfo_Delta1(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, N_LL, B_field, Delta1s, Delta1_steps, eps)