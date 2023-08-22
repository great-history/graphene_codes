function [LL_K_m, LL_Kp_m] = monolayer_LLs_solver(gamma0, U1, U2, N_LL, B_fields, B_steps, eps)
    dims_m = 2 * N_LL + 3;
    LL_K_m = zeros(B_steps, dims_m);
    LL_Kp_m = zeros(B_steps, dims_m);

    % 先将第一个磁场先算出来，作为下一步的一个参考标准
    B_field = B_fields(1);
    E0 = (sqrt(2 * B_field) * gamma0 * 0.142 * 3) / (2 * 25.66);
    [HK_m, HKp_m] = construct_HK_monolayer(E0, U1, U2, N_LL);
    [eigvec_HK_m_last, eigval_HK_m] = eig(HK_m);
    [eigvec_HKp_m_last, eigval_HKp_m] = eig(HKp_m);
    
    eigval_HK_m_diag_last = diag(eigval_HK_m);
    eigval_HKp_m_diag_last = diag(eigval_HKp_m);

    LL_K_m(1,:) = eigval_HK_m_diag_last;
    LL_Kp_m(1,:) = eigval_HKp_m_diag_last;
    
    for B_index = 2:B_steps
        B_field = B_fields(B_index);
        % construct Hamiltonian
        E0 = (sqrt(2 * B_field) * gamma0 * 0.142 * 3) / (2 * 25.66);
        [HK_m, HKp_m] = construct_HK_monolayer(E0, U1, U2, N_LL);
        
        % call the eig sovler
        [eigvec_HK_m_now, eigval_HK_m] = eig(HK_m);
        [eigvec_HKp_m_now, eigval_HKp_m] = eig(HKp_m);

        eigval_HK_m_diag_now = diag(eigval_HK_m);
        eigval_HKp_m_diag_now = diag(eigval_HKp_m);

        
        [eigvec_HK_m_now, eigval_HK_m_diag_now, eig_num_K_m] = helper_re_order_states(eigvec_HK_m_last, eigval_HK_m_diag_last, ...
                                                                                      eigvec_HK_m_now, eigval_HK_m_diag_now, dims_m, eps);
        [eigvec_HKp_m_now, eigval_HKp_m_diag_now, eig_num_Kp_m] = helper_re_order_states(eigvec_HKp_m_last, eigval_HKp_m_diag_last, ...
                                                                                         eigvec_HKp_m_now, eigval_HKp_m_diag_now, dims_m, eps);
        
        % push into the LLs
        LL_K_m(B_index,:) = eigval_HK_m_diag_now;
        LL_Kp_m(B_index,:) = eigval_HKp_m_diag_now;
        
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
        
    end
end


% % example
% gamma0 = 3100;
% gamma1 = 390;
% gamma2 = -28;
% gamma3 = 315;
% gamma4 = 41;
% gamma5 = 50;
% 
% delta = 46;
% Delta1 = 0;
% Delta2 = 5.7;
% 
% % monolayer
% U5 = Delta2 - gamma2 / 2;
% U6 = Delta2 + delta - gamma5 / 2;
% 
% N_LL = 30;
% 
% B_start = 0.05;
% B_end = 14;
% B_steps = 1000;
% B_fields = linspace(B_start, B_end, B_steps);
% 
% eps = 0.005 * gamma0;  % 取一个能够大致描述能隙大小的值
% tic
% [LL_K_m, LL_Kp_m] = monolayer_LLs_solver(gamma0, U5, U6, N_LL, B_fields, B_steps, eps);
% toc
%
% figure
% axis([B_start B_end -100 100])
% hold on
% 
% dims_m = 2 * N_LL + 3;
% 
% % plot monolayer-like LL
% for i = 1:dims_m
%     plot(B_fields, LL_K_m(:,i), 'r', 'LineWidth', 0.5)
%     plot(B_fields, LL_Kp_m(:,i), 'r--', 'LineWidth', 1.0)
% end
% 
% grid on