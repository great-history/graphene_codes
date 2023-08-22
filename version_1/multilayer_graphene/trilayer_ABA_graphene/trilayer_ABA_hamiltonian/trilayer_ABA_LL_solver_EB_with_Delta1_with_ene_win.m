function [eigvals_LL_K, eigvals_LL_Kp, eig_info_HK_select_cell, eig_info_HKp_select_cell] = ...
            trilayer_ABA_LL_solver_EB_with_Delta1_with_ene_win(v0, v3, v4, gamma1, delta, gamma2, gamma5, Delta1, Delta2, B_fields_list, B_steps, LL_index_cutoff, dims_m, dims_b, ene_ub, ene_lb)
    h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
    dims = dims_b + dims_m;
    
    Ham_LL_K = zeros(dims);
    Ham_LL_Kp = zeros(dims);
    eigvals_LL_K = zeros(B_steps, dims);
    eigvals_LL_Kp = zeros(B_steps, dims);
    
    % 用来存放energy window内的所有的LL本征态和本征能量
    eig_info_HK_select_cell = cell(B_steps, 4); % 第一个维度放eigvec_select, 第二个维度放eigval_select，第三个维度放本征态的个数, 第四个维度存放各个本征态对应的LL_index
    eig_info_HKp_select_cell = cell(B_steps, 4); % 第一个维度放eigvec_select, 第二个维度放eigval_select，第三个维度放本征态的个数, 第四个维度存放各个本征态对应的LL_index
    
    for B_index = 1:B_steps
        %% parameters set up
        B_field = B_fields_list(B_index);
        mag_length = 25.66 / sqrt(B_field); % 以nm为单位
        
        x0 = h_bar * v0 / mag_length * 10^9; % 以eV为单位
        x3 = h_bar * v3 / mag_length * 10^9; % 以eV为单位
        x4 = h_bar * v4 / mag_length * 10^9; % 以eV为单位

        %% construct Hamiltonian @ valley K
        [Ham_m_LL_K, Ham_b_LL_K, D_K] = construct_trilayer_ABA_LL_six_bands(x0, x3, x4, gamma1, delta, gamma2, gamma5, Delta1, Delta2, +1, LL_index_cutoff, dims_m, dims_b);
        [Ham_m_LL_Kp, Ham_b_LL_Kp, D_Kp] = construct_trilayer_ABA_LL_six_bands(x0, x3, x4, gamma1, delta, gamma2, gamma5, Delta1, Delta2, -1, LL_index_cutoff, dims_m, dims_b);

        % helper_check_hermite(Ham_LL_K, 1e-8);

        %% 对哈密顿量进行对角化

        % call the eig sovler
        Ham_LL_K(1:dims_m, 1:dims_m) = Ham_m_LL_K;
        Ham_LL_K((dims_m + 1):end, (dims_m + 1):end) = Ham_b_LL_K;
        Ham_LL_K((dims_m + 1):end, 1:dims_m) = D_K;
        Ham_LL_K(1:dims_m, (dims_m + 1):end) = D_K';

        [eigvec_HK_now, eigval_HK] = eig(Ham_LL_K);
        eigval_HK_diag_now = diag(eigval_HK);

        Ham_LL_Kp(1:dims_m, 1:dims_m) = Ham_m_LL_Kp;
        Ham_LL_Kp((dims_m + 1):end, (dims_m + 1):end) = Ham_b_LL_Kp;
        Ham_LL_Kp((dims_m + 1):end, 1:dims_m) = D_Kp;
        Ham_LL_Kp(1:dims_m, (dims_m + 1):end) = D_Kp';

        [eigvec_HKp_now, eigval_HKp] = eig(Ham_LL_Kp);
        eigval_HKp_diag_now = diag(eigval_HKp);

        %% push into the LLs
        eigvals_LL_K(B_index, :) = eigval_HK_diag_now;
        eigvals_LL_Kp(B_index, :) = eigval_HKp_diag_now;

        %% 找出位于energy window [ene_lb, ene_ub]范围内的所有本征态和相应的本征能量
        [eigvec_HK_select_array, eigval_HK_select_list] = select_LLs_by_ene_window(eigvec_HK_now, eigval_HK_diag_now, ene_ub, ene_lb, dims);
        eig_info_HK_select_cell{B_index, 1} = eigvec_HK_select_array;
        eig_info_HK_select_cell{B_index, 2} = eigval_HK_select_list;
        eig_info_HK_select_cell{B_index, 3} = length(eigval_HK_select_list);

        [eigvec_HKp_select_array, eigval_HKp_select_list] = select_LLs_by_ene_window(eigvec_HKp_now, eigval_HKp_diag_now, ene_ub, ene_lb, dims_b);
        eig_info_HKp_select_cell{B_index, 1} = eigvec_HKp_select_array;
        eig_info_HKp_select_cell{B_index, 2} = eigval_HKp_select_list;
        eig_info_HKp_select_cell{B_index, 3} = length(eigval_HKp_select_list);

    end
end