function plot_LLs_with_config()
    % 单层双层全部画(颜色一样)
    % eigvals_LL_cell = cell(2, 4);
    % eigvals_LL_cell{1,1} = eigvals_LL_K;
    % eigvals_LL_cell{1,2} = dims; % 维数 dims : 单层双层全部画
    % eigvals_LL_cell{1,3} = 'b--'; % line color
    % eigvals_LL_cell{1,4} = 0.5; % line width
    % 
    % eigvals_LL_cell{2,1} = eigvals_LL_Kp;
    % eigvals_LL_cell{2,2} = dims;
    % eigvals_LL_cell{2,3} = 'b';
    % eigvals_LL_cell{2,4} = 1.0;

    % 单层双层全部画(颜色不同)
    eigvals_LL_cell = cell(4, 4); % mono_K, bi_K, mono_Kp, bi_Kp
    eigvals_LL_cell{1,1} = eigvals_LL_K(:, 1:dims_m);
    eigvals_LL_cell{1,2} = dims_m; % 维数 dims : 单层双层全部画
    eigvals_LL_cell{1,3} = 'b'; % line color
    eigvals_LL_cell{1,4} = 0.5; % line width

    eigvals_LL_cell{2,1} = eigvals_LL_Kp(:, 1:dims_m);
    eigvals_LL_cell{2,2} = dims_m; % 维数 dims : 单层双层全部画
    eigvals_LL_cell{2,3} = 'b--'; % line color
    eigvals_LL_cell{2,4} = 0.5; % line width

    eigvals_LL_cell{3,1} = eigvals_LL_K(:, (dims_m + 1):dims);
    eigvals_LL_cell{3,2} = dims_b;
    eigvals_LL_cell{3,3} = 'r';
    eigvals_LL_cell{3,4} = 0.5;

    eigvals_LL_cell{4,1} = eigvals_LL_Kp(:, (dims_m + 1):dims);
    eigvals_LL_cell{4,2} = dims_b;
    eigvals_LL_cell{4,3} = 'r--';
    eigvals_LL_cell{4,4} = 0.5;

    % % 只画单层
    % eigvals_LL_cell = cell(2, 4);
    % eigvals_LL_cell{1,1} = eigvals_LL_K(:, 1:dims_m); % eigvals_LL_K(B_index, 1:dims_m) = eigval_HK_m_diag_now;
    % eigvals_LL_cell{1,2} = dims_m;
    % eigvals_LL_cell{1,3} = 'b--'; % line color
    % eigvals_LL_cell{1,4} = 0.5; % line width
    % 
    % eigvals_LL_cell{2,1} = eigvals_LL_Kp(:, 1:dims_m);
    % eigvals_LL_cell{2,2} = dims_m;
    % eigvals_LL_cell{2,3} = 'b';
    % eigvals_LL_cell{2,4} = 1.0;

    % % 只画双层
    % eigvals_LL_cell = cell(1, 4);
    % eigvals_LL_cell{1,1} = eigvals_LL_K(:, (dims_m + 1):dims); 
    % eigvals_LL_cell{1,2} = dims_b;
    % eigvals_LL_cell{1,3} = 'b--'; % line color
    % eigvals_LL_cell{1,4} = 0.5; % line width
    % 
    % eigvals_LL_cell{2,1} = eigvals_LL_Kp(:, (dims_m + 1):dims);
    % eigvals_LL_cell{2,2} = dims_b;
    % eigvals_LL_cell{2,3} = 'b';
    % eigvals_LL_cell{2,4} = 1.0;
    % 
    save_path = "";
    E_bottom = -0.100;
    E_top = 0.100;
    fig0 = plot_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eigvals_LL_cell, save_path);
    % fig1 = plot_LLs_weights(eigvec_HK_now(:,2 * LL_index_cutoff - 3 :2 * LL_index_cutoff + 3)); % 对本征态在LL basis下的权重作图
    % fig2 = plot_select_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eig_info_HK_b_select_cell, "r", save_path);
    % fig3 = plot_select_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eig_info_HKp_b_select_cell, "r", save_path);
    % fig4 = plot_select_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eig_info_HK_m_select_cell, "b", save_path);
    % fig5 = plot_select_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eig_info_HKp_m_select_cell, "b", save_path);
end