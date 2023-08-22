function trilayer_ABA_LLs_EB_function(hopping_params, B_fields_list, LL_index_cutoff)
    %% parameters set up
    %% 基本参数
    h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
    one_electron_coulomb = 1.602176634 * 10^(-19); % 单位是C
    epsilon_0 = 8.85 * 10^(-12); % 单位是F/m
    epsilon_bn = 6.6;
    d_interlayer = 0.335; % 单位是nm
    a_intralayer = 0.246; % 单位是nm
    
    gamma0 = hopping_params(1);
    gamma1 = hopping_params(2);
    gamma2 = hopping_params(3);
    gamma3 = hopping_params(4);  
    gamma4 = hopping_params(5);
    gamma5 = hopping_params(6);
    
    delta = hopping_params(7);
    Delta2 = hopping_params(8);
    Delta1 = 0.0; 
    
    % 将gamma_i转换为v_i，同时某些gamma_i需要扩大sqrt(2)倍!!!!
    gamma1 = gamma1 * sqrt(2);
    v0 = sqrt(3) / 2 * gamma0 * (a_intralayer * 10^(-9)) / h_bar;
    v3 = sqrt(3) / 2 * sqrt(2) * gamma3 * (a_intralayer * 10^(-9)) / h_bar; % 不要忘了乘以sqrt(2) !!!
    v4 = sqrt(3) / 2 * sqrt(2) * gamma4 * (a_intralayer * 10^(-9)) / h_bar; % 不要忘了乘以sqrt(2) !!!
    
    B_start = B_fields_list(1);
    B_end = B_fields_list(end);
    B_steps = length(B_fields_list);
    % B_fields_list = linspace(B_start, B_end, B_steps);
    
    %% 第一部分：Landau level as a function of E & B
    %% 朗道能级参数
    % LL_index_cutoff = 30;
    dims_m = 2 * LL_index_cutoff + 1;
    dims_b = 4 * LL_index_cutoff;
    dims = dims_b + dims_m;
    % energy window
    ene_ub = 0.1; % 100meV
    ene_lb = - 0.1; % -100meV
    
    %% 计算LL fan diagram
    tic
    [eigvals_LL_K, eigvals_LL_Kp, eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, eig_info_HK_m_select_cell, eig_info_HKp_m_select_cell] = ...
                trilayer_ABA_LL_solver_nB_without_Delta1_with_ene_win(v0, v3, v4, gamma1, delta, gamma2, gamma5, Delta2, B_fields_list, B_steps, LL_index_cutoff, dims_m, dims_b, ene_ub, ene_lb);
    toc
    
    %% 作图
    flag_plot_LL = true;
    if flag_plot_LL
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

    %% 对selected_LL eig info进行重新排列（从大磁场开始）
    % ene of LLL
    [LL_K_m_0, LL_Kp_m_0, LL_K_b_0, LL_Kp_b_0] = trilayer_LLL_ene(gamma2, gamma5, delta, Delta2);
    LL_K_b_1 = LL_K_b_0;
    LL_Kp_b_1 = LL_Kp_b_0;

    LL_K_m0_index_cell = cell(B_steps, 1);
    LL_Kp_m0_index_cell = cell(B_steps, 1);
    LL_K_b0_index_cell = cell(B_steps, 1);
    LL_K_b1_index_cell = cell(B_steps, 1);
    LL_Kp_b0_index_cell = cell(B_steps, 1);
    LL_Kp_b1_index_cell = cell(B_steps, 1);

    % 如果不考虑gamma3的影响，bilayer-like branch的0th LL和monolayer-like branch的0th LL都是有解析形式的并且是完全极化的
    weight = 0.8; % 如果不考虑gamma3的影响，那么weight应该是1，因此我们只需要关心gamma3的大小即可，即weight实际上应该是gamma3的函数, weight(gamma3 = 0) = 1
    ene_eps = 0.005; % 5meV的误差

    % 先确定单层对应的LL_K_m0, LL_Kp_m0, 
    [eig_info_HK_m_select_cell, eig_info_HKp_m_select_cell] = sort_LLs_m(eig_info_HK_m_select_cell, eig_info_HKp_m_select_cell, ...
                                                                         LL_K_m_0, LL_Kp_m_0, LL_index_cutoff, B_steps, ene_eps, weight);

    % 寻找LL_K_b0和LL_K_b1的指标
    weight_K_b0 = weight;
    weight_Kp_b0 = weight;
    flag_slope_K_b0 = false;
    flag_slope_Kp_b0 = false;
    [eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, eigval_LL_K_b0_list, eigval_LL_Kp_b0_list] = ...
                        sort_LLs_b(eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, LL_K_b_0, LL_Kp_b_0, LL_index_cutoff, B_steps, ene_eps, weight_K_b0, weight_Kp_b0);

    %% 验证是否找对了
    figure
    hold on
    plot(B_fields_list, eigval_LL_Kp_b0_list)
    plot(B_fields_list, eigval_LL_K_b0_list)

    %% 确定了所有LL的index之后可以确定交叉点(crossing point)了
    % 确定 LL_K_b2, LL_K_b3, LL_K_b4, LL_Kp_b2, LL_Kp_b3, LL_Kp_b4 与 LL_K_m0, LL_Kp_m0 的交叉点的位置
    % 其中 LL_K_m0, LL_Kp_m0 是严格的直线，因此在考虑disorder之后就是一个矩形
    ene_width_disorder = 0.002; % 假设有2meV的disorder带来的能量展开

    [LL_K_m0_crossing_points_list, LL_Kp_m0_crossing_points_list, poly_LL_K_m0, poly_LL_Kp_m0, eigval_LL_K_b_list, eigval_LL_Kp_b_list] = ...
                trilayer_ABA_LLs_find_cross_points(eig_info_HK_b_select_cell, eig_info_HKp_b_select_cell, LL_K_m_0, LL_Kp_m_0, ene_width_disorder, B_steps, B_fields_list);
    figure
    hold on
    plot(poly_LL_K_m0)
    plot(poly_LL_Kp_m0)
    for ii = 1:4
        plot(B_fields_list, eigval_LL_K_b_list(ii, :), 'r');
        plot(B_fields_list, eigval_LL_Kp_b_list(ii, :), 'r--');

        plot(LL_K_m0_crossing_points_list(1, ii), LL_K_m_0 - ene_width_disorder, 'g*') % left crossing points
        plot(LL_K_m0_crossing_points_list(2, ii), LL_K_m_0 + ene_width_disorder, 'g*') % left crossing points

        plot(LL_Kp_m0_crossing_points_list(1, ii), LL_Kp_m_0 - ene_width_disorder, 'b*') % left crossing points
        plot(LL_Kp_m0_crossing_points_list(2, ii), LL_Kp_m_0 + ene_width_disorder, 'b*') % left crossing points
    end
end