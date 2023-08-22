% 使用Hartree damping来进行自洽计算
%% 导入数据
% 从test_data / model_data 中导入相关的数据
% load
%% PRB 85，165139(2012)
zeroth_LL_info_cell = cell(6, 4); % 第一个指标存放名称(LL_K_m0 // LL_Kp_m0 // LL_K_b0 // LL_Kp_b0 // LL_K_b1 // LL_Kp_b1)
                                  % 第二个指标是 本征值(eigvals)
                                  % 第三个指标是 本征态各个component的分量(在基组phi_k下)
                                  % 第四个指标是 本征态各个component的分量(在基组alpha下)
if 1
    zeroth_LL_info_cell{1,1} = "LL_K_m0";
    zeroth_LL_info_cell{1,2} = 0.03; % 30meV = delta - gamma5 / 2
    zeroth_LL_info_cell{1,3} = [0, 1; 0, 0];
    zeroth_LL_info_cell{1,4} = get_LL_m_components_each_alpha(model.T_mat, zeroth_LL_info_cell{1,3}, 1);

    zeroth_LL_info_cell{2,1} = "LL_Kp_m0";
    zeroth_LL_info_cell{2,2} = 0.01; % 10mV = - gamma2 / 2
    zeroth_LL_info_cell{2,3} = [1, 0; 0, 0];
    zeroth_LL_info_cell{2,4} = get_LL_m_components_each_alpha(model.T_mat, zeroth_LL_info_cell{2,3}, 1);

    zeroth_LL_info_cell{3,1} = "LL_K_b0";
    zeroth_LL_info_cell{3,2} = 0.0; % 0meV = - 2 * Delta2
    zeroth_LL_info_cell{3,3} = [0, 0, 0, 1; 0, 0, 0, 0];
    zeroth_LL_info_cell{3,4} = get_LL_b_components_each_alpha(model.T_mat, zeroth_LL_info_cell{3,3}, 1);

    zeroth_LL_info_cell{4,1} = "LL_Kp_b0";
    zeroth_LL_info_cell{4,2} = -0.01; % -10mV = gamma2 / 2 + Delta2
    zeroth_LL_info_cell{4,3} = [1, 0, 0, 0; 0, 0, 0, 0];
    zeroth_LL_info_cell{4,4} = get_LL_b_components_each_alpha(model.T_mat, zeroth_LL_info_cell{4,3}, 1);

    zeroth_LL_info_cell{5,1} = "LL_K_b1";
    zeroth_LL_info_cell{5,2} = 0.0; % 0meV = - 2 * Delta2
    zeroth_LL_info_cell{5,3} = [0, 0, 0, 0; 0, 0, 0, 1];
    zeroth_LL_info_cell{5,4} = get_LL_b_components_each_alpha(model.T_mat, zeroth_LL_info_cell{5,3}, 1);

    zeroth_LL_info_cell{6,1} = "LL_Kp_b1";
    zeroth_LL_info_cell{6,2} = -0.01; % -10mV = gamma2 / 2 + Delta2
    zeroth_LL_info_cell{6,3} = [0, 0, 0, 0; 1, 0, 0, 0];
    zeroth_LL_info_cell{6,4} = get_LL_b_components_each_alpha(model.T_mat, zeroth_LL_info_cell{6,3}, 1);
end

dim_subspace = 6;
LL_index_max = 1; % 只需考虑到LL_index_max = 1即可

% LL_index_max为1，说明对于零朗道能级，即使考虑了gamma3的影响也只需要考虑到LL_index=1即可

% 得到transform_mat_cell
transform_mat_cell = cell(LL_index_max + 1, 1);
for n = 0:LL_index_max
    transform_mat = zeros(6, dim_subspace);
    for ii = 1:6
        transform_mat(:, ii) = zeroth_LL_info_cell{ii, 4}(n + 1, :, end);
    end
    transform_mat_cell{n + 1} = transform_mat;
end

%% 计算所需的exchange_integral
num_LL = LL_index_max + 1;
B_field = 10;
d_interlayer = 0.335; % 可调参数 0.335
d_interlayer_list = [0, d_interlayer, 2 * d_interlayer]; % 以nm为单位
kd_interlayer_list = 2 * d_interlayer_list / 25.6 * sqrt(B_field); % 无量纲值
epsilon_bn = 1.0; % 可调参数 6.6

% 开启多核
flag_parfor = false;
[exchange_integrals_list, indice_list, indice_cell] = calc_exchange_integrals(LL_index_max, kd_interlayer_list, flag_parfor);

% 得到energy scale
[E_zeeman, E_F, E_exchange, E_H] = model.get_energy_scales(epsilon_bn, B_field, d_interlayer);
% E_H = E_H * 2;
ene_unit = E_exchange * sqrt(pi / 2);
E_LL_list = zeros(1, dim_subspace);
for ii = 1:dim_subspace
    E_LL_list(ii) = zeroth_LL_info_cell{ii, 2}(end);
end

% 得到S_tensor, 这样就不用每次都计算exchange coeff，也不需要把密度矩阵切换到LL basis下
% <a|H_F|b> = - E_F * \sum_{c,d}rho(c,d) * S(a,d,c,b)
S_tensor = get_S_tensor(exchange_integrals_list, transform_mat_cell, model.dist_index_array, indice_cell, LL_index_max, dim_subspace);

% 得到R_mat用来得到hartree term
R_mat = get_R_mat(transform_mat_cell, LL_index_max, dim_subspace);

%% 对所有的filling factor进行迭代
ene_hf_array = zeros(2 * dim_subspace, 2 * dim_subspace);
for filling_factor = 1:2*dim_subspace % filling_factor = 12; % (整数)填充数，填充数主要影响的就是密度矩阵的构造
    iteration_max = 100;
    num_record = 100;
    steps = 0;
    error_list = zeros(iteration_max / num_record + 1, 1); % 存放误差
    
    % 在初始阶段使用damping，效果可能会好一些
    damping_factor = 0.3;
    ground_state_ene = inf; % 无限大
    for kk = 1:10 % 进行10次迭代，找出能量最低的态
        error_max = 1.0;
        % 随机构造一个密度矩阵(in the basis of zeroth LL)
        weight = rand();
        density_matrix_down_old = weight * construct_random_density_matrix(filling_factor, dim_subspace, "complex");
        density_matrix_up_old = (1 - weight) * construct_random_density_matrix(filling_factor, dim_subspace, "complex");
        % density_matrix_up_old = zeros(dim_subspace);
        
        while error_max > 1e-8
        % for ii = 1:iteration_max
            %% 由density_matrix_old得到H_hf
            Delta_mid = get_delta_mid_method2(density_matrix_down_old, density_matrix_up_old, transform_mat_cell, LL_index_max);
            Delta_mid = real(Delta_mid);

            H_hf_down_new = construct_ABA_trilayer_H_hf_spin(density_matrix_down_old, S_tensor, R_mat, ...
                                                             E_LL_list, E_exchange, E_H, Delta_mid, E_zeeman, -1, dim_subspace);
            H_hf_up_new = construct_ABA_trilayer_H_hf_spin(density_matrix_up_old, S_tensor, R_mat, ...
                                                           E_LL_list, E_exchange, E_H, Delta_mid, E_zeeman, +1, dim_subspace);
            if ~helper_check_hermite(H_hf_down_new, 1e-16)
                disp("H_hf_down没有厄密性")
            end

    %         if error_max < 1e-16 % 在error_max很小时退出迭代，否则会报错(线性方程组奇异性太大)
    %             disp("收敛提前结束")
    %             break
    %         end

            %% 得到新的密度矩阵
            [eigvecs_down_new, eigvals_down_new] = eig(H_hf_down_new); % MATLAB eig默认是从大到小进行排序
            [eigvecs_up_new, eigvals_up_new] = eig(H_hf_up_new); % MATLAB eig默认是从大到小进行排序
            eigvals_new = [real(diag(eigvals_down_new)); real(diag(eigvals_up_new))];
            [~, index_list] = sort(eigvals_new, 'ascend'); 
            
            density_matrix_down_new = zeros(dim_subspace);
            density_matrix_up_new = zeros(dim_subspace);
            for jj = 1:filling_factor
                if index_list(jj) > dim_subspace
                    density_matrix_up_new = density_matrix_up_new + eigvecs_up_new(:, index_list(jj) - dim_subspace) * eigvecs_up_new(:, index_list(jj) - dim_subspace)';
                else
                    density_matrix_down_new =  density_matrix_down_new + eigvecs_down_new(:, index_list(jj)) * eigvecs_down_new(:, index_list(jj))'; 
                end
            end

            %% 计算误差
            error_matrix = abs(density_matrix_up_new - density_matrix_up_old) + abs(density_matrix_down_new - density_matrix_down_old);
            error_max = max(max(error_matrix));

            %% 为下一步做准备
            density_matrix_up_old = density_matrix_up_new * damping_factor + density_matrix_up_old * (1 - damping_factor);
            density_matrix_down_old = density_matrix_down_new * damping_factor + density_matrix_down_old * (1 - damping_factor);

            steps = steps + 1;
            if mod(steps, iteration_max) == 0
                disp(["误差为", error_max])
            end

        end
        
        disp(["填充数", filling_factor, "已经收敛"])
        if sum(eigvals_new(index_list(1:filling_factor))) < ground_state_ene
            ground_state_ene = sum(eigvals_new(index_list(1:filling_factor)));
            % 保存能量到ene_hf_array
            ene_hf_array(:, filling_factor) = eigvals_new(index_list) / ene_unit;
        end
    end
end


figure
poly_width = 0.2;
poly_height = 0.005;
hold on
for f_factor = 1:(2 * dim_subspace)
    x_list = [f_factor - poly_width  f_factor - poly_width  f_factor + poly_width  f_factor + poly_width];
    for ii = 1:(2 * dim_subspace)
        ene_hf = ene_hf_array(ii, f_factor);
        y_list = [ene_hf - poly_height  ene_hf + poly_height  ene_hf + poly_height  ene_hf - poly_height];
        pgon = polyshape(x_list, y_list);
        plot(pgon)
    end
end
f_factor = 0;
x_list = [f_factor - poly_width  f_factor - poly_width  f_factor + poly_width  f_factor + poly_width];
for ii = 1:(2 * dim_subspace)
    ene_hf = E_LL_list(mod(ii, 2) + fix(ii / 2)) / ene_unit;
    y_list = [ene_hf - poly_height  ene_hf + poly_height  ene_hf + poly_height  ene_hf - poly_height];
    % pgon = polyshape(x_list, y_list);
    % plot(pgon)
    fill(x_list, y_list, 'r')
end
