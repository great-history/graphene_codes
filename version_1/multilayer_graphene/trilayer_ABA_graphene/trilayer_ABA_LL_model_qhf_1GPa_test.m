%% 该脚本用来模拟1GPa下时zeroth LL的splitting

%% 导入数据
% 从test_data / model_data 中导入相关的数据
% load
%% 找出monolayer和bilayer的zeroth LL
% 找出在较高磁场(B=4.5T/10T/14T)下的LL_K_m_0 // LL_Kp_m_0在相互作用下导致的quantum hall ferromagnetism
model.flag_EB = true;
model.flag_ED = false;
field1 = 5.0;
field_index1 = model.find_closet_field_index(field1);
field2 = 10.0;
field_index2 = model.find_closet_field_index(field2);
field3 = 14.0;
field_index3 = model.find_closet_field_index(field3);

% zeroth_LL_info_cell = model.get_zeroth_LL_info_cell(field_index2 - 9, field_index2);
LL_K_m_index_bound = [0, 0];
LL_Kp_m_index_bound = [0, 0];
LL_K_b_index_bound = [-1, 2];
LL_Kp_b_index_bound = [-1, 2];
num_LL_K_m  = LL_K_m_index_bound(2) - LL_K_m_index_bound(1) + 1;
num_LL_Kp_m = LL_Kp_m_index_bound(2) - LL_Kp_m_index_bound(1) + 1;
num_LL_K_b  = LL_K_b_index_bound(2) - LL_K_b_index_bound(1) + 1;
num_LL_Kp_b = LL_Kp_b_index_bound(2) - LL_Kp_b_index_bound(1) + 1;
dim_subspace = num_LL_K_m + num_LL_Kp_m + num_LL_K_b + num_LL_Kp_b;
zeroth_LL_info_cell = model.get_lowest_LL_info_cell(field_index2 - 9, field_index2, LL_K_m_index_bound, LL_Kp_m_index_bound, LL_K_b_index_bound, LL_Kp_b_index_bound);

%% 确定出LL_index_max
threshold = 0.90; % 只有权重大于threshold的LL_index才会被记录
LL_index_max = 0; % 单层零朗道能级完全极化为|0>
for ii = 3:6
    norm_sum = 0.0;
    for n = 0:model.LL_index_cutoff
        eigvec_alpha_component = zeroth_LL_info_cell{ii, 4}(n + 1, :, end);
        norm_sum = norm_sum + sum(abs(eigvec_alpha_component).^2);
        if sqrt(norm_sum) >= threshold
            if n > LL_index_max
                LL_index_max = n;
            end
            
            break
        end
    end
end

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
B_field = field2;
d_interlayer = 0.300; % 可调参数 0.335
d_interlayer_list = [0, d_interlayer, 2 * d_interlayer]; % 以nm为单位
kd_interlayer_list = 2 * d_interlayer_list / 25.6 * sqrt(B_field); % 无量纲值
epsilon_bn = 6.6; % 可调参数 5.0 / 6.6 / 7.0 / 8.6 / 10.0

[exchange_integrals_list, indice_list, indice_cell] = calc_exchange_integrals(LL_index_max, kd_interlayer_list, flag_parfor);

% 得到energy scale
[E_zeeman, E_F, E_exchange, E_H] = model.get_energy_scales(epsilon_bn, B_field, d_interlayer);
E_H = E_H * 2;
% E_H = 0 * E_H;
ene_unit = E_exchange * sqrt(pi / 2);
E_LL_list = zeros(1, dim_subspace);
for ii = 1:dim_subspace
    E_LL_list(ii) = zeroth_LL_info_cell{ii, 2}(end);
end

S_tensor = get_S_tensor(exchange_integrals_list, transform_mat_cell, model.dist_index_array, indice_cell, LL_index_max, dim_subspace);
R_mat = get_R_mat(transform_mat_cell, LL_index_max, dim_subspace);

%% 对所有的filling factor进行迭代
ene_hf_array = zeros(2 * dim_subspace, 2 * dim_subspace);
occu_ene_array = zeros(2 * dim_subspace, 2 * dim_subspace);

% for filling_factor = 1:2*dim_subspace % filling_factor = 12; % (整数)填充数，填充数主要影响的就是密度矩阵的构造
filling_start = 13; % 13   1
filling_end = 16;   % 16   4
for filling_factor = filling_start:filling_end % filling_factor = 12; % (整数)填充数，填充数主要影响的就是密度矩阵的构造
    iteration_max = 10;
    num_record = 100;
    steps = 0;
    % error_list = zeros(iteration_max / num_record + 1, 1); % 存放误差
    
    % 在初始阶段使用damping，效果可能会好一些
    damping_factor = 0.3;
    ground_state_ene = inf; % 无限大
    for kk = 1:10 % 进行10次迭代，找出能量最低的态
        error_max = 1.0;
        % 随机构造一个密度矩阵(in the basis of zeroth LL)
        weight = rand();
        density_matrix_down_old = weight * construct_random_density_matrix(filling_factor, dim_subspace, "complex");
        density_matrix_up_old = (1 - weight) * construct_random_density_matrix(filling_factor, dim_subspace, "complex");
        
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
            % ene_hf_array(:, filling_factor) = eigvals_new(index_list);
            
            for jj = 1:filling_factor
                if index_list(jj) > dim_subspace
                    occu_ene_array(jj, filling_factor) = + 1; % 代表自旋向上
                else
                    occu_ene_array(jj, filling_factor) = - 1; % 代表自旋向上
                end
            end
        end
        
    end
end

%% 作图
% 按照自旋向上/向下来区分
figure
title(['dim subspace = ', num2str(dim_subspace), ', epsilon bn = ', num2str(epsilon_bn)])
poly_width = 0.2;
poly_height = 0.005;
hold on
cnp_f_factor = 1 - LL_Kp_b_index_bound(1) + 1 + 1 - LL_K_b_index_bound(1);
cnp_f_factor = cnp_f_factor * 2;
% for f_factor = 1 - cnp_f_factor:(2 * dim_subspace) - cnp_f_factor % f_factor is the real filling factor with regard to the cnp
for f_factor = filling_start:filling_end % f_factor is the real filling factor with regard to the cnp
    x_list = [f_factor - poly_width  f_factor - poly_width  f_factor + poly_width  f_factor + poly_width];
    for ii = 1:(2 * dim_subspace)
        ene_hf = ene_hf_array(ii, f_factor);
        y_list = [ene_hf - poly_height  ene_hf + poly_height  ene_hf + poly_height  ene_hf - poly_height];
        % pgon = polyshape(x_list, y_list);
        % plot(pgon)
        if ii <= f_factor
            if occu_ene_array(ii, f_factor) == 1 % spin up
                p = fill(x_list, y_list, 'r');
                p.EdgeColor = 'none';  
            elseif occu_ene_array(ii, f_factor) == -1 % spin down
                p = fill(x_list, y_list, 'b');
                p.EdgeColor = 'none';  
            end
              
        else
            p = fill(x_list, y_list, 'g');
            p.EdgeColor = 'none';    
        end
        
    end
end

% f_factor = - cnp_f_factor;
% x_list = [f_factor - poly_width  f_factor - poly_width  f_factor + poly_width  f_factor + poly_width];
% for ii = 1:(2 * dim_subspace)
%     ene_hf = E_LL_list(mod(ii, 2) + fix(ii / 2)) / ene_unit;
%     % ene_hf = E_LL_list(mod(ii, 2) + fix(ii / 2));
%     y_list = [ene_hf - poly_height  ene_hf + poly_height  ene_hf + poly_height  ene_hf - poly_height];
%     % pgon = polyshape(x_list, y_list);
%     % plot(pgon)
%     p = fill(x_list, y_list, 'b');
%     p.EdgeColor = 'none';
% end