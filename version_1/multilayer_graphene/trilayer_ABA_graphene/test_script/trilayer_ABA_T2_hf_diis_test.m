%% 添加路径
addpath('.\interaction_effects\');
addpath('..\hartree_fock_package\');
addpath('..\utils\')

%% 加载数据
B_field = 1.25;
LL_index_max = 12;
% 加载T2 states相关的信息
file_path = ['.\trilayer_ABA_data\T2_states_info_B_', num2str(roundn(B_field, -2)), 'T_LL_index_max_', num2str(LL_index_max), '.mat'];
load(file_path);
% 加载交换积分
file_path = ['.\trilayer_ABA_data\exchange_integrals_B_', num2str(roundn(B_field, -2)), 'T_LL_index_max_', num2str(LL_index_max), '.mat'];
load(file_path);

%% 基本参数
LL_index_max = 12;
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
one_electron_coulomb = 1.602176634 * 10^(-19); % 单位是C
epsilon_0 = 8.85 * 10^(-12) * 10^(-9); % 单位是F/nm
epsilon_bn = 5; % 无量纲，the in-plane dielectric constant of BN
d_interlayer = 0.335; % 单位是nm
a_intralayer = 0.246; % 单位是nm

%% energy scale
num_LL = LL_index_max + 1;
mag_length = 25.6 / sqrt(B_field); % 以nm为单位
E_zeeman = 0.116 * B_field / 1000; % 单位是eV
E_F = one_electron_coulomb / (4 * pi * epsilon_0 * mag_length); % 单位是 eV fock interaction strength
E_exchange = E_F / epsilon_bn; % 由于hBN存在介电常数（这里是in-plane dielectric constant），这里取为6.6，它可以有效减小交换相互作用的强度
E_H = d_interlayer / (2 * mag_length) * E_F; % 单位是 eV hartree interaction strength, 注意在ABC trilayer那里是(2 * d / l_b)

%% 其它矩阵
% 分别用1，2，3代表intralayer, nearest interlayer, next-nearest interlayer
dist_index_array = zeros(6,6);
for alpha = 1:6
    layer_alpha = mod(alpha, 2) + fix(alpha / 2);
    for beta = 1:6
        layer_beta = mod(beta, 2) + fix(beta / 2);
        if abs(layer_beta - layer_alpha) == 1
            dist_index_array(alpha, beta) = 2;
        elseif abs(layer_beta - layer_alpha) == 2
            dist_index_array(alpha, beta) = 3;
        else % 对应 alpha == beta
            dist_index_array(alpha, beta) = 1;
        end
    end
end

% 将本征态从phi_k basis转换到alpha basis
% 在eigvec_phi_k_component_Delta1_array中已经将本征态按照
% [|0,phi1>,|1,phi1>,|2,phi1>,|3,phi1>,...] // [|0,phi2>,|1,phi2>,|2,phi2>,|3,phi2>,...] // 
% [|0,phi3>,|1,phi3>,|2,phi3>,|3,phi3>,...] // [|0,phi4>,|1,phi4>,|2,phi4>,|3,phi4>,...] //
% [|0,phi5>,|1,phi5>,|2,phi5>,|3,phi5>,...] // [|0,phi6>,|1,phi6>,|2,phi6>,|3,phi6>,...] //

% transform_mat = zeros(num_LL, 6, 3); % 第二个指标是子格指标，第三个指标是本征态指标
transform_mat_cell = cell(num_LL, 1);
for n = 0:LL_index_max
    transform_mat_cell{n + 1} = T_mat * T2_eigvecs_cell{n + 1};
end

% 得到S_tensor, 这样就不用每次都计算exchange coeff，也不需要把密度矩阵切换到LL basis下
% <a|H_F|b> = - E_F * \sum_{c,d}rho(c,d) * S(a,d,c,b)
S_tensor = get_S_tensor_T2(exchange_integrals_list, transform_mat_cell, dist_index_array, indice_cell, LL_index_max);
% 厄密化
for a = 1:3
    for d = 1:3
        for c = 1:3
            for b = 1:3
                if abs(S_tensor(a, d, c, b) - conj(S_tensor(b, c, d, a))) >= 1e-16
                    disp("不厄密")
                    S_tensor(a, d, c, b)
                    conj(S_tensor(b, c, d, a))
                    S_tensor(a, d, c, b) - conj(S_tensor(b, c, d, a))
                end
                
                S_tensor(a, d, c, b) = (S_tensor(a, d, c, b) + conj(S_tensor(b, c, d, a))) / 2;
            end
        end
    end
end

% 得到R_mat用来得到hartree term
R_mat = get_R_mat_T2(transform_mat_cell, LL_index_max);
% R_mat
R_mat = refine_H_hf(R_mat, 3);
% R_mat

%% 随机构造一个密度矩阵(in the basis of triplet T2 states)
dim_subspace = 3; % ABA trilayer graphene 的 T2 states只有三个，并且都在K valley处
filling_factor = 1; % (整数)填充数，填充数主要影响的就是密度矩阵的构造
% density_matrix_down_temp = diag([1,0,0]);
% density_matrix_up_temp = zeros(3);
% 构造方法三
% density_matrix_down_old = construct_random_density_matrix(1, dim_subspace, "complex");
% density_matrix_up_old = zeros(3);
weight = rand();
density_matrix_down_old = weight * construct_random_density_matrix(1, dim_subspace, "complex");
density_matrix_up_old = (1 - weight) * construct_random_density_matrix(1, dim_subspace, "complex");

%% 初始化DIIS相关的矩阵(使用DIIS用以加速收敛，防止出现震荡)
num_diis = 6;
H_hf_down_diis_cell = cell(num_diis, 1);
H_hf_up_diis_cell = cell(num_diis, 1);
error_vector_cell = cell(num_diis, 1);
for ii = 1:num_diis
    %% 计算第二层上的电荷密度Delta_mid
    Delta_mid = get_delta_mid_method2(density_matrix_down_old, density_matrix_up_old, transform_mat_cell, LL_index_max);
    Delta_mid = real(Delta_mid);
    
    % H_hf_down_new = construct_ABA_trilayer_H_hf_T2_spin(density_matrix_down_old, S_tensor, R_mat, ...
    %                                                     T2_eigvals_list, E_zeeman, E_H, Delta_mid, E_exchange);
    
    [H_hf_up_new, H_hf_down_new] = construct_ABA_trilayer_H_hf_T2(density_matrix_up_old, density_matrix_down_old, S_tensor, R_mat, ...
                                                                  T2_eigvals_list, E_zeeman, E_H, Delta_mid, E_exchange);
    if ~helper_check_hermite(H_hf_down_new, 1e-16)
        disp("H_hf_down没有厄密性")
    end
    
    % 构造哈密顿量之后直接计算error_vector，然后再对H_hf进行对角化
    % error_vector = get_error_vector(H_hf_down_new, density_matrix_down_old, dim_subspace, "vector");
    error_vector_down = get_error_vector(H_hf_down_new, density_matrix_down_old, dim_subspace, "matrix");
    error_vector_up = get_error_vector(H_hf_up_new, density_matrix_up_old, dim_subspace, "matrix");
    
    error_vector_cell{mod(ii, num_diis) + 1} = [error_vector_down, zeros(dim_subspace);zeros(dim_subspace), error_vector_up];
    H_hf_down_diis_cell{mod(ii, num_diis) + 1} = H_hf_down_new;
    H_hf_up_diis_cell{mod(ii, num_diis) + 1} = H_hf_up_new;
    
    % 对角化自旋向下块
    [eigvecs_down_new, eigvals_down_new] = eig(H_hf_down_new); % MATLAB eig默认是从大到小进行排序
    [eigvecs_up_new, eigvals_up_new] = eig(H_hf_up_new); % MATLAB eig默认是从大到小进行排序
    eigvals_new = [real(diag(eigvals_down_new)), real(diag(eigvals_up_new))];
    [~, index_list] = sort(eigvals_new, 'ascend'); 
    
    for jj = 1:filling_factor
        if index_list(jj) > dim_subspace
            density_matrix_up_old = eigvecs_up_new(:, index_list(jj) - dim_subspace) * eigvecs_down_new(:, index_list(jj) - dim_subspace)';
        else
            density_matrix_down_old = eigvecs_down_new(:, index_list(jj)) * eigvecs_down_new(:, index_list(jj))'; 
        end
    end
    % density_matrix_down_old = eigvecs_down_new(:, index_list(1)) * eigvecs_down_new(:, index_list(1))' + eigvecs_down_new(:, index_list(2)) * eigvecs_down_new(:, index_list(2))';
end

%% get error
error_max = 0.0;
for ii = 1:num_diis
    error = max(max(abs(error_vector_cell{ii})));
    if error > error_max
        error_max = error;
    end
end

%% 开始DIIS iteration
num_steps = 1000;
num_record = 10;
steps = 0;
error_list = zeros(num_steps / num_record + 1, 1); % 存放误差

% 在初始阶段使用damping，效果可能会好一些
% damp_part_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.0];
damp_part = 0.3;

for ii = 1:num_steps
    %% 由density_matrix_old得到H_hf
    Delta_mid = get_delta_mid_method2(density_matrix_down_old, density_matrix_up_old, transform_mat_cell, LL_index_max);
    Delta_mid = real(Delta_mid);
    H_hf_down_new = construct_ABA_trilayer_H_hf_T2_spin(density_matrix_down_old, S_tensor, R_mat, ...
                                                        T2_eigvals_list, E_zeeman, E_H, Delta_mid, E_exchange);
                                                    
    %% 将H_hf存放到cell中保存
    H_hf_diis_cell{mod(ii, num_diis) + 1} = H_hf_down_new;
    
    %% 计算error vector并保存
    % error_vector = get_error_vector(H_hf_down_new, density_matrix_down_old, dim_subspace, "vector");
    error_vector = get_error_vector(H_hf_down_new, density_matrix_down_old, dim_subspace, "matrix");
    error_vector_cell{mod(ii, num_diis) + 1} = error_vector;
    error = max(max(abs(error_vector)));
    if error < error_max
        error_max = error;
    end
    
    if error_max < 1e-12 % 在error_max很小时退出迭代，否则会报错(线性方程组奇异性太大)
        disp("收敛结束")
        break
    end
    
    %% extrapolate H_hf
    H_hf_down_new = get_H_hf_by_diis(H_hf_diis_cell, error_vector_cell, num_diis, dim_subspace);
    
    %% 得到新的密度矩阵
    [eigvecs_down_new, eigvals_down_new] = eig(H_hf_down_new); % MATLAB eig默认是从大到小进行排序
    eigvals_new = [real(diag(eigvals_down_new))];
    [~, index_list] = sort(eigvals_new, 'ascend');
    
    density_matrix_down_old = eigvecs_down_new(:, index_list(1)) * eigvecs_down_new(:, index_list(1))';
    % density_matrix_down_old = eigvecs_down_new(:, index_list(1)) * eigvecs_down_new(:, index_list(1))' + eigvecs_down_new(:, index_list(2)) * eigvecs_down_new(:, index_list(2))';
end

% for ii = 1:num_steps
%     %% 确定系数c(j)
%     H_hf_down_new = get_H_hf_by_diis(H_hf_diis_cell, error_vector_cell, num_diis, dim_subspace);
%     
%     %% 得到error_vector
%     error_vector = get_error_vector(H_hf_down_new, density_matrix_down_old, dim_subspace);
%     error_vector_cell{mod(ii, num_diis) + 1} = error_vector;
%     
%     %% 对角化H_hf得到新的密度矩阵
%     [eigvecs_down_new, eigvals_down_new] = eig(H_hf_down_new); % MATLAB eig默认是从大到小进行排序
%     eigvals_new = [real(diag(eigvals_down_new))];
%     [~, index_list] = sort(eigvals_new, 'ascend');
%     density_matrix_down_new = eigvecs_down_new(:, index_list(1)) * eigvecs_down_new(:, index_list(1))';
%     
%     %% 得到误差
%     error = 0.0;
%     for jj = 1:num_diis
%         if max(abs(error_vector_cell{jj})) > error
%             error = max(abs(error_vector_cell{jj}));
%         end
%     end
%     
%     %% 计算第二层上的电荷密度Delta_mid
%     Delta_mid = get_delta_mid_method2(density_matrix_down_old, density_matrix_up_old, transform_mat_cell, LL_index_max);
%     Delta_mid = real(Delta_mid);
%     
%     %% 构造新的Hartree-Fock Hamiltonian来代替旧的Hartree-Fock Hamiltonian
%     H_hf_down = construct_ABA_trilayer_H_hf_T2_spin(density_matrix_down_new, S_tensor, R_mat, ...
%                                                     T2_eigvals_list, E_zeeman, E_H, Delta_mid, E_F);
%     if ~helper_check_hermite(H_hf_down, 1e-16)
%         disp("H_hf_down没有厄密性")
%     end
%     
%     H_hf_diis_cell{mod(ii, num_diis) + 1} = H_hf_down;
%     density_matrix_down_old = density_matrix_down_new;
%     
% %     steps = steps + 1;
% %     if steps == 1
% %         error_list(1) = error;
% %     end
% %     if ~mod(steps, num_record)
% %         error_list(steps / num_record + 1) = error;
% %         fprintf("进行到第%d步，误差为%10.36f\n", steps, error)
% %     end
% end

%% 作图
%% visualizing the symmetry broken states
[eigvecs_down_temp, eigvals_down_temp] = eig(H_hf_down_new);
eigvals_temp = [real(diag(eigvals_down_temp))];
[~, index_list] = sort(eigvals_temp, 'ascend');
eigvec_down_occ = eigvecs_down_temp(:, index_list(1));
eigvec_in_LL_basis = zeros(LL_index_max + 1, 6);
for n = 0:LL_index_max
    transform_mat = transform_mat_cell{n + 1};
    eigvec_in_LL_basis(n+1, :) = eigvec_down_occ(1) * transform_mat(:, 1) + eigvec_down_occ(2) * transform_mat(:, 2) + eigvec_down_occ(3) * transform_mat(:, 3);
end

dims_x = 201;
dims_y = 201;
x_list = linspace(-10, 10, dims_x); % 以mag_length为单位
y_list = linspace(-10, 10, dims_y); % 以mag_length为单位
[xx_mesh, yy_mesh] = meshgrid(x_list, y_list);
probability_density_mesh_after = zeros(dims_x, dims_y, 6);
for n = 0:LL_index_max
    % 得到波函数
    wave_packet_mesh = get_localized_wave_packet(xx_mesh, yy_mesh, dims_x, dims_y, n);
    % 得到probability density
    for alpha = 1:6
        probability_density_mesh_after(:, :, alpha) = probability_density_mesh_after(:, :, alpha) + eigvec_in_LL_basis(n+1, alpha) * wave_packet_mesh;
    end
end

probability_density_mesh_before_cell = cell(3,1);
for a = 1:3
    probability_density_mesh_before = zeros(dims_x, dims_y, 6);
    for n = 0:LL_index_max
        % 得到波函数
        wave_packet_mesh = get_localized_wave_packet(xx_mesh, yy_mesh, dims_x, dims_y, n);
        transform_mat = transform_mat_cell{n + 1};
        % 得到probability density
        for alpha = 1:6
            probability_density_mesh_before(:, :, alpha) = probability_density_mesh_before(:, :, alpha) + transform_mat(alpha, 1) * wave_packet_mesh;
        end
    end
    probability_density_mesh_before_cell{a} = probability_density_mesh_before;
end

probability_density_mesh_after = (abs(probability_density_mesh_after)).^2;
probability_density_mesh_after = sum(probability_density_mesh_after, 3);
for a = 1:3
    probability_density_mesh_before = probability_density_mesh_before_cell{a};
    probability_density_mesh_before = (abs(probability_density_mesh_before)).^2;
    probability_density_mesh_before = sum(probability_density_mesh_before, 3);
    probability_density_mesh_before_cell{a} = probability_density_mesh_before;
end

figure
im1 = imagesc(x_list, y_list, probability_density_mesh_after);
colorbar
figure
for a = 1:3
    subplot(1, 3, a)
    im1 = imagesc(x_list, y_list, probability_density_mesh_before_cell{a});
    colorbar
end

%% 保存数据
% filename = '_gully_symmetry_breaking.mat';
% save_path = ['D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\test_data\gully_symmetry\', datestr(datetime, 'yy-mm-dd-HH-MM-SS'), filename];
% save(save_path);