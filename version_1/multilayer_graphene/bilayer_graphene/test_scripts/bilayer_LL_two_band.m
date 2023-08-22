%% test
% v = 0.8 * 10^6; % 来自gamma0，单位是m/s
% gamma1 = 0.39;
% % gamma3 = 0.9; % 0.0 / 0.1 * gamma0
% % v3 = 0.1 * v; % 0.0 * v / 0.1 * v
% v3_steps = 300;
% v3_list = linspace(0, 0.3, v3_steps);
% v4 = 0;
% delta = 0;
% u = 0.0; % 0.0 / 0.04 / 0.065 / 0.07 / 0.08
% h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
% single_electron = 1.6021892 * 10^(-19); % 以C为单位
% 
% LL_index_cutoff = 100;
% dims = 2 * LL_index_cutoff;
% eigvals_LL_K = zeros(v3_steps, dims);
% eigvals_LL_Kp = zeros(v3_steps, dims);
% 
% B_field = 0.1;
% mag_length = 25.66 / sqrt(B_field); % 以nm为单位
% x0 = h_bar * v / mag_length * 10^9;
% x4 = h_bar * v4 / mag_length * 10^9;
% ene_cyclo = 2 * x0 ^ 2 / gamma1;
% %%
% tic
% for v3_index = 1:v3_steps
%     %% parameters set up
%     % v3 = v3_list(v3_index) * v;
%     x3 = x0 * v3_list(v3_index);
%     
%     %% construct Hamiltonian @ valley K
%     Ham_LL_K = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, +1, LL_index_cutoff, dims);
%     Ham_LL_Kp = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, -1, LL_index_cutoff, dims);
%     
%     % helper_check_hermite(Ham_LL_K, 1e-8);
%     
%     %% 对哈密顿量进行对角化
%     % call the eig sovler
%     [eigvec_HK_now, eigval_HK] = eig(Ham_LL_K);
%     eigval_HK_diag_now = diag(eigval_HK);
%     
%     [eigvec_HKp_now, eigval_HKp] = eig(Ham_LL_Kp);
%     eigval_HKp_diag_now = diag(eigval_HKp);
%     
%     % [eigvec_HK_now, eigval_HK_diag_now, eig_num_K] = helper_re_order_states(eigvec_HK_last, eigval_HK_diag_last, ...
%     %                                                                         eigvec_HK_now, eigval_HK_diag_now, dims, eps);
% 
%     % push into the LLs
%     eigvals_LL_K(v3_index, :) = eigval_HK_diag_now;
%     eigvals_LL_Kp(v3_index, :) = eigval_HKp_diag_now;
%     
% end
% toc
% 
% %% 进行作图
% figure
% % subplot(2,1,1)
% axis([v3_list(1) v3_list(end) -25 25])
% hold on
% 
% % plot bilayer-like LL
% for i = 1:dims
%     % plot(v3_list, eigvals_LL_K(:,i), 'b', 'LineWidth', 0.5)
%     plot(v3_list, eigvals_LL_Kp(:,i) / ene_cyclo, 'b--', 'LineWidth', 1.0)
% end

%% parameters set up
v = 1 * 10^6; % 来自gamma0，单位是m/s
gamma1 = 0.40;
% gamma3 = 0.1 * gamma0; % 0.0 / 0.1 * gamma0
v3 = 0.1 * v; % 0.0 * v / 0.1 * v
v4 = 0;
delta = 0;
u = 0.2 * gamma1; % 0.0 / 0.04 / 0.065 / 0.07 / 0.08
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位

B_start = 0.0;
B_end = 6;
B_steps = 601;
B_fields_list = linspace(B_start, B_end, B_steps);

LL_index_cutoff = 100;
dims = 2 * LL_index_cutoff;
eigvals_LL_K = zeros(B_steps, dims);
eigvals_LL_Kp = zeros(B_steps, dims);

tic
for B_index = 1:B_steps
    %% parameters set up
    B_field = B_fields_list(B_index);
    mag_length = 25.66 / sqrt(B_field); % 以nm为单位
    x0 = h_bar * v / mag_length * 10^9;
    x3 = h_bar * v3 / mag_length * 10^9;
    x4 = h_bar * v4 / mag_length * 10^9;
    
    %% construct Hamiltonian @ valley K
    Ham_LL_K = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, +1, LL_index_cutoff, dims);
    Ham_LL_Kp = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, -1, LL_index_cutoff, dims);
    
    % helper_check_hermite(Ham_LL_K, 1e-8);
    
    %% 对哈密顿量进行对角化
    % call the eig sovler
    [eigvec_HK_now, eigval_HK] = eig(Ham_LL_K);
    eigval_HK_diag_now = diag(eigval_HK);
    
    [eigvec_HKp_now, eigval_HKp] = eig(Ham_LL_Kp);
    eigval_HKp_diag_now = diag(eigval_HKp);
    
    % [eigvec_HK_now, eigval_HK_diag_now, eig_num_K] = helper_re_order_states(eigvec_HK_last, eigval_HK_diag_last, ...
    %                                                                         eigvec_HK_now, eigval_HK_diag_now, dims, eps);

    % push into the LLs
    eigvals_LL_K(B_index, :) = eigval_HK_diag_now;
    eigvals_LL_Kp(B_index, :) = eigval_HKp_diag_now;
    
end
toc

%% 进行作图
addpath(".\plot_funcs\")
eigvals_LL_cell = cell(2, 4);
eigvals_LL_cell{1,1} = eigvals_LL_K / gamma1;
eigvals_LL_cell{1,2} = dims; % 维数
eigvals_LL_cell{1,3} = 'b'; % line color
eigvals_LL_cell{1,4} = 0.5; % line width

eigvals_LL_cell{2,1} = eigvals_LL_Kp / gamma1;
eigvals_LL_cell{2,2} = dims;
eigvals_LL_cell{2,3} = 'b--';
eigvals_LL_cell{2,4} = 1.0;

save_path = ['D:\matlab\graphene-package\multilayer_graphene\bilayer_graphene\bilayer_figs\bilayer_LLs_v3_', num2str(roundn(v3 / v, -2)), '.jpg'];
current_fig = plot_LLs(B_start, B_end, -0.125, 0.125, B_fields_list, eigvals_LL_cell, save_path);

% save_path = ['D:\matlab\graphene-package\multilayer_graphene\bilayer_graphene\bilayer_figs\bilayer_LLs_Kp_valence_band_v3_', num2str(roundn(v3 / v, -2)), '.jpg'];
% current_fig = plot_LLs(B_start, B_end, -0.125, -0.075, B_fields_list, eigvals_LL_cell, save_path);
% save_path = ['D:\matlab\graphene-package\multilayer_graphene\bilayer_graphene\bilayer_figs\bilayer_LLs_Kp_conduction_band_v3_', num2str(roundn(v3 / v, -2)), '.jpg'];
% current_fig = plot_LLs(B_start, B_end, 0.075, 0.125, B_fields_list, eigvals_LL_cell, save_path);

%% 选取最靠近CNP的6条朗道能级(3条位于导带，3条位于价带)
% K valley 取价带能量绝对值最小的三个, Kp valley 取导带能量绝对值最小的三个
negative_index_list = find(eigvals_LL_K(end, :) < 0);
negative_ene_vals = eigvals_LL_K(end, eigvals_LL_K(end, :) < 0);
[~, sort_index_list] = sort(negative_ene_vals, 'descend');
eigvals_LLL_K = eigvals_LL_K(:, negative_index_list(sort_index_list(1:3)));

positive_index_list = find(eigvals_LL_Kp(end, :) > 0);
positive_ene_vals = eigvals_LL_Kp(end, eigvals_LL_Kp(end, :) > 0);
[~, sort_index_list] = sort(positive_ene_vals, 'ascend');
eigvals_LLL_Kp = eigvals_LL_Kp(:, positive_index_list(sort_index_list(1:3)));

eigvals_LLL_cell = cell(2, 4);
eigvals_LLL_cell{1,1} = eigvals_LLL_K / gamma1;
eigvals_LLL_cell{1,2} = 3; % 维数
eigvals_LLL_cell{1,3} = 'b'; % line color
eigvals_LLL_cell{1,4} = 0.5; % line width

eigvals_LLL_cell{2,1} = eigvals_LLL_Kp / gamma1;
eigvals_LLL_cell{2,2} = 3;
eigvals_LLL_cell{2,3} = 'b--';
eigvals_LLL_cell{2,4} = 1.0;

save_path = ['D:\matlab\graphene-package\multilayer_graphene\bilayer_graphene\bilayer_figs\bilayer_LLLs_cnp_v3_', num2str(roundn(v3 / v, -2)), '.jpg'];
current_fig = plot_LLs(B_start, B_end, -0.125, 0.125, B_fields_list, eigvals_LLL_cell, save_path);

