addpath('.\classes\')
% 创建模型
% hopping_params = output_value_array(1,:);
gamma0 = 3100;
gamma1 = 390;
gamma2 = -28;
gamma3 = 315;  
gamma4 = 41;
gamma5 = 50;
delta = 46;

Delta1 = 0; % 0 / 25 / 50 / 150 / 180 / 250
Delta2 = 0;

hopping_params(1) = gamma0;
hopping_params(2) = gamma1;
hopping_params(3) = gamma2;
hopping_params(4) = gamma3;  
hopping_params(5) = gamma4;
hopping_params(6) = gamma5;
hopping_params(7) = delta;
hopping_params(8) = Delta2;

model = trilayer_ABA_class(hopping_params);

% 计算LL EB
Delta1 = 0; 
B_start = 1.0;
B_end = 10;
B_steps = 451;
B_fields_list = linspace(B_start, B_end, B_steps);
LL_index_cutoff = 30;
ene_ub = 0.1;
ene_lb = -0.1;
ene_eps = 0.002;
model = model.trilayer_ABA_LLs_EB(Delta1, B_start, B_end, B_steps, LL_index_cutoff, ene_ub, ene_lb, ene_eps);

% 作图LL EB
line_info_cell = cell(4,2);
line_info_cell{1,1} = 'b';
line_info_cell{1,2} = 0.5;
line_info_cell{2,1} = 'r--';
line_info_cell{2,2} = 1.0;
line_info_cell{3,1} = 'b--';
line_info_cell{3,2} = 0.5;
line_info_cell{4,1} = 'b';
line_info_cell{4,2} = 1.0;
save_path = "";
fig0 = model.representative_lines_plot(line_info_cell, ene_lb, ene_ub, save_path);

% 如果不考虑gamma3的影响，bilayer-like branch的0th LL和monolayer-like branch的0th LL都是有解析形式的并且是完全极化的
weight = 0.8; % 如果不考虑gamma3的影响，那么weight应该是1，因此我们只需要关心gamma3的大小即可，即weight实际上应该是gamma3的函数, weight(gamma3 = 0) = 1
ene_eps = 0.005; % 5meV的误差
% ene of LLL
[LL_K_m_0, LL_Kp_m_0, LL_K_b_0, LL_Kp_b_0] = trilayer_LLL_ene(hopping_params(3), hopping_params(6), hopping_params(7), hopping_params(8));
LL_K_b_1 = LL_K_b_0;
LL_Kp_b_1 = LL_Kp_b_0;
% 先确定单层对应的LL_K_m0, LL_Kp_m0, 
[model.eig_info_HK_m_EB_select_cell, model.eig_info_HKp_m_EB_select_cell] = sort_LLs_m(model.eig_info_HK_m_EB_select_cell, model.eig_info_HKp_m_EB_select_cell, ...
                                                                                       LL_K_m_0, LL_Kp_m_0, LL_index_cutoff, B_steps, ene_eps, weight);
                                                                                   
% 寻找LL_K_b0和LL_K_b1的指标
weight_K_b0 = weight;
weight_Kp_b0 = weight;
flag_slope_K_b0 = false;
flag_slope_Kp_b0 = false;
[model.eig_info_HK_b_EB_select_cell, model.eig_info_HKp_b_EB_select_cell, model.eigval_LL_K_b0_list, model.eigval_LL_Kp_b0_list] = ...
                    sort_LLs_b(model.eig_info_HK_b_EB_select_cell, model.eig_info_HKp_b_EB_select_cell, LL_K_b_0, LL_Kp_b_0, LL_index_cutoff, B_steps, ene_eps, weight_K_b0, weight_Kp_b0);                                                                                   

% 定出cnp_index1和cnp_index2
index1 = model.eig_info_HK_b_EB_select_cell{model.B_steps, 4} == 0;
eigval_temp1 = model.eig_info_HK_b_EB_select_cell{model.B_steps, 2}(index1);
index2 = model.eig_info_HK_b_EB_select_cell{B_steps, 4} == 1;
eigval_temp2 = model.eig_info_HK_b_EB_select_cell{B_steps, 2}(index2);

eigvals_LL_sort = sort([model.eigvals_LL_K_EB(B_steps, :), model.eigvals_LL_Kp_EB(B_steps, :)], 'ascend');
cnp_index1 = find(eigvals_LL_sort == eigval_temp1); % 也适用于ED，只要hopping_params相同即可
cnp_index2 = find(eigvals_LL_sort == eigval_temp2);                
model = model.trilayer_ABA_LLs_EB_cnp_shift(cnp_index1, cnp_index2);

[model, dos_HK_m_EB_mat, dos_HKp_m_EB_mat, dos_HK_b_EB_mat, dos_HKp_b_EB_mat] = model.trilayer_ABA_LLs_nB(0.001, 2, 201);

% 作图 LL EB 2D imges
dos_EB_mat_cell = cell(4, 1);
dos_EB_mat_cell{1} = dos_HK_m_EB_mat;
dos_EB_mat_cell{2} = dos_HKp_m_EB_mat;
dos_EB_mat_cell{3} = dos_HK_b_EB_mat;
dos_EB_mat_cell{4} = dos_HKp_b_EB_mat;
save_path = "";
fig1 = plot_some_imgs_in_same_xy(2, 2, model.B_fields_list, model.energy_list, dos_EB_mat_cell, save_path);
fig2 = model.representative_imgs_plot();

% 计算LL ED
B_field = 3.0;
Delta1_steps = 241;
Delta1_start = 0.0;
Delta1_end = 0.4 * 0.1; % 与实验数据匹配
Delta1_list = linspace(Delta1_start, Delta1_end, Delta1_steps); % 以eV为单位

model = model.trilayer_ABA_LLs_ED(B_field, Delta1_start, Delta1_end, Delta1_steps, LL_index_cutoff, ene_ub, ene_lb, ene_eps);
% 
line_info_cell = cell(4,2);
line_info_cell{1,1} = 'b--';
line_info_cell{1,2} = 0.5;
line_info_cell{2,1} = 'b';
line_info_cell{2,2} = 1.0;
fig0 = model.representative_lines_plot(line_info_cell, ene_lb, ene_ub, save_path);

model = model.trilayer_ABA_LLs_ED_cnp_shift(cnp_index1, cnp_index2);
[model, dos_HK_ED_mat, dos_HKp_ED_mat] = model.trilayer_ABA_LLs_nD(0.001, 2, 201);

% 作图LL ED images
dos_ED_mat_cell = cell(2, 1);
dos_ED_mat_cell{1} = dos_HK_ED_mat;
dos_ED_mat_cell{2} = dos_HKp_ED_mat;
save_path = "";
fig3 = plot_some_imgs_in_same_xy(1, 2, model.Delta1_list, model.energy_list, dos_ED_mat_cell, save_path);
fig4 = model.representative_imgs_plot();

%% 保存数据
filename = '_1GPa_model.mat';
save_path = ['D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\test_data\model_data\', datestr(datetime, 'yy-mm-dd-HH-MM-SS'), filename];
save(save_path);

%%
% % 找出在较高磁场(B=4.5T/10T/14T)下的LL_K_m_0 // LL_Kp_m_0在相互作用下导致的quantum hall ferromagnetism
% model.flag_EB = true;
% model.flag_ED = false;
% field1 = 4.5;
% field_index1 = model.find_closet_field_index(field1);
% field2 = 10.0;
% field_index2 = model.find_closet_field_index(field2);
% 
% LL_K_m0_eigvals_list = zeros(B_steps, 1);
% LL_Kp_m0_eigvals_list = zeros(B_steps, 1);
% LL_K_b0_eigvals_list = zeros(B_steps, 1);
% LL_Kp_b0_eigvals_list = zeros(B_steps, 1);
% LL_K_b1_eigvals_list = zeros(B_steps, 1);
% LL_Kp_b1_eigvals_list = zeros(B_steps, 1);
% for B_index = 1:B_steps
%     LL_K_m0_index = find(model.eig_info_HK_m_EB_select_cell{B_index, 4} == 0);
%     LL_K_m0_eigvals_list(B_index) = model.eig_info_HK_m_EB_select_cell{B_index, 2}(LL_K_m0_index);
%     
%     LL_Kp_m0_index = find(model.eig_info_HKp_m_EB_select_cell{B_index, 4} == 0);
%     LL_Kp_m0_eigvals_list(B_index) = model.eig_info_HKp_m_EB_select_cell{B_index, 2}(LL_Kp_m0_index);
%     
%     LL_K_b0_index = find(model.eig_info_HK_b_EB_select_cell{B_index, 4} == 0);
%     LL_K_b0_eigvals_list(B_index) = model.eig_info_HK_b_EB_select_cell{B_index, 2}(LL_K_b0_index);
%     
%     LL_Kp_b0_index = find(model.eig_info_HKp_b_EB_select_cell{B_index, 4} == 0);
%     LL_Kp_b0_eigvals_list(B_index) = model.eig_info_HKp_b_EB_select_cell{B_index, 2}(LL_Kp_b0_index);
%     
%     LL_K_b1_index = find(model.eig_info_HK_b_EB_select_cell{B_index, 4} == 1);
%     LL_K_b1_eigvals_list(B_index) = model.eig_info_HK_b_EB_select_cell{B_index, 2}(LL_K_b1_index);
%     
%     LL_Kp_b1_index = find(model.eig_info_HKp_b_EB_select_cell{B_index, 4} == 1);
%     LL_Kp_b1_eigvals_list(B_index) = model.eig_info_HKp_b_EB_select_cell{B_index, 2}(LL_Kp_b1_index);
% end
% fig8 = figure;
% hold on
% plot(B_fields_list, LL_K_m0_eigvals_list, 'r--')
% plot(B_fields_list, LL_Kp_m0_eigvals_list)
% plot(B_fields_list, LL_Kp_b0_eigvals_list)
% plot(B_fields_list, LL_K_b0_eigvals_list)
% plot(B_fields_list, LL_Kp_b1_eigvals_list)
% plot(B_fields_list, LL_K_b1_eigvals_list)

% LL_K_m_0 // LL_Kp_m_0是完全极化的，并且只和|0>有关


% field3 = 14.0;
% field_index3 = model.find_closet_field_index(field3);


% fig4 = figure;
% find(model.eig_info_HK_m_EB_select_cell{model.B_steps, 4} == 0)
