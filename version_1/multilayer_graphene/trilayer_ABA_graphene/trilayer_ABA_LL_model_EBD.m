addpath('.\classes\')
% 创建模型
hopping_params = output_value_array(1,:);
model = trilayer_ABA_class(hopping_params);

%% 计算LL EB
Delta1 = 0.0; 
B_start = 0.5;
B_end = 10;
B_steps = 501;
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
                                                                                   
%% 计算LL ED
B_field = 10.0;
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
fig1 = model.representative_lines_plot(line_info_cell, ene_lb, ene_ub, save_path);