addpath('.\classes\')
% 创建模型
load('test_data\23-04-07-14-01-46_0GPa_without_Delta1.mat')
hopping_params_before = output_value_array(1,:);
model_before = trilayer_ABA_class(hopping_params_before);
load('test_data\23-04-07-20-23-52_1GPa_without_Delta1.mat')
hopping_params_after = output_value_array(1,:);
model_after = trilayer_ABA_class(hopping_params_after);

% 创建hdf5文件并保存数据集
h5create('./code_test/hopping_parameters.h5', '/hopping_parameters_before', [1, 8]);
h5write('./code_test/hopping_parameters.h5', '/hopping_parameters_before', hopping_params_before);
h5create('./code_test/hopping_parameters.h5', '/hopping_parameters_after', [1, 8]);
h5write('./code_test/hopping_parameters.h5', '/hopping_parameters_after', hopping_params_after);

%% 计算LL EB
Delta1 = 0.0; 
B_start = 0.1;
B_end = 14;
B_steps = 501;
B_fields_list = linspace(B_start, B_end, B_steps);
LL_index_cutoff = 30;
ene_ub = 0.1;
ene_lb = -0.1;
ene_eps = 0.002;
model_before = model_before.trilayer_ABA_LLs_EB(Delta1, B_start, B_end, B_steps, LL_index_cutoff, ene_ub, ene_lb, ene_eps);
model_after = model_after.trilayer_ABA_LLs_EB(Delta1, B_start, B_end, B_steps, LL_index_cutoff, ene_ub, ene_lb, ene_eps);

% 作图LL EB
line_info_cell = cell(4,2);
line_info_cell{1,1} = 'r';
line_info_cell{1,2} = 1.0;
line_info_cell{2,1} = 'r--';
line_info_cell{2,2} = 1.0;
line_info_cell{3,1} = 'b--';
line_info_cell{3,2} = 1.0;
line_info_cell{4,1} = 'b';
line_info_cell{4,2} = 1.0;
save_path = "";

% //////////////////////////////////////////////////// 0GPa
eigvals_LL_cell_before = cell(4, 4); % 分为 m_K / m_Kp / b_K / b_Kp 四种
eigvals_LL_cell_before{1,1} = model_before.eigvals_LL_K_EB(:, 1:model_before.dims_m); % m_K
eigvals_LL_cell_before{1,2} = model_before.dims_m; % 维数
eigvals_LL_cell_before{1,3} = line_info_cell{1, 1}; % line color 'b--'
eigvals_LL_cell_before{1,4} = line_info_cell{1, 2}; % line width

eigvals_LL_cell_before{2,1} = model_before.eigvals_LL_Kp_EB(:, 1:model_before.dims_m); % m_Kp
eigvals_LL_cell_before{2,2} = model_before.dims_m;
eigvals_LL_cell_before{2,3} = line_info_cell{2, 1};
eigvals_LL_cell_before{2,4} = line_info_cell{2, 2};

eigvals_LL_cell_before{3,1} = model_before.eigvals_LL_K_EB(:, (1 + model_before.dims_m):end);
eigvals_LL_cell_before{3,2} = model_before.dims_b; % 维数
eigvals_LL_cell_before{3,3} = line_info_cell{3, 1}; % line color 'b--'
eigvals_LL_cell_before{3,4} = line_info_cell{3, 2}; % line width

eigvals_LL_cell_before{4,1} = model_before.eigvals_LL_Kp_EB(:, (1 + model_before.dims_m):end);
eigvals_LL_cell_before{4,2} = model_before.dims_b;
eigvals_LL_cell_before{4,3} = line_info_cell{4, 1};
eigvals_LL_cell_before{4,4} = line_info_cell{4, 2};

% //////////////////////////////////////////////////// 1GPa
eigvals_LL_cell_after = cell(4, 4); % 分为 m_K / m_Kp / b_K / b_Kp 四种
eigvals_LL_cell_after{1,1} = model_after.eigvals_LL_K_EB(:, 1:model_after.dims_m); % m_K
eigvals_LL_cell_after{1,2} = model_after.dims_m; % 维数
eigvals_LL_cell_after{1,3} = line_info_cell{1, 1}; % line color 'b--'
eigvals_LL_cell_after{1,4} = line_info_cell{1, 2}; % line width

eigvals_LL_cell_after{2,1} = model_after.eigvals_LL_Kp_EB(:, 1:model_after.dims_m); % m_Kp
eigvals_LL_cell_after{2,2} = model_after.dims_m;
eigvals_LL_cell_after{2,3} = line_info_cell{2, 1};
eigvals_LL_cell_after{2,4} = line_info_cell{2, 2};

eigvals_LL_cell_after{3,1} = model_after.eigvals_LL_K_EB(:, (1 + model_after.dims_m):end);
eigvals_LL_cell_after{3,2} = model_after.dims_b; % 维数
eigvals_LL_cell_after{3,3} = line_info_cell{3, 1}; % line color 'b--'
eigvals_LL_cell_after{3,4} = line_info_cell{3, 2}; % line width

eigvals_LL_cell_after{4,1} = model_after.eigvals_LL_Kp_EB(:, (1 + model_after.dims_m):end);
eigvals_LL_cell_after{4,2} = model_after.dims_b;
eigvals_LL_cell_after{4,3} = line_info_cell{4, 1};
eigvals_LL_cell_after{4,4} = line_info_cell{4, 2};

%% 作图
current_fig = figure;
subplot(1,2,1)
% axis([B_start B_end 0.075 0.125])
axis([B_start B_end -0.04 0.04])
hold on

% plot LL
for index = 1:size(eigvals_LL_cell_before, 1)
    eigvals_LL_current = eigvals_LL_cell_before{index, 1};
    dims = eigvals_LL_cell_before{index, 2};
    line_color = eigvals_LL_cell_before{index, 3};
    line_width = eigvals_LL_cell_before{index, 4};

    for i = 1:dims
        plot(B_fields_list, eigvals_LL_current(:,i), line_color, 'LineWidth', line_width)
    end
end
set(gca, 'XDir', 'normal')
set(gca,'xtick',[0 2 4 6 8 10 12 14]);
set(gca,'xticklabel',[]);
set(gca, 'YDir', 'normal')
set(gca,'ytick',[-0.04, -0.02, 0.00, 0.02, 0.04]);
set(gca,'yticklabel',[]);

subplot(1,2,2)
% axis([B_start B_end 0.075 0.125])
axis([B_start B_end -0.04 0.04])
hold on

% plot LL
for index = 1:size(eigvals_LL_cell_after, 1)
    eigvals_LL_current = eigvals_LL_cell_after{index, 1};
    dims = eigvals_LL_cell_after{index, 2};
    line_color = eigvals_LL_cell_after{index, 3};
    line_width = eigvals_LL_cell_after{index, 4};

    for i = 1:dims
        plot(B_fields_list, eigvals_LL_current(:,i), line_color, 'LineWidth', line_width)
    end
end
set(gca, 'XDir', 'normal')
set(gca,'xtick',[0 2 4 6 8 10 12 14]);
set(gca,'xticklabel',[]);
set(gca, 'YDir', 'normal')
set(gca,'ytick',[-0.04, -0.02, 0.00, 0.02, 0.04]);
set(gca,'yticklabel',[]);

%% MATLAB模式化作图
% figure
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　画布(figure)和坐标轴(axes)的尺寸设置 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(gcf, 'unit', 'inch', 'position', [0,0,2.75, 3]) % figure
% set(gca, 'Position', [0.18, 0.15, 0.75, 0.75]) % axes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　坐标轴的方向                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(gca, 'XDir', 'normal')
% set(gca, 'YDir', 'normal')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　坐标轴标签的命                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xlabel('B(T)', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial')
% ylabel('E(meV)', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　坐标轴的范围                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis([B_start B_end -0.04 0.04])
% hold on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　坐标轴上的标度                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(gca,'xtick', [2 4 6 8 10 12 14]);
% set(gca,'xticklabel', [2 4 6 8 10 12 14]);
% set(gca,'ytick', [-0.04, -0.02, 0.00, 0.02, 0.04]);
% set(gca,'yticklabel', [-0.04, -0.02, 0.00, 0.02, 0.04]);
% set(gca, 'TickLength', [0.01, 0.03])  % 第一个是Minor Tick，第二个是Major Tick，类似于Origin
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　绘制曲线                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % mK
% index = 1;
% eigvals_LL_current = eigvals_LL_cell_after{index, 1};
% line_color = 'b';
% line_width = 2.0;
% h1 = plot(B_fields_list, eigvals_LL_current(:,index), line_color, 'LineWidth', line_width);
% % mKp
% index = 2;
% eigvals_LL_current = eigvals_LL_cell_after{index, 1};
% line_color = 'b';
% line_width = 2.0;
% h2 = plot(B_fields_list, eigvals_LL_current(:,index), line_color, 'LineWidth', line_width);
% % bK
% index = 3;
% eigvals_LL_current = eigvals_LL_cell_after{index, 1};
% line_color = 'r';
% line_width = 2.0;
% h3 = plot(B_fields_list, eigvals_LL_current(:,index), line_color, 'LineWidth', line_width);
% % bKp
% index = 1;
% eigvals_LL_current = eigvals_LL_cell_after{index, 1};
% dims = eigvals_LL_cell_after{index, 2};
% line_color = 'r';
% line_width = 2.0;
% h4 = plot(B_fields_list, eigvals_LL_current(:,index), line_color, 'LineWidth', line_width);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　设置图例(Legend)                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend([h1(1), h2(1), h3(1), h4(1)], 'm-K', 'm-Kp', 'b-K', 'b-Kp')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　将坐标轴的框置于顶层                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(gca, 'Layer', 'top')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　框样式的设置                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(gca, 'Layer', 'top')
% set(gca, 'Box', 1) % 框的四边
% set(gcf, 'color', 'none') % 背景透明化 
% set(gca, 'color', 'none')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　(手动)保存成eps格式                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % print(['figure_name', '.eps'])
% % legend
% % h1 = plot(B_fields_list, eigvals_LL_current, 'b', 'LineWidth', 1);
% % legend([h1(1), h1(2)],'asd','sd')