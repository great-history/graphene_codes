function current_fig = plot_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eigvals_LL_cell, save_path)
    % 画出朗道扇形图(横坐标为磁场(T),纵坐标为能量(eV))
    % eigvals_LL_cell包含了想要画的数据，之所以用cell是因为有可能数据大小不一致
    current_fig = figure;
    % axis([B_start B_end 0.075 0.125])
    axis([B_start B_end E_bottom E_top])
    hold on

    % plot LL
    for index = 1:size(eigvals_LL_cell, 1)
        eigvals_LL_current = eigvals_LL_cell{index, 1};
        dims = eigvals_LL_cell{index, 2};
        line_color = eigvals_LL_cell{index, 3};
        line_width = eigvals_LL_cell{index, 4};
        
        for i = 1:dims
            plot(B_fields_list, eigvals_LL_current(:,i), line_color, 'LineWidth', line_width)
        end
        
    end
    
    if ~(save_path == "")
        saveas(gcf, save_path); %保存当前窗口的图像
    end
    
end