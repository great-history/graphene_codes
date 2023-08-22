function current_fig = plot_select_LLs(B_start, B_end, E_bottom, E_top, B_fields_list, eig_info_select_cell, color, save_path)
    % 画出朗道扇形图(横坐标为磁场(T),纵坐标为能量(eV))
    % eigvals_LL_cell包含了想要画的数据，之所以用cell是因为有可能数据大小不一致
    current_fig = figure;
    % axis([B_start B_end 0.075 0.125])
    axis([B_start B_end E_bottom E_top])
    hold on
    sz = 1;
    % plot selected LL
    for index = 1:length(B_fields_list)
        B_field = B_fields_list(index);
        eigvals_LL_current = eig_info_select_cell{index, 2};
        dims = eig_info_select_cell{index, 3};
        scatter(B_field * ones(dims, 1), eigvals_LL_current, sz, color)
    end

    if ~(save_path == "")
        saveas(gcf, save_path); %保存当前窗口的图像
    end
end