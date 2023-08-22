function current_fig = plot_some_imgs_in_same_xy(row, col, x_list, y_list, data_cell, save_path)
    % 所有的图共用相同的x,y坐标
    current_fig = figure; 
    
    for ii = 0:(row - 1)
        for jj = 1:col
            subplot(row, col, ii * row + jj)
            imagesc(x_list, y_list, (data_cell{ii * row + jj})');
            set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
            colorbar
        end
    end
    
    if ~(save_path == "")
        saveas(gcf, save_path); %保存当前窗口的图像
    end
end