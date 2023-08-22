function current_fig = plot_LLs_k_dist(phi_k_mesh, kx_list, ky_list, num_img, save_path)
    kx_min = kx_list(1);
    kx_max = kx_list(end);
    ky_min = ky_list(1);
    ky_max = ky_list(end);
    num_points = length(kx_list);
    
    current_fig = figure;
    set(gcf,'position',[400 400 400 * num_img 400])
    % set(gcf,'position',[250 300 1000 900])
    for lambda = 1:num_img
        subplot(1, num_img, lambda)
        im1 = imagesc(kx_list, ky_list, reshape(phi_k_mesh(lambda, :, :), num_points, num_points));
        axis([kx_min kx_max ky_min ky_max])
        box on
        shading interp
        hold on

        %kx_valley = k_valley_list(lambda, 1);
        %ky_valley = k_valley_list(lambda, 2);
        %scatter(kx_valley, ky_valley, 'filled', 'd')

        set(gca, 'YDir', 'normal')
        xlabel('$\frac{\hbar v k_x}{\gamma_1}$','interpreter','latex', 'FontSize', 12);
        ylabel('$\frac{\hbar v k_y}{\gamma_1}$','interpreter','latex', 'FontSize', 12);
        colormap jet
        colorbar
        title('$|\phi(\vec{k})|^2$', 'interpreter', 'latex', 'FontSize', 14);
    end
    
    if ~(save_path == "")
        saveas(gcf, save_path); %保存当前窗口的图像
    end
end