function current_fig = plot_LLs_weights(eigvecs_LL)
    % 画出朗道能级在每个分量下的权重
    current_fig = figure;
    hold on
    for ii = 1:size(eigvecs_LL, 2)
        plot(abs(eigvecs_LL(:, ii)).^2)
    end
end