function [carrier_density_vecs, DOS_nB_matrix] = helper_get_LLs_nB(DOS_EB_matrix, carrier_density_matrix, carrier_density_ub, carrier_density_lb, n_points, B_steps)
    carrier_density_vecs = linspace(carrier_density_lb, carrier_density_ub, n_points);

    % 开始作插值

    DOS_nB_matrix = zeros(n_points, B_steps);

    % tic
    for B_index = 1:B_steps
        x = carrier_density_matrix(:, B_index);
        v = DOS_EB_matrix(:,B_index);
        xq = carrier_density_vecs;

        vq2 = interp1(x,v,xq, 'linear'); % 不要使用spline去插值，因为边界会由问题
        DOS_nB_matrix(:, B_index) = vq2;
    end
    % toc 
end