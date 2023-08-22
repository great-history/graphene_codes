function [carrier_density_vecs, DOS_nD_matrix] = helper_get_LLs_nD(DOS_ED_matrix, carrier_density_matrix, carrier_density_ub, carrier_density_lb, n_points, Delta1_steps)
    carrier_density_vecs = linspace(carrier_density_lb, carrier_density_ub, n_points);

    % 开始作插值

    DOS_nD_matrix = zeros(n_points, Delta1_steps);

    % tic
    for Delta1_index = 1:Delta1_steps
        x = carrier_density_matrix(:, Delta1_index);
        v = DOS_ED_matrix(:, Delta1_index);
        xq = carrier_density_vecs;

        vq2 = interp1(x,v,xq, 'linear'); % 不要使用spline去插值，因为边界会由问题
        DOS_nD_matrix(:, Delta1_index) = vq2;
    end
    % toc 
end