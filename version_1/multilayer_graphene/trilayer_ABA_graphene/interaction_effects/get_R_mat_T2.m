function R_mat = get_R_mat_T2(transform_mat_cell, LL_index_max)
    R_mat = zeros(3);
    layer_mat = diag([-1, -1, 1, 1, -1, -1]);
    for n = 0:LL_index_max
        Wa = (transform_mat_cell{n + 1})';
        Wb = transform_mat_cell{n + 1};
        R_mat = R_mat + Wa * layer_mat * Wb;
    end
end