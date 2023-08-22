function R_mat = get_R_mat(transform_mat_cell, LL_index_max, dim_subspace)
    R_mat = zeros(dim_subspace);
    layer_mat = diag([-1, -1, 1, 1, -1, -1]);
    for n = 0:LL_index_max
        Wa = (transform_mat_cell{n + 1})';
        Wb = transform_mat_cell{n + 1};
        R_mat = R_mat + Wa * layer_mat * Wb;
    end
end