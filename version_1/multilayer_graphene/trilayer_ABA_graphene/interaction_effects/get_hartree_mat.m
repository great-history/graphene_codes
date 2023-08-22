function hartree_mat = get_hartree_mat(Delta_mid, E_H, transform_mat_cell, LL_index_max)
    % 开始计算矩阵元<a|U_H|b>
    layer_potential_mat = diag([-1, -1, 1, 1, -1, -1]);
    hartree_mat = zeros(3);
    for n = 0:LL_index_max
        transform_mat_left = (transform_mat_cell{n + 1})';
        transform_mat_right = transform_mat_cell{n + 1};
        hartree_mat = hartree_mat + (transform_mat_left * layer_potential_mat * transform_mat_right);
    end
    hartree_mat = E_H / 2 * real(Delta_mid) * hartree_mat;
end