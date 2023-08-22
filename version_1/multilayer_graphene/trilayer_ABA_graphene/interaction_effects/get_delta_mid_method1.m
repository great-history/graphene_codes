function Delta_mid = get_delta_mid_method1( density_matrix_up_temp_LL, density_matrix_down_temp_LL, LL_index_max)
    % 第一种方式 % 与第一种方法得到的结果应该是一致的，推荐使用第二种，算起来会快一些
    Delta_mid = 0.0;
    for n = 0:LL_index_max
        index_start = 6 * n;
        Delta_mid = Delta_mid + density_matrix_down_temp_LL(index_start + 3, index_start + 3);
        Delta_mid = Delta_mid + density_matrix_down_temp_LL(index_start + 4, index_start + 4);
        Delta_mid = Delta_mid + density_matrix_up_temp_LL(index_start + 3, index_start + 3);
        Delta_mid = Delta_mid + density_matrix_up_temp_LL(index_start + 4, index_start + 4);
    end
end