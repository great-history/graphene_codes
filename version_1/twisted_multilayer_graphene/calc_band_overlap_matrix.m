%% 计算overlap_matrix among bands
function overlap_matrix = calc_band_overlap_matrix(eig_vecs_temp, eig_vecs_last, num_band_temp, num_band_last)
    overlap_matrix = zeros(num_band_temp, num_band_last);
    for i = 1:num_band_temp
        vec_temp = eig_vecs_temp(:, i);
        for j = 1:num_band_last
            vec_last = eig_vecs_last(:, j);
            overlap_matrix(i,j) = (abs(dot(vec_temp, vec_last)))^2;
        end
    end
end