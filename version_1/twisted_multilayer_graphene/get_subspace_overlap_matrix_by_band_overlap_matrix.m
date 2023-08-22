%% 计算overlap_matrix among subspaces
function ssp_overlap_matrix = get_subspace_overlap_matrix_by_band_overlap_matrix(band_overlap_matrix, subspace_range_temp, subspace_range_last, num_ssp_temp, num_ssp_last)
    ssp_overlap_matrix = zeros(num_ssp_temp, num_ssp_last);
    for i = 1:num_ssp_temp
        temp_start = subspace_range_temp(i,1);
        temp_end = subspace_range_temp(i,2);

        for j = 1:num_ssp_last
            last_start = subspace_range_last(j,1);
            last_end = subspace_range_last(j,2);
            
            inner_dot_sum = sum(sum(band_overlap_matrix(temp_start:temp_end, last_start:last_end)));
            ssp_overlap_matrix(i, j) = inner_dot_sum;
        end
    end
end