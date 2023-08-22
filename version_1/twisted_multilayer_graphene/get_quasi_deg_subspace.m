function [subspace_index_cell, subspace_range_cell] = get_quasi_deg_subspace(eig_vals_K_low_ene_cell, num_band_list, num_k, ene_eps)
    %% 确定准简并子空间
    subspace_index_cell = cell(num_k, 1);

    for k_index = 1:num_k
        num_band = num_band_list(k_index);
        subspace_index_list = zeros(num_band, 1);

        eig_vals_K_low_ene_diff = diff(eig_vals_K_low_ene_cell{k_index});

        count_subspace = 1;
        subspace_index_list(1) = count_subspace;

        for band_index = 2:num_band
            if eig_vals_K_low_ene_diff(band_index - 1) < ene_eps
                subspace_index_list(band_index) = count_subspace;
            else
                count_subspace = count_subspace + 1;
                subspace_index_list(band_index) = count_subspace;
            end
        end

        subspace_index_cell{k_index} = subspace_index_list;
    end

    subspace_range_cell = cell(num_k, 1);

    for k_index = 1:num_k
        num_band = num_band_list(k_index);
        num_subspace = subspace_index_cell{k_index}(end);
        sp_diff_list = diff(subspace_index_cell{k_index});
        subspace_range_list = zeros(num_subspace, 2);

        sp_index = 1;
        subspace_range_list(sp_index, 1) = 1;
        for band_index = 2:num_band
            if ~(sp_diff_list(band_index - 1) == 0)
                subspace_range_list(sp_index, 2) = band_index - 1;
                sp_index = sp_index + 1;
                subspace_range_list(sp_index, 1) = band_index;
            end
        end

        subspace_range_list(sp_index, 2) = num_band;
        subspace_range_cell{k_index} = subspace_range_list;
    end
end