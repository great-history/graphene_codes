function [eig_vals_K_low_ene_gssp_cache, eig_vecs_K_low_ene_gssp_cache, dim_gssp_list, num_gssp] = get_gssp(eig_vecs_K_low_ene_cell, eig_vals_K_low_ene_cell, num_band_list, ...
                                                                                                            subspace_index_cell, subspace_range_cell, num_k, dims, overlap_eps)
    %% 确定彼此有能隙的子空间
    subspace_overlap_info_cell = cell(num_k, 1); % 存放subspace相互重叠的信息
    % k_1
    num_ssp_temp = subspace_index_cell{1}(end);
    subspace_overlap_info_cell{1} = cell(num_ssp_temp, 1);
    for ii = 1:num_ssp_temp
        subspace_overlap_info_cell{1}{ii} = ii;
    end

    % 从k_2到k_end
    % overlap_eps = 0.01; % 如果重叠小于overlap_eps的话就被认为是有重叠的，一般overlap_eps取为0.1 / 0.01，如果取1e-3 / 1e-4则会把所有子空间都Mix在一起

    for k_index = 2:num_k
        % calculate overlap matrix among bands
        num_band_temp = num_band_list(k_index);
        num_band_last = num_band_list(k_index - 1);
        eig_vecs_temp = eig_vecs_K_low_ene_cell{k_index};
        eig_vecs_last = eig_vecs_K_low_ene_cell{k_index - 1};

        band_overlap_matrix = calc_band_overlap_matrix(eig_vecs_temp, eig_vecs_last, num_band_temp, num_band_last);

        % calculate overlap matrix among subspaces
        num_ssp_temp = subspace_index_cell{k_index}(end);
        num_ssp_last = subspace_index_cell{k_index - 1}(end);
        subspace_range_temp = subspace_range_cell{k_index};
        subspace_range_last = subspace_range_cell{k_index - 1};
        ssp_overlap_matrix = get_subspace_overlap_matrix_by_band_overlap_matrix(band_overlap_matrix, ...
                                                                                subspace_range_temp, subspace_range_last, num_ssp_temp, num_ssp_last);

        % 得到子空间重叠的信息
        overlap_info_temp = cell(num_ssp_temp, 1);
        overlap_info_last = subspace_overlap_info_cell{k_index - 1};

        for ssp_index = 1:num_ssp_temp
            [overlap_list, new_order] = sort(ssp_overlap_matrix(ssp_index, :), 'descend');

            this_info = [];
            for ii = 1:num_ssp_last
                if overlap_list(ii) > overlap_eps
                    this_info = union(this_info, overlap_info_last{new_order(ii)});
                else
                    break
                end
            end

            overlap_info_temp{ssp_index} = this_info;
        end

        subspace_overlap_info_cell{k_index} = overlap_info_temp;
    end

    % 从subspace_overlap_info_cell{end}(即为overlap_info_temp)确定出彼此有空隙的子空间
    gssp_info_cell = {}; % gssp : gapped subspace
    ssp_start = 0;
    for ssp_index = 1:num_ssp_temp  % gssp_temp有可能是空的，即[]
        if ~isempty(overlap_info_temp{ssp_index})
            gssp_last = overlap_info_temp{ssp_index};
            ssp_start = ssp_index;
            break
        end
    end

    if ~(ssp_start == 0)
        count = 1;
    end

    for ssp_index = ssp_start:num_ssp_temp
        gssp_temp = overlap_info_temp{ssp_index};
        if isempty(gssp_temp) % gssp_temp有可能是空的，即[]
            continue
        end

        if isempty(intersect(gssp_last, gssp_temp))
            gssp_info_cell{count,1} = gssp_last;
            count = count + 1;
            gssp_last = gssp_temp;
        else
            gssp_last = union(gssp_last, gssp_temp);
        end

        if ssp_index == num_ssp_temp
            gssp_info_cell{count,1} = gssp_last;
        end
    end
    num_gssp = count;

    % 分类：将所有能带按照gssp进行分类,第一个gssp和最后一个gssp可以扔掉，因为有外界的能带Mix进来
    eig_vals_K_low_ene_gssp_cell = cell(num_k, num_gssp - 2);
    eig_vecs_K_low_ene_gssp_cell = cell(num_k, num_gssp - 2);
    for k_index = 1:num_k
        eig_vals_K_low_ene_temp = eig_vals_K_low_ene_cell{k_index};
        eig_vecs_K_low_ene_temp = eig_vecs_K_low_ene_cell{k_index};
        overlap_info_temp = subspace_overlap_info_cell{k_index};
        ssp_range_temp = subspace_range_cell{k_index};
        num_ssp_temp = length(overlap_info_temp);

        for gssp_index = 2:(num_gssp - 1)
            gssp_temp = gssp_info_cell{gssp_index};

            % 初始化
            eig_vals_K_low_ene_gssp_temp = [];
            eig_vecs_K_low_ene_gssp_temp = [];

            for ssp_index = 1:num_ssp_temp
                ssp_temp = overlap_info_temp{ssp_index};
                if ~(isempty(intersect(gssp_temp, ssp_temp)))
                    for band_index = ssp_range_temp(ssp_index, 1):ssp_range_temp(ssp_index, 2)
                        eig_vals_K_low_ene_gssp_temp = [eig_vals_K_low_ene_gssp_temp, eig_vals_K_low_ene_temp(band_index)];
                        eig_vecs_K_low_ene_gssp_temp = [eig_vecs_K_low_ene_gssp_temp, eig_vecs_K_low_ene_temp(:, band_index)];
                    end
                end
            end

            eig_vals_K_low_ene_gssp_cell{k_index, gssp_index - 1} = eig_vals_K_low_ene_gssp_temp;
            eig_vecs_K_low_ene_gssp_cell{k_index, gssp_index - 1} = eig_vecs_K_low_ene_gssp_temp;
        end
    end


    num_gssp = num_gssp - 2;
    % convert cell to mat
    eig_vals_K_low_ene_gssp_cache = cell(num_gssp, 1);
    eig_vecs_K_low_ene_gssp_cache = cell(num_gssp, 1);

    dim_gssp_list = zeros(num_gssp, 1);
    for gssp_index = 1:num_gssp
       dim_temp = length(eig_vals_K_low_ene_gssp_cell{1, gssp_index});
       dim_gssp_list(gssp_index) = dim_temp;

       eig_vals_K_low_ene_gssp_list = cell2mat(eig_vals_K_low_ene_gssp_cell(:, gssp_index)); % dims = (num_k, dim_temp)

       eig_vecs_K_low_ene_gssp_list = zeros(num_k, dims, dim_temp); 
       for k_index = 1:num_k
           eig_vecs_K_low_ene_gssp_list(k_index, :, :) = eig_vecs_K_low_ene_gssp_cell{k_index, gssp_index};
       end

       eig_vals_K_low_ene_gssp_cache{gssp_index} = eig_vals_K_low_ene_gssp_list;
       eig_vecs_K_low_ene_gssp_cache{gssp_index} = eig_vecs_K_low_ene_gssp_list;
    end

end