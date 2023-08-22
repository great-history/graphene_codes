function dos_list = get_dos_list(eig_enes, ene_list, sigma)
    % sigma以eV为单位
    
    num_ene_point = length(ene_list);
    dos_list = zeros(size(ene_list));
    if length(size(eig_enes_K)) == 3  % 二维k空间
        num_band = size(eig_enes, 3); % 能带数目
        dim_x = size(eig_enes, 2);
        dim_y = size(eig_enes, 1);
        
        for i = 1:num_ene_point
            val = 0.0;
            ene_val = ene_list(i);
            for ii = 1:dim_y
                for jj = 1:dim_x
                    for kk = 1:num_band
                        eig_ene = eig_enes(ii, jj, kk);
                        val = val + 1 / (sigma * sqrt(2 * pi)) * exp(-(ene_val - eig_ene)^2 / (2 * sigma^2));
                    end
                end
            end
            dos_list(i) = val;
        end
        
    else % length(size(eig_enes_K)) == 2
        num_band = size(eig_enes, 2); % 能带数目
        dim = size(eig_enes, 1);
        
        for i = 1:num_ene_point
            val = 0.0;
            ene_val = ene_list(i);
            for ii = 1:dim
                for kk = 1:num_band
                    eig_ene = eig_enes(ii, kk);
                    val = val + 1 / (sigma * sqrt(2 * pi)) * exp(-(ene_val - eig_ene)^2 / (2 * sigma^2));
                end
            end
            dos_list(i) = val;
        end
    end
    
end