function dos_list = dos_func(type, ene_list, ene_num_points, w_bandwidth)
    %% get dos list
    % common parameters
    % type = "linear" / "high_order_VHS" / "asymmetry" / "continuum_model"; % 选取DOS函数的类型
    dos_list = zeros(1, ene_num_points);
    
    % different types of DOS
    % Linear DOS
    if type == "linear"
        for ii = 1:ene_num_points
            ene_current = ene_list(ii);
            if abs(ene_current) < w_bandwidth
                dos_list(ii) = 2 * abs(ene_current) / (w_bandwidth)^2;
            else
                dos_list(ii) = 0;
            end
        end
    end

end