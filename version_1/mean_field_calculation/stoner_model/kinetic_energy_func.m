function kinetic_list = kinetic_energy_func(type, ene_list, dos_list, ene_num_points, w_bandwidth)
    % get the kinetic energy at zero temperature
    kinetic_list = zeros(1, ene_num_points);
    
    % different types of DOS
    % Linear DOS
    if type == "linear" % 有解析表达式
        for ii = 1:ene_num_points
            ene_current = ene_list(ii);
            if ene_current <= - w_bandwidth
                kinetic_list(ii) = - 2 / 3 * w_bandwidth;
            elseif ene_current >= w_bandwidth
                kinetic_list(ii) = 2 / 3 * w_bandwidth;
            else
                kinetic_list(ii) = 2 / 3 * (ene_current)^3 / (w_bandwidth)^2;
            end
        end
    else
        % 使用梯形法近似积分
        kinetic_list = cumtrapz(ene_list, ene_list .* dos_list);
    end
end