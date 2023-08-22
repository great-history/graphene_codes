function filling_list = filling_factor_func(type, ene_list, dos_list, ene_num_points, w_bandwidth)
    filling_list = zeros(1, ene_num_points);
    
    % different types of DOS
    % Linear DOS
    if type == "linear" % 有解析表达式
        for ii = 1:ene_num_points
            ene_current = ene_list(ii);
            %             if ene_current <= - w_bandwidth
            %                 filling_list(ii) = - 1;
            %             elseif (- w_bandwidth < ene_current) && (ene_current <= 0)
            %                 filling_list(ii) = - (ene_current / w_bandwidth)^2;
            %             elseif (ene_current < w_bandwidth) && (ene_current > 0)
            %                 filling_list(ii) = (ene_current / w_bandwidth)^2;
            %             else
            %                 filling_list(ii) = 1;
            %             end
            filling_list(ii) = filling_factor_func_linear(ene_current, w_bandwidth);
        end
    else
        % 使用梯形法近似积分
        filling_list = cumtrapz(ene_list, dos_list);
        filling_list = filling_list - 1.0; % ensure the filling factore at the lowest energy is - 1
    end

end