function wave_packet_mesh = get_localized_wave_packet(xx_mesh, yy_mesh, dims_x, dims_y, LL_index)
    % xx_mesh, yy_mesh都是以mag_length为单位
    wave_packet_mesh = zeros(dims_x, dims_y);
    for ii = 1:dims_x
        for jj = 1:dims_y
            x_temp = xx_mesh(ii, jj);
            y_temp = yy_mesh(ii, jj);
            
            factor = 1 / sqrt(factorial(LL_index)) * ((x_temp - 1j * y_temp) / sqrt(2))^LL_index;
            mag_factor = exp( - (x_temp^2 + y_temp^2) / 4 );
            phase_factor = exp( 1j * x_temp * y_temp / 2 );
            wave_packet_mesh(ii, jj) = factor * mag_factor * phase_factor;
        end
    end

end

% 
% function probability_density_mesh = get_probability_density_mesh(xx_mesh, yy_mesh, dims_x, dims_y, mag_length, LL_index)
% 
% end

% dims_x = 201;
% dims_y = 201;
% x_list = linspace(-2, 2, dims_x);
% y_list = linspace(-2, 2, dims_y);
% [xx_mesh, yy_mesh] = meshgrid(x_list, y_list);
% wave_packet_mesh = get_localized_wave_packet(xx_mesh, yy_mesh, dims_x, dims_y, 0)