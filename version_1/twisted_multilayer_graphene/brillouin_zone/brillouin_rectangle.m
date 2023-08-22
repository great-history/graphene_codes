function [akx_mesh, aky_mesh, k_K_list] = brillouin_rectangle(ld_corner, width, np_width, height, np_height, ak_norm)
    %% 撒点
    bz_x_array = linspace(0, width, np_width);
    bz_y_array = linspace(0, height, np_height);
    bz_x_array = bz_x_array + ld_corner(1);
    bz_y_array = bz_y_array + ld_corner(2);
    
    [akx_mesh, aky_mesh] = meshgrid(bz_x_array, bz_y_array);
    
    %% 旋转120°和240°得到其它的内部点
    Rotation_C3 = zeros(2);
    Rotation_C3(1,1) = -1/2;
    Rotation_C3(1,2) = -sqrt(3)/2;
    Rotation_C3(2,1) = sqrt(3)/2;
    Rotation_C3(2,2) = -1/2;
    
    %% 找出K和K'点
    k_K_list = [];
    k_K = [-sqrt(3)/2; 1/2];
    k_Kp = [-sqrt(3)/2; -1/2];
    k_K_list = [k_K_list, k_K, k_Kp, Rotation_C3 * k_K, Rotation_C3 * k_Kp, Rotation_C3 * Rotation_C3 * k_K, Rotation_C3 * Rotation_C3 * k_Kp];

    %% 但是对于转角石墨烯，最后还要再平移一个矢量
    k_K_list = k_K_list + [sqrt(3)/2;1/2];
    
    %% 最后还要乘以ak_norm
    akx_mesh = akx_mesh * ak_norm;
    aky_mesh = aky_mesh * ak_norm;
    k_K_list = k_K_list * ak_norm;
end