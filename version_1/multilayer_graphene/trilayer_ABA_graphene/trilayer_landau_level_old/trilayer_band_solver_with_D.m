function [eig_enes, k1s, k2s] = trilayer_band_solver_with_D(N_pb, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, m1, m2, k1_start, k2_start)
    % N_pb:周期性边界条件 gamma0:最近邻跃迁强度 m1和m2定义了多大的布里渊区, k1_start和k2_start定义了起点(不一定是原点)
    % 周期性边界条件
    N1 = N_pb * m1;
    N2 = N_pb * m2; % 假定不同方向的周期性条件是一样的，否则下面作图时可能会出错

    % 生成k-points
    k1s = zeros(N1,1);
    k2s = zeros(N2,1);
    for i = 1:N1
        k1s(i) = (i-1) / N_pb + k1_start;
    end

    for i = 1:N2
        k2s(i) = (i-1) / N_pb + k2_start;
    end

    % 计算相因子
    [f_k_mat, g_k_mat] = get_phase_factors(N1, N2, k1s, k2s);
    
    eig_enes = zeros(N1, N2, 6);
    for k2_ind = 1:N2
        for k1_ind = 1:N1
            f_k = f_k_mat(k1_ind, k2_ind);
            
            H_t_k = construct_trilayer_Ham_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, f_k);

            % get the eigenvals and eigenvecs
            [eigvecs, eigvals] = eig(H_t_k);
            % re-ordering(to do)

            % push eigvals into eig_enes
            eig_enes(k1_ind, k2_ind,:) = diag(eigvals);
        end
    end
end


function [f_k_mat, g_k_mat] = get_phase_factors(N1, N2, k1s, k2s)
    f_k_mat = zeros(N1,N2); % 第一个维度是沿着a1方向，第二个维度是沿着a2方向
    g_k_mat = zeros(N1,N2);

    Im = 1i;
    for k2_ind = 1:N2
        for k1_ind = 1:N1
            % \vec{\delta_1} = (-1/3, -2/3) // \vec{\delta_2} = (2/3, 1/3) // \vec{\delta_3} = (-1/3, 1/3)
            k1 = k1s(k1_ind);
            k2 = k2s(k2_ind);

            ang1 = 2 * pi * (- k1 / 3 - 2 * k2 / 3);
            ang2 = 2 * pi * (2 * k1 / 3 + k2 / 3);
            ang3 = 2 * pi * (- k1 / 3 + k2 / 3);
            f_k = exp(- Im * ang1) + exp(- Im * ang2) + exp(- Im * ang3);

            f_k_mat(k1_ind, k2_ind) = f_k;
            g_k_mat(k1_ind, k2_ind) = real(f_k)^2 + imag(f_k)^2 - 3;
        end
    end
end


% % example
% gamma0 = 3100;
% gamma1 = 355;
% gamma2 = -20.7904;
% gamma3 = 346.7153;  
% gamma4 = 91.2578;
% gamma5 = 35.7223;
% 
% D_field = 0;
% delta = 1.5+(gamma5-gamma2)/2;
% Delta1 = D_field * 0.1 * 1000;
% Delta2 = 1.8;
% 
% N_pb = 200;
% N1 = N_pb;
% N2 = N_pb;
% 
% [tri_eig_enes, k1s, k2s] = trilayer_band_solver_with_D(N_pb, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, 1, 1, 0, 0);
% 
% % 转到直角坐标系
% [kxs, kys] = convert_k_space(N1, N2, k1s, k2s);
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % example of trilayer on an extended Brillouin Zone
% figure  % 这个图是FBZ为honeycomb的情形，也是graphene最最常用的
% [eig_enes_super, k1s_super, k2s_super] = trilayer_band_solver_with_D(N_pb, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, 2, 2, -1, -1);
% 
% % 做一个扩大的布里渊区:转到直角坐标系
% [b_zone_x, b_zone_y] = convert_k_space(2 * N1, 2 * N2, k1s_super, k2s_super);
% 
% % eig_enes_super(1:N1, 1:N2, :) = eig_enes;
% % eig_enes_super((N1+1):2*N1, 1:N2, :) = eig_enes;
% % eig_enes_super(1:N1, (N2+1):2*N2, :) = eig_enes;
% % eig_enes_super((N1+1):2*N1, (N2+1):2*N2, :) = eig_enes;
% figure
% mesh(b_zone_x,b_zone_y,eig_enes_super(:,:,1))
% hold on
% mesh(b_zone_x,b_zone_y,eig_enes_super(:,:,2))
% colorbar
% 
% figure
% h1 = surf(b_zone_x, b_zone_y, eig_enes_super(:,:,1), 'FaceAlpha',0.8);
% hold on
% h2 = surf(b_zone_x, b_zone_y, eig_enes_super(:,:,2), 'FaceAlpha',0.8);
% % set(h2, 'edgecolor','k');
% % shading interp
% h1.EdgeColor = 'none';
% h1.LineStyle = '-';
% h2.EdgeColor = 'none';
% h2.LineStyle = '-';