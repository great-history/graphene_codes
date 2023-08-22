function [eig_enes, k1s, k2s] = monolayer_band_solver(N_pb, gamma0, m1, m2, k1_start, k2_start)
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
    
    eig_enes = zeros(N1, N2, 2);
    % H_m_k = zeros(2,2);
    for k2_ind = 1:N2
        for k1_ind = 1:N1
            % construct temporary monolayer Hamiltonian for each k = (k1, k2)
            f_k = f_k_mat(k1_ind, k2_ind);
            H_m_k = construct_monolayer_Ham(0, 0, gamma0, f_k);

            % get the eigenvals and eigenvecs
            [eigvecs, eigvals] = eig(H_m_k);
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
% N_pb = 200;
% N1 = N_pb;
% N2 = N_pb;
% 测第一布里渊区
% [eig_enes, k1s, k2s] = monolayer_band_solver(N_pb, gamma0, 1, 1, 0, 0);

% % figure  % 这个图是FBZ为honeycomb的情形，也是graphene最最常用的
% [eig_enes_super, k1s_super, k2s_super] = monolayer_band_solver(N_pb, gamma0, 2, 2, -1, -1);
% 
% % 做一个扩大的布里渊区:
% [b_zone_x, b_zone_y] = convert_k_space(2 * N1, 2 * N2, k1s_super, k2s_super);

% eig_enes_super(1:N1, 1:N2, :) = eig_enes;
% eig_enes_super((N1+1):2*N1, 1:N2, :) = eig_enes;
% eig_enes_super(1:N1, (N2+1):2*N2, :) = eig_enes;
% eig_enes_super((N1+1):2*N1, (N2+1):2*N2, :) = eig_enes;
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
% 
% K_points = zeros(7,2); % 共有6个Dirac points(即K，K_prime)需要存放
% vec_b1 = [1, 1 / sqrt(3)];  % 用直角坐标系表示，省略了系数 2 * pi / a
% vec_b2 = [1, - 1 / sqrt(3)];
% K_points(1,:) = 1 / 3 * vec_b1 + 1 / 3 * vec_b2;
% K_points(4,:) = - 1 / 3 * vec_b1 - 1 / 3 * vec_b2;
% 
% K_points(2,:) = 2 / 3 * vec_b1 - 1 / 3 * vec_b2;
% K_points(5,:) = - 2 / 3 * vec_b1 + 1 / 3 * vec_b2;
% 
% K_points(3,:) = 1 / 3 * vec_b1 - 2 / 3 * vec_b2;
% K_points(6,:) = - 1 / 3 * vec_b1 + 2 / 3 * vec_b2;
% 
% K_points(7,:) = K_points(1,:);
% 
% figure
% pcolor(b_zone_x,b_zone_y,eig_enes_super(:,:,1));
% shading('interp');
% hold on
% contour(b_zone_x,b_zone_y,eig_enes_super(:,:,1),'LineColor','b');
% hold on
% plot(K_points(:,1),K_points(:,2),'b--o', 'MarkerSize',10, 'MarkerEdgeColor','b', 'MarkerFaceColor',[0.5,0.5,0.5])
% colorbar
% xlabel('k_x \rightarrow');
% % xticks([-6*pi/3 -4*pi/3 -2*pi/3 0 2*pi/3 4*pi/3 6*pi/3]);
% % xticklabels({'-6\pi/3' '-4\pi/3' '-2\pi/3' '0' '2\pi/3' '4*\pi/3' '6*\pi/3'})
% % 
% % ylabel('k_y \rightarrow');
% % yticks([-6*pi/3 -4*pi/3 -2*pi/3 0 2*pi/3 4*pi/3 6*pi/3]);
% % yticklabels({'-6\pi/3' '-4\pi/3' '-2\pi/3' '0' '2\pi/3' '4*\pi/3' '6*\pi/3'})
% 
% zlabel('E(\bf{k}) in eV \rightarrow')
% 
% title('Contour plot of low-energy band: Numerical')
% legend('Valence band', 'Contour');
% 
% figure
% pcolor(b_zone_x,b_zone_y,eig_enes_super(:,:,2));
% shading('interp');
% hold on
% contour(b_zone_x,b_zone_y,eig_enes_super(:,:,2),'LineColor','b');
% colorbar
% xlabel('k_x \rightarrow');
% % xticks([-6*pi/3 -4*pi/3 -2*pi/3 0 2*pi/3 4*pi/3 6*pi/3]);
% % xticklabels({'-6\pi/3' '-4\pi/3' '-2\pi/3' '0' '2\pi/3' '4*\pi/3' '6*\pi/3'})
% % 
% % ylabel('k_y \rightarrow');
% % yticks([-6*pi/3 -4*pi/3 -2*pi/3 0 2*pi/3 4*pi/3 6*pi/3]);
% % yticklabels({'-6\pi/3' '-4\pi/3' '-2\pi/3' '0' '2\pi/3' '4*\pi/3' '6*\pi/3'})
% 
% zlabel('E(\bf{k}) in eV \rightarrow')
% 
% title('Contour plot of high-energy band: Numerical')
% legend('Conductance band', 'Contour');