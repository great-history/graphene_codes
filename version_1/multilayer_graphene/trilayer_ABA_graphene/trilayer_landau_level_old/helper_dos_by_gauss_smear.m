function [ene_vecs, dos_vecs] = dos_by_gauss_smear(eig_enes, sigma, N1, N2, ene_lb, ene_ub, ene_points)
    % sigma的值 50 20, 10, 5 值越小越不光滑，值越大越光滑，但不能太大
    ene_vecs = linspace(ene_lb, ene_ub, ene_points);
    dos_vecs = zeros(ene_points,1);
    
    size_eig_enes = size(eig_enes);
    for i = 1:ene_points
        val = 0.0;
        for jj = 1:N1
            for kk = 1:N2
                for ii = 1:size_eig_enes(3)
                    val = val + 1 / (sigma * sqrt(2 * pi)) * exp(-(ene_vecs(i) - eig_enes(jj, kk, ii))^2 / (2 * sigma^2));
                end
            end
        end
        dos_vecs(i) = val / (N1 * N2);
    end
end

% % example of trilayer plot
% gamma0 = 3100;
% gamma1 = 355;
% gamma2 = -20.7904;
% gamma3 = 346.7153;  
% gamma4 = 91.2578;
% gamma5 = 35.7223;
% 
% D_field = 1;
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
% %% plot DOS
% figure
% % histogram(eig_enes, 200)
% % hold on
% % 使用gauss smear
% sigma = 50;
% ene_points = 1000;
% ene_lb = min(min(min(tri_eig_enes)));
% ene_ub = max(max(max(tri_eig_enes)));
% [ene_vecs, dos_vecs] = dos_by_gauss_smear(tri_eig_enes, sigma, N1, N2, ene_lb+1000, ene_ub-2000, ene_points);
% plot(ene_vecs, dos_vecs)