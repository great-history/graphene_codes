function H_t_k = construct_trilayer_Ham_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, f_k)
    H_t_k = zeros(6,6);

    H_m_k = construct_monolayer_Ham(Delta2 - gamma2 / 2, Delta2 + delta - gamma5 / 2, gamma0, f_k);
    on_sites = [Delta2 + gamma2 / 2, -2*Delta2, -2*Delta2 + delta, Delta2 + delta + gamma5 /2];
    H_b_k = construct_bilayer_Ham(on_sites, gamma0, sqrt(2) * gamma1, sqrt(2) * gamma3, sqrt(2) * gamma4, f_k);
    
    H_t_k(1:2,1:2) = H_m_k;
    H_t_k(3:6,3:6) = H_b_k;

    H_t_k(1,3) = Delta1;
    H_t_k(3,1) = Delta1;
    H_t_k(2,4) = Delta1;
    H_t_k(4,2) = Delta1;
    
    hem1 = helper_check_hermite(H_t_k,1e-8);
    if hem1 == 0
        disp("the Hamiltonian is not hermitian")
    end
end


function H_b_k = construct_bilayer_Ham(on_sites, gamma0, sqrt2_gamma1, sqrt2_gamma3, sqrt2_gamma4, f_k)
    H_b_k = zeros(4,4);
    H_b_k(1,1) = on_sites(1);
    H_b_k(2,2) = on_sites(2);
    H_b_k(3,3) = on_sites(3);
    H_b_k(4,4) = on_sites(4);
    
    H_b_k(3,4) = sqrt2_gamma1;
    H_b_k(4,3) = sqrt2_gamma1;

    x0 = gamma0 * f_k;
    x3 = sqrt2_gamma3 * f_k;
    x4 = sqrt2_gamma4 * f_k;
    
    H_b_k(1,4) = - x0;
    H_b_k(4,1) = - conj(x0);
    H_b_k(2,3) = - conj(x0);
    H_b_k(3,2) = - x0;
    
    H_b_k(1,2) = - conj(x3);
    H_b_k(2,1) = - x3;
    
    H_b_k(1,3) = x4;
    H_b_k(3,1) = conj(x4);
    H_b_k(2,4) = conj(x4);
    H_b_k(4,2) = x4;
    
    hem = helper_check_hermite(H_b_k, 1e-8);
    if hem == 0
       disp("monolayer Ham is not hermitian")
    end
end

function H_t_k = construct_trilayer_Ham_with_D_old(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, f_k)
    % 这种构造哈密顿量的方式是原始的，没有进行幺正变换前的样子
    H_t_k = zeros(6,6);

    H_m_k = construct_monolayer_Ham(0, 0, gamma0, f_k);
    
    V_k = zeros(2,2);
    V_k(1,1) = gamma4 * f_k;
    V_k(1,2) = - gamma3 * conj(f_k);
    V_k(2,1) = gamma1;
    V_k(2,2) = gamma4 * f_k;

    W_k = zeros(2,2);
    W_k(1,1) = gamma2 / 2;
    W_k(2,2) = gamma5 / 2;

    H_t_k(1:2, 1:2) = H_m_k;
    H_t_k(3:4, 3:4) = H_m_k;
    H_t_k(5:6, 5:6) = H_m_k;
    H_t_k(1:2, 3:4) = V_k;
    H_t_k(3:4, 1:2) = V_k';
    H_t_k(3:4, 5:6) = V_k';
    H_t_k(5:6, 3:4) = V_k;
    H_t_k(1:2, 5:6) = W_k;
    H_t_k(5:6, 1:2) = W_k';

    H_delta = diag([0, delta, delta, 0, 0, delta]);
    H_Delta1 = diag([Delta1, Delta1, 0, 0, -Delta1, -Delta1]);
    H_Delta2 = diag([Delta2, Delta2, -2 * Delta2, -2*Delta2, Delta2, Delta2]);

    H_t_k = H_t_k + H_delta + H_Delta1 + H_Delta2;

    hem1 = check_hermite(H_t_k,1e-8);
    if hem1 == 0
        disp("the Hamiltonian is not hermitian")
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
% 
% %% plot DOS
% figure
% histogram(tri_eig_enes, 200)
% hold on
% 
% figure
% % 使用gauss smear
% sigma = 30;
% ene_points = 1000;
% ene_lb = min(min(min(tri_eig_enes)));
% ene_ub = max(max(max(tri_eig_enes)));
% [ene_vecs, dos_vecs] = dos_by_gauss_smear(tri_eig_enes, sigma, N1, N2, ene_lb+1000, ene_ub-2000, ene_points);
% plot(ene_vecs, dos_vecs)
% 
% % figure  % 这个图是FBZ为平行四边形的情形
% figure
% for i = 1:6
%     mesh(kxs,kys,tri_eig_enes(:,:,i))
%     hold on
% end
% colorbar