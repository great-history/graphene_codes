% parameters adopted from the PRB paper:<<New Dirac points and multiple Landau level crossings in biased trilayer graphene>>
gamma0 = 3100;
gamma1 = 390;
gamma2 = -28;
gamma3 = 315;  
gamma4 = 41;
gamma5 = 50;
delta = 46;

Delta1 = 25; % 0 / 25 / 50 / 150 / 180 / 250
Delta2 = 0;

N_pb = 200;
N1 = N_pb;
N2 = N_pb;
err = 200;  % 这里的error一定要设的足够大才行
[tri_eig_enes, k1s, k2s] = trilayer_band_solver_with_D(N_pb, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, 1, 1, 0, 0);

% 转到直角坐标系
[kxs, kys] = helper_convert_k_space(N1, N2, k1s, k2s);

% plot DOS
% figure % 1
% histogram(tri_eig_enes, 200)

figure % 2
% 使用gauss smear
sigma = 30;
ene_points = 1000;
ene_lb = min(min(min(tri_eig_enes)));
ene_ub = max(max(max(tri_eig_enes)));
[ene_vecs, dos_vecs] = helper_dos_by_gauss_smear(tri_eig_enes, sigma, N1, N2, ene_lb+1000, ene_ub-2000, ene_points);
plot(ene_vecs, dos_vecs)

figure % 4
for i = 1:6
    mesh(kxs,kys,tri_eig_enes(:,:,i))
    hold on
end
colorbar



% % effective model of trilayer graphene
% gamma0 = 3100;
% gamma1 = 390;
% gamma2 = -28;
% gamma3 = 315;  
% gamma4 = 41;
% gamma5 = 50;
% delta = 46;
% 
% Delta1 = 25; % 0 / 25 / 50 / 150 / 180 / 250
% Delta2 = 0;
% 
% % k_start / k_end 取为 [-0.1,-0.1], [0.1,0.1] 和 取为 [-1,-1], [1,1]
% tic
% Nx = 201;
% Ny = 201;
% [eig_enes_K, eig_enes_Kp, akxs, akys] = trilayer_effective_band_solver_with_D(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, [-0.125,-0.125], [0.125,0.125], err, true);
% toc
% 
% figure % 1
% [ene_vecs, dos_vecs] = helper_dos_by_gauss_smear(eig_enes_K, sigma, 201, 201, -2000, 2000, ene_points);
% plot(ene_vecs, dos_vecs, 'b-o')
% hold on 
% [ene_vecs, dos_vecs] = helper_dos_by_gauss_smear(eig_enes_Kp, sigma, 201, 201, -2000, 2000, ene_points);
% plot(ene_vecs, dos_vecs, 'r--')
% 
% 
% [aakkxs, aakkys]=meshgrid(akxs, akys);
% figure % 2
% for i = 1:6
%     mesh(aakkxs,aakkys,eig_enes_K(:,:,i))
%     hold on
% end
% colorbar
% 
% figure % 3
% for i = 1:2
%     for j = 1:3
%         index = (i-1) * 3 + j;
%         subplot(2,3,index)
%         
%         pcolor(aakkxs,aakkys,eig_enes_K(:,:,index));
%         shading('interp');
%         hold on
%         contour(aakkxs,aakkys,eig_enes_K(:,:,index),'LineColor','b');
%         hold on
%         contour(aakkxs,aakkys,eig_enes_Kp(:,:,index),'--','LineColor','r');
%         hold on
%         colorbar
%         xlabel('k_x \rightarrow');
% 
%         zlabel('E(\bf{k}) in eV \rightarrow')
% 
%         title('Contour plot of low-energy band: Numerical')
%     end
% end
% 
% figure % 4
% center_x = round((1 + Nx) / 2);
% center_y = round((1 + Ny) / 2);
% subplot(1,2,1);
% for i = 1:6
%     plot(akxs,eig_enes_K(:,center_y,i),'LineWidth', 2)
%     hold on;
%     plot(akys,eig_enes_K(center_x,:,i),'--','LineWidth', 2)
%     hold on
% end
% grid on;
% xlim([-0.125,0.125]) % 无量纲量ka在-0.1到0.1之间
% ylim([-50,50]) % 30meV
% xlabel('ka \rightarrow');
% ylabel('E(\bf{k}) in eV \rightarrow');
% xticks([-1 -1/2 0 1/2 1]);
% xticklabels({'-1' '-1/2' '0' '1/2' '1'})
% title('effective model E_k versus k_x; for k_y = 0 at K valley')
% hold on
% 
% subplot(1,2,2);
% for i = 1:6
%     plot(akxs,eig_enes_Kp(:,center_y,i),'LineWidth', 2)
%     hold on;
%     plot(akys,eig_enes_Kp(center_x,:,i),'--','LineWidth', 2)
%     hold on
% end
% grid on;
% xlim([-0.125,0.125]) % 无量纲量ka在-0.1到0.1之间
% ylim([-50,50]) % 30meV
% xlabel('ka \rightarrow');
% ylabel('E(\bf{k}) in eV \rightarrow');
% xticks([-1 -1/2 0 1/2 1]);
% xticklabels({'-1' '-1/2' '0' '1/2' '1'})
% title('effective model E_k versus k_x; for k_y = 0 at Kp valley')