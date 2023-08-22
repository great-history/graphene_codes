gamma0 = 3.1;
% N_pb = 200;
% N1 = N_pb;
% N2 = N_pb;
% [eig_enes, k1s, k2s] = monolayer_band_solver(N_pb, gamma0, 1, 1, 0, 0);
Nx = 201;
Ny = 201;
err = 0.2;
re_order = true;
% 画图时需要注意的一点是：
% 如果存在能隙闭合的情况(尤其是对于Dirac Point这种奇异点)，那么画线型和画三维有点不一样
% 画线性时(比如graphene band structure vs k / Landau level vs B)需要把re-order打开，因为这样才能把每条能带的演化轨迹看清楚
% 画三维时就不需要把re-order打开了，如果打开，那么在一些奇异点附近会出现问题：能带上下震荡

[eig_enes_K, eig_enes_Kp, akxs, akys] = monolayer_effective_band_solver(Nx, Ny, 0, 0, gamma0, [-0.1, -0.1], [0.1, 0.1], err, re_order);
% 转到直角坐标系
% [kxs, kys] = convert_k_space(N1, N2, k1s, k2s);


%% plot DOS of monolayer graphene
figure % 1
% histogram(eig_enes_K, 200)
% hold on

% 使用gauss smear
sigma = gamma0 * 0.0005;
ene_points = 1000;
[ene_vecs, dos_vecs] = helper_dos_by_gauss_smear(eig_enes_K, sigma, Nx, Ny, -gamma0 * 0.05, gamma0 * 0.05, ene_points);
plot(ene_vecs, dos_vecs, 'b--o')
hold on 
[ene_vecs, dos_vecs] = helper_dos_by_gauss_smear(eig_enes_Kp, sigma, Nx, Ny, -gamma0 * 0.05, gamma0 * 0.05, ene_points);
plot(ene_vecs, dos_vecs, 'r-*')

% % plot band structure of monolayer graphene using effective model
figure % 2
subplot(1,2,1);
for i = 1:Nx
    plot(akxs,eig_enes_K(i,:,1),'LineWidth', 2)
    hold on;
    plot(akxs,eig_enes_K(i,:,2),'LineWidth', 2)
    hold on;
end
grid on;
xlabel('k_x \rightarrow');
ylabel('E(\bf{k}) in eV \rightarrow');
xticks([-0.1 -0.05 0 0.05 0.1]);
xticklabels({'-0.1' '-0.05' '0' '0.05' '0.1'})
title('effective model E_k versus k_x; for k_y = 0 at K valley')
legend('Valence band', 'Conduction band');
hold on

subplot(1,2,2);
plot(akxs,eig_enes_Kp(101,:,1),'LineWidth', 2)  % 要把re-order关掉
grid on;
hold on;
plot(akxs,eig_enes_Kp(101,:,2),'LineWidth', 2)
xlabel('k_x \rightarrow');
ylabel('E(\bf{k}) in eV \rightarrow');
xticks([-1 -1/2 0 1/2 1]);
xticklabels({'-1' '-1/2' '0' '1/2' '1'})
title('effective model E_k versus k_x; for k_y = 0 at Kp valley')
legend('Valence band', 'Conduction band');

[aakkxs, aakkys]=meshgrid(akxs, akys);

figure % 4
mesh(aakkxs,aakkys,eig_enes_K(:,:,1))
hold on
mesh(aakkxs,aakkys,eig_enes_K(:,:,2))
colorbar

figure % 5
mesh(akxs,akys, eig_enes_Kp(:,:,1))
hold on
mesh(akxs,akys, eig_enes_Kp(:,:,2))
colorbar