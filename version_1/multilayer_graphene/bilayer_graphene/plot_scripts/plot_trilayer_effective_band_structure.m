% effective model of trilayer graphene
gamma0 = 3100;
gamma1 = 390;
gamma2 = -28;
gamma3 = 315;  
gamma4 = 41;
gamma5 = 50;
delta = 46;

Delta1 = 25; % 0 / 25 / 50 / 150 / 180 / 250
Delta2 = 0;

err = gamma0 * 0.011;
% k_start / k_end 取为 [-0.1,-0.1], [0.1,0.1] 和 取为 [-1,-1], [1,1]
tic
Nx = 201;
Ny = 201;
[eig_enes_K, eig_enes_Kp, akxs, akys] = trilayer_effective_band_solver_with_D(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, [-0.125,-0.125], [0.125,0.125], err, false);
toc


figure % 4
center_x = round((1 + Nx) / 2);
center_y = round((1 + Ny) / 2);
subplot(1,2,1);
for i = 1:6
    plot(akxs,eig_enes_K(:,center_y,i),'LineWidth', 2)
    hold on;
    plot(akys,eig_enes_K(center_x,:,i),'--','LineWidth', 2)
    hold on
end
grid on;
xlim([-0.125,0.125]) % 无量纲量ka在-0.1到0.1之间
ylim([-50,50]) % 30meV
xlabel('ka \rightarrow');
ylabel('E(\bf{k}) in eV \rightarrow');
xticks([-1 -1/2 0 1/2 1]);
xticklabels({'-1' '-1/2' '0' '1/2' '1'})
title('effective model E_k versus k_x; for k_y = 0 at K valley')
hold on

subplot(1,2,2);
for i = 1:6
    plot(akxs,eig_enes_Kp(:,center_y,i),'LineWidth', 2)
    hold on;
    plot(akys,eig_enes_Kp(center_x,:,i),'--','LineWidth', 2)
    hold on
end
grid on;
xlim([-0.125,0.125]) % 无量纲量ka在-0.1到0.1之间
ylim([-50,50]) % 30meV
xlabel('ka \rightarrow');
ylabel('E(\bf{k}) in eV \rightarrow');
xticks([-1 -1/2 0 1/2 1]);
xticklabels({'-1' '-1/2' '0' '1/2' '1'})
title('effective model E_k versus k_x; for k_y = 0 at Kp valley')