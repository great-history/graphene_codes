%% 该脚本计算的是双层石墨烯四能带模型下的朗道能级图
% 该脚本分为两个部分
% 第一部分：band calculation in continuum model  
%          ref : Ising Superconductivity and Nematicity in Bernal Bilayer Graphene with Strong Spin Orbit Coupling
% 第二部分：estimation of gap by self-consistent calculation

addpath('.\bilayer_hamitonians\')
addpath('.\utils\')
%% parameters set up
gamma0 = 3.16;
gamma1 = 0.381;
gamma3 = 0.380;
gamma4 = 0.140;
delta = 0;
%  u = 0.2 * gamma1; % 0.0 / 0.04 / 0.065 / 0.07 / 0.08
% u的作用：high LL中电荷主要集中在第二层(u > 0)或第一层(u < 0)，而LLL中电荷主要集中在第二层(u > 0)或第一层(u <
% 0) ———— 由此也可推测出写代码时假定了u > 0时第一层电势比第二层低
% 当 | u | 比较小时，high LL在两层上的电荷分布几乎是差不多的，而lowest LL会存在层极化
% 当 | u | 比较大时，high LL 和 lowest LL都会存在层极化， LLL的层极化可能会被改变
% 当 | u | 非常非常大时：
u = 0.1;
lambda = 0.0016;
ext_mat_K_up = u / 2 * diag([-1, -1, 1, 1]) + lambda * diag([1,1,0,0]); % Ising SOC only acts on A1, B1 sites
ext_mat_K_down = u / 2 * diag([-1, -1, 1, 1]) - lambda * diag([1,1,0,0]);

ext_mat_Kp_up = u / 2 * diag([-1, -1, 1, 1]) + lambda * diag([1,1,0,0]);
ext_mat_Kp_down = u / 2 * diag([-1, -1, 1, 1]) - lambda * diag([1,1,0,0]);

%% band calculation
tic
dim_x = 201;
dim_y = 1;
akx_mesh = linspace(-0.10, 0.10, dim_x);
aky_mesh = 0.0 * ones(size(akx_mesh));
eig_enes_K_up = zeros(dim_y, dim_x, 4);
eig_enes_K_down = zeros(dim_y, dim_x, 4);
eig_enes_Kp_up = zeros(dim_y, dim_x, 4);
eig_enes_Kp_down = zeros(dim_y, dim_x, 4);
for i = 1:dim_x
    for j = 1:dim_y
        akx = akx_mesh(j, i);
        aky = aky_mesh(j, i);

        [HK_ham, HKp_ham] = construct_bilayer_continuum_model(gamma0, gamma1, gamma3, gamma4, delta, zeros(4), zeros(4), akx, aky);
        HK_ham_up = HK_ham + ext_mat_K_up;
        HK_ham_down = HK_ham + ext_mat_K_down;
        HKp_ham_up = HKp_ham + ext_mat_Kp_up;
        HKp_ham_down = HKp_ham + ext_mat_Kp_down;
        
        [~, D] = eig(HK_ham_up);
        eig_enes_K_up(j, i, :) = diag(D);
        
        [~, D] = eig(HK_ham_down);
        eig_enes_K_down(j, i, :) = diag(D);
        
        [~, D] = eig(HKp_ham_up);
        eig_enes_Kp_up(j, i, :) = diag(D);
        
        [~, D] = eig(HKp_ham_down);
        eig_enes_Kp_down(j, i, :) = diag(D);
    end
end
toc


figure % 4
center_x = round((1 + dim_x) / 2);
center_y = round((1 + dim_y) / 2);
subplot(1,2,1)
for i = 1:4
    plot(akx_mesh, eig_enes_K_up(center_y, :, i), 'LineWidth', 2)
    hold on;
    plot(akx_mesh, eig_enes_K_down(center_y, :, i), '--', 'LineWidth', 2)
    hold on
end
grid on;
xlim([-0.10,0.10]) % 无量纲量ka在-0.1到0.1之间
ylim([-0.055, -0.040]) % 30meV
xlabel('ka \rightarrow');
ylabel('E(\bf{k}) in meV \rightarrow');
xticks([-0.05 0 0.05]);
xticklabels({'-0.05' '0' '0.05'})
title('effective model E_k versus k_x; for k_y = 0 at K valley')

subplot(1,2,2)
for i = 1:4
    plot(akx_mesh, eig_enes_Kp_up(center_y, :, i), 'LineWidth', 2)
    hold on;
    plot(akx_mesh, eig_enes_Kp_down(center_y, :, i), '--', 'LineWidth', 2)
    hold on
end
grid on;
xlim([-0.10,0.10]) % 无量纲量ka在-0.1到0.1之间
ylim([-0.055, -0.040]) % 30meV
xlabel('ka \rightarrow');
ylabel('E(\bf{k}) in meV \rightarrow');
xticks([-0.05 0 0.05]);
xticklabels({'-0.05' '0' '0.05'})
title('effective model E_k versus k_x; for k_y = 0 at Kp valley')