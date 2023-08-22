% % effective model of trilayer graphene
% gamma0 = 3100;
% gamma1 = 390;
% gamma2 = -28;
% gamma3 = 315;  
% gamma4 = 41;
% gamma5 = 50;
% delta = 46;
% 
% Delta1 = 180; % 0 / 25 / 50 / 150 / 180 / 250
% Delta2 = 0;
% 
% % k_start / k_end 取为 [-0.1,-0.1], [0.1,0.1] 和 取为 [-1,-1], [1,1]
% tic
% err = 0.2;
% re_order = true;
% Nx = 201;
% Ny = 201;
% % Ny = 1;
% 
% % [eig_enes_K, eig_enes_Kp, akxs, akys] = ...
% %     trilayer_effective_band_solver_with_D(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, [-0.125,0], [0.125,0], err, re_order);
% [eig_enes_K, eig_enes_Kp, akxs, akys] = ...
%     trilayer_effective_band_solver_with_D(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, [-0.125,-0.125], [0.125,0.125], err, re_order);
% toc

% figure % 1
% [ene_vecs, dos_vecs] = dos_by_gauss_smear(eig_enes_K, sigma, 201, 201, -2000, 2000, ene_points);
% plot(ene_vecs, dos_vecs, 'b-o')
% hold on 
% [ene_vecs, dos_vecs] = dos_by_gauss_smear(eig_enes_Kp, sigma, 201, 201, -2000, 2000, ene_points);
% plot(ene_vecs, dos_vecs, 'r--')

% [aakkxs, aakkys]=meshgrid(akxs, akys);
% figure % 2
% for i = 1:6
%     mesh(aakkxs,aakkys,eig_enes_K(:,:,i))
%     hold on
% end
% colorbar

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


% for i = 3:4
%     plot(akxs,eig_enes_Kp(:,65,i),'LineWidth', 2)
%     hold on;
%     plot(akxs,eig_enes_Kp(:,70,i),'LineWidth', 2)
%     hold on
%     plot(akxs,eig_enes_Kp(:,75,i),'LineWidth', 2)
%     hold on
%     plot(akxs,eig_enes_Kp(:,85,i),'LineWidth', 2)
%     hold on
%     plot(akxs,eig_enes_Kp(:,90,i),'LineWidth', 2)
%     hold on
%     plot(akxs,eig_enes_Kp(:,95,i),'LineWidth', 2)
%     hold on
%     plot(akxs,eig_enes_Kp(:,101,i),'LineWidth', 2)
%     hold on;
%     plot(akxs,eig_enes_Kp(:,105,i),'LineWidth', 2)
%     hold on
%     plot(akxs,eig_enes_Kp(:,110,i),'LineWidth', 2)
%     hold on
%     plot(akxs,eig_enes_Kp(:,115,i),'LineWidth', 2)
%     hold on
%     plot(akxs,eig_enes_Kp(:,120,i),'LineWidth', 2)
%     hold on
%     plot(akxs,eig_enes_Kp(:,125,i),'LineWidth', 2)
%     hold on
% end
% grid on;
% xlim([-0.125,0.125]) % 无量纲量ka在-0.1到0.1之间
% ylim([-50,50]) % 30meV

test_trilayer_graphene_linecut(-0.10, 0.10, 201);  
function test_trilayer_graphene_linecut(aky_start, aky_end, Ny)
    % effective model of trilayer graphene
    gamma0 = 3100;
    gamma1 = 390;
    gamma2 = -28;
    gamma3 = 315;  
    gamma4 = 41;
    gamma5 = 50;
    delta = 46;

    Delta1 = 180; % 0 / 25 / 50 / 150 / 180 / 250
    Delta2 = 0;

    % k_start / k_end 取为 [-0.1,-0.1], [0.1,0.1] 和 取为 [-1,-1], [1,1]
    tic
    err = 0.005 * gamma0;
    re_order = true;
    Nx = 201;

    [eig_enes_K, eig_enes_Kp, akxs, akys] = ...
        trilayer_effective_band_solver_with_D(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, [-0.125,aky_start], [0.125,aky_end], err, re_order);
    toc
    
    figure
    for idx = 1:Ny
        for i = 1:6
            plot(akxs,eig_enes_K(:,idx,i),'LineWidth', 2)
            hold on;
        end
    end
    grid on;
    xlim([-0.125,0.125]) % 无量纲量ka在-0.1到0.1之间
    ylim([-50,50]) % 30meV
    
    [aakkxs, aakkys]=meshgrid(akxs, akys);
    figure % 2
    for i = 1:6
        mesh(aakkxs,aakkys,eig_enes_K(:,:,i))
        hold on
    end
    colorbar

end