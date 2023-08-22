%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算能带结构 附加 计算每个本征态的单层/双层的成分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% effective model of trilayer graphene
% 0GPa
% gamma0 = 3100;
% gamma1 = 384.5;
% gamma2 = -19.4;
% gamma3 = 351;  
% gamma4 = 65.3;
% gamma5 = 63.6;
% delta = 41.5;
% 
% Delta1 = 0; % 0 / 25 / 50 / 150 / 180 / 250
% Delta2 = 2.1;

% 1GPa
% gamma0 = 3100;
% gamma1 = 378.4;
% gamma2 = -23.4;
% gamma3 = 366.5;  
% gamma4 = 61.2;
% gamma5 = 83.3;
% delta = 53.4;
% 
% Delta1 = 0; % 0 / 25 / 50 / 150 / 180 / 250
% Delta2 = 3.1;

gamma0 = 3100;
gamma1 = 390;
gamma2 = -28;
gamma3 = 315;  
gamma4 = 41;
gamma5 = 50;
delta = 46;

Delta1 = 25; % 0 / 25 / 50 / 150 / 180 / 250
Delta2 = 3.1;

% k_start / k_end 取为 [-0.1,-0.1], [0.1,0.1] 和 取为 [-1,-1], [1,1]
tic
Nx = 201;
Ny = 201;
ak_start = [-0.125,-0.125];
ak_end = [0.125,0.125];

% [eig_enes_K, eig_enes_Kp, akxs, akys] = trilayer_effective_band_solver_with_D(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, [-0.125,-0.125], [0.125,0.125]);
toc

akx_list = linspace(ak_start(1), ak_end(1), Nx);
aky_list = linspace(ak_start(2), ak_end(2), Ny);

% 存放本征值
eig_enes_K = zeros(Ny, Nx, 6);
eig_enes_Kp = zeros(Ny, Nx, 6);
red_RGB = [1 0 0];
blue_RGB = [0 0 1];
eig_clc_K = zeros(Ny, Nx, 6, 3); % 每个本征态对应的颜色 red [] 对应 单层 blue 对应 双层 []
eig_clc_Kp = zeros(Ny, Nx, 6, 3);

for i = 1:Nx
    akx = akx_list(i);
    for j = 1:Ny
        aky = aky_list(j);

        [HK_m_ham, HKp_m_ham, HK_b_ham, HKp_b_ham] = construct_trilayer_ABA_six_bands_continuum_model(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, akx, aky);
        if Delta1 == 0
            [~, D] = eig(HK_m_ham);
            eig_enes_K(j, i, 1:2) = diag(D);
            
            [~, D] = eig(HK_b_ham);
            eig_enes_K(j, i, 3:6) = diag(D);
            
            eig_clc_K(j, i, 1:2, :) = [red_RGB; red_RGB];
            eig_clc_K(j, i, 3:6, :) = [blue_RGB; blue_RGB; blue_RGB; blue_RGB];
            
            [~, D] = eig(HKp_m_ham);
            eig_enes_Kp(j, i, 1:2) = diag(D);

            [~, D] = eig(HKp_b_ham);
            eig_enes_Kp(j, i, 3:6) = diag(D);
            
            eig_clc_Kp(j, i, 1:2, :) = [red_RGB; red_RGB];
            eig_clc_Kp(j, i, 3:6, :) = [blue_RGB; blue_RGB; blue_RGB; blue_RGB];
            
        else
            Delta1_mat = zeros(2, 4);
            Delta1_mat(1, 1) = Delta1;
            Delta1_mat(2, 2) = Delta1;

            HK_ham = [HK_m_ham, Delta1_mat; Delta1_mat', HK_b_ham];
            HKp_ham = [HKp_m_ham, Delta1_mat; Delta1_mat', HKp_b_ham];

            [V, D] = eig(HK_ham);
            eig_enes_K(j, i, :) = diag(D);
            
            for k = 1:6
                red_weight = sum(abs(V(1:2, k)).^2);
                blue_weight = sum(abs(V(3:6, k)).^2);
                eig_clc_K(j, i, k, :) = [red_RGB * red_weight + blue_RGB * blue_weight];
            end
            
            [V, D] = eig(HKp_ham);
            eig_enes_Kp(j, i, :) = diag(D);
            
            for k = 1:6
                red_weight = sum(abs(V(1:2, k)).^2);
                blue_weight = sum(abs(V(3:6, k)).^2);
                eig_clc_Kp(j, i, k, :) = [red_RGB * red_weight + blue_RGB * blue_weight];
            end
        end

    end
end


figure % 4
set(gca,'FontName','Times New Roman','FontSize',7);

center_x = round((1 + Nx) / 2);
center_y = round((1 + Ny) / 2);
select_direction = "y"; % it can be "x" or "y"
select_valley = "K"; % it can be "K" or "Kp"

%示例代码
switch [select_direction, select_valley]
    case ["x", "K"]
        if Delta1 == 0
             % 单层
            for i = 1:2
                plot(akx_list, eig_enes_K(center_y, :, i),'Color', red_RGB, 'LineWidth', 2)
                hold on
            end

            % 双层
            for i = 3:6
                plot(akx_list, eig_enes_K(center_y, :, i),'Color', blue_RGB, 'LineWidth', 2)
                hold on
            end
        else % Delta1不等于0
            % 单层
            for i = 1:6
                for j = 1:(Nx - 1)
                   plot(akx_list(j:j+1), eig_enes_K(center_y, j:j+1, i),'Color', 1 / 2 * (eig_clc_K(center_y, j, i, :) + eig_clc_K(center_y, j, i, :)), 'LineWidth', 2)
                   hold on;
                end
            end
        end
    case ["y", "K"]
        if Delta1 == 0
             % 单层
            for i = 1:2
                plot(akx_list, eig_enes_K(:, center_x, i),'Color', red_RGB, 'LineWidth', 2)
                hold on
            end

            % 双层
            for i = 3:6
                plot(akx_list, eig_enes_K(:, center_x, i),'Color', blue_RGB, 'LineWidth', 2)
                hold on
            end
        else % Delta1不等于0
            % 单层
            for i = 1:6
                for j = 1:(Nx - 1)
                   plot(aky_list(j:j+1), eig_enes_K(j:j+1, center_x, i),'Color', 1 / 2 * (eig_clc_K(j, center_x, i, :) + eig_clc_K(j, center_x, i, :)), 'LineWidth', 2)
                   hold on;
                end
            end
        end
    case ["x", "Kp"]
        if Delta1 == 0
             % 单层
            for i = 1:2
                plot(akx_list, eig_enes_Kp(center_y, :, i),'Color', red_RGB, 'LineWidth', 2)
                hold on
            end

            % 双层
            for i = 3:6
                plot(akx_list, eig_enes_Kp(center_y, :, i),'Color', blue_RGB, 'LineWidth', 2)
                hold on
            end
        else % Delta1不等于0
            % 单层
            for i = 1:6
                for j = 1:(Nx - 1)
                   plot(akx_list(j:j+1), eig_enes_Kp(center_y, j:j+1, i),'Color', 1 / 2 * (eig_clc_Kp(center_y, j, i, :) + eig_clc_Kp(center_y, j, i, :)), 'LineWidth', 2)
                   hold on;
                end
            end
        end
    case ["y", "Kp"]
        if Delta1 == 0
             % 单层
            for i = 1:2
                plot(akx_list, eig_enes_Kp(:, center_x, i),'Color', red_RGB, 'LineWidth', 2)
                hold on
            end

            % 双层
            for i = 3:6
                plot(akx_list, eig_enes_Kp(:, center_y, i),'Color', blue_RGB, 'LineWidth', 2)
                hold on
            end
        else % Delta1不等于0
            % 单层
            for i = 1:6
                for j = 1:(Nx - 1)
                   plot(aky_list(j:j+1), eig_enes_Kp(j:j+1, center_y, i),'Color', 1 / 2 * (eig_clc_Kp(j, center_y, i, :) + eig_clc_Kp(j, center_y, i, :)), 'LineWidth', 2)
                   hold on;
                end
            end
        end
    otherwise
        disp('No such type !!!')
end

% % 单层
% for i = 1:2
%     plot(akx_list, eig_enes_K(:,center_y, i),'r', 'LineWidth', 2)
%     hold on;
%     % plot(aky_list, eig_enes_Kp(:,center_y, i),'r--','LineWidth', 2)
%     % plot(aky_list, eig_enes_K(center_x, :, i),'r--','LineWidth', 2)
%     hold on
% end
% 
% % 双层
% for i = 3:6
%     plot(akx_list, eig_enes_K(:,center_y, i), 'b', 'LineWidth', 2)
%     hold on;
%     % plot(aky_list, eig_enes_Kp(:,center_y,i), 'b--', 'LineWidth', 2)
%     % plot(aky_list, eig_enes_K(center_x, :,i), 'b--', 'LineWidth', 2)
%     hold on
% end

grid off;
xlim([-0.125,0.125]) % 无量纲量ka在-0.1到0.1之间
ylim([-50,50]) % 30meV
% xlabel('ka \rightarrow');
% ylabel('E(\bf{k}) in eV \rightarrow');
% xticks([-1 -1/2 0 1/2 1]);
% xticklabels({'-1' '-1/2' '0' '1/2' '1'})
% title('effective model E_k versus k_x; for k_y = 0 at K valley')