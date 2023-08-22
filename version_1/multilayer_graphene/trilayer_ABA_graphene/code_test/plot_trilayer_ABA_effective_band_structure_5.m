%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 运行前请先运行qhe_LL_plot.m 以获得跃迁参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算能带结构 附加 计算每个本征态的层分布
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma0 = hopping_params_after(1) * 1000; % 转化为meV单位
gamma1 = hopping_params_after(2) * 1000;
gamma2 = hopping_params_after(3) * 1000;
gamma3 = hopping_params_after(4) * 1000;  
gamma4 = hopping_params_after(5) * 1000;
gamma5 = hopping_params_after(6) * 1000;
delta = hopping_params_after(7) * 1000;
Delta2 = hopping_params_after(8) * 1000;

Delta1 = 100; % 0 / 25 / 50 / 150 / 180 / 250

% k_start / k_end 取为 [-0.1,-0.1], [0.1,0.1] 和 取为 [-1,-1], [1,1]
tic
Nx = 101;
Ny = 101;
ak_start = [-0.1, -0.1];
ak_end = [0.1, 0.1];
% ak_start = [-0.125,-0.125];
% ak_end = [0.125,0.125];

% [eig_enes_K, eig_enes_Kp, akxs, akys] = trilayer_effective_band_solver_with_D(Nx, Ny, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, [-0.125,-0.125], [0.125,0.125]);
toc

akx_list = linspace(ak_start(1), ak_end(1), Nx);
aky_list = linspace(ak_start(2), ak_end(2), Ny);

% 存放本征值
eig_enes_K = zeros(Ny, Nx, 6);
eig_enes_Kp = zeros(Ny, Nx, 6);
% 颜色
red_RGB = [1 0 0];
blue_RGB = [0 0 1];
eig_clc_K = zeros(Ny, Nx, 6, 3); % 每个本征态对应的颜色 red [] 对应 单层 blue 对应 双层 []
eig_clc_Kp = zeros(Ny, Nx, 6, 3);
% 直方图
eig_layer_K = zeros(Ny, Nx, 6, 3); % 每个本征态对应的层分布
eig_layer_Kp = zeros(Ny, Nx, 6, 3);

for i = 1:Nx
    akx = akx_list(i);
    for j = 1:Ny
        aky = aky_list(j);

        [HK_m_ham, HKp_m_ham, HK_b_ham, HKp_b_ham] = construct_trilayer_ABA_six_bands_continuum_model(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, akx, aky);
        if Delta1 == 0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % K valley 的计算
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Vm_K, D] = eig(HK_m_ham);
            eig_enes_K(j, i, 1:2) = diag(D);
            
            [Vb_K, D] = eig(HK_b_ham);
            eig_enes_K(j, i, 3:6) = diag(D);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 计算颜色
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            eig_clc_K(j, i, 1:2, :) = [red_RGB; red_RGB];
            eig_clc_K(j, i, 3:6, :) = [blue_RGB; blue_RGB; blue_RGB; blue_RGB];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 计算层分布
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            V_K = [Vm_K, zeros(2,4);zeros(4,2), Vb_K];
            % layer 1 sites
            % prob_mK_layer1_mat = transpose(Vm_K) * [1, 0; 0, 1] * 1 / sqrt(2);
            % prob_bK_layer1_mat = transpose(Vb_K) * [1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            prob_K_layer1_mat = transpose(V_K) * [1, 0; 0, 1; 1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            
            % layer 2 sites
            % prob_bK_layer2_mat = transpose(Vb_K) * [0, 0; 0, 0; 1, 0; 0, 1];
            prob_K_layer2_mat = transpose(V_K) * [0, 0; 0, 0; 0, 0; 0, 0; 1, 0; 0, 1];
            
            % layer 3 sites
            % prob_mK_layer3_mat = - transpose(Vm_K) * [1, 0; 0, 1] * 1 / sqrt(2);
            % prob_bK_layer3_mat = transpose(Vb_K) * [1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            prob_K_layer3_mat = transpose(V_K) * [-1, 0; 0, -1; 1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            
            % 存入矩阵中
            for k = 1:6
                % layer 1
                eig_layer_K(j, i, k, 1) = (abs(prob_K_layer1_mat(k, 1)))^2 + (abs(prob_K_layer1_mat(k, 2)))^2;
                % layer 2
                eig_layer_K(j, i, k, 2) = (abs(prob_K_layer2_mat(k, 1)))^2 + (abs(prob_K_layer2_mat(k, 2)))^2;
                % layer 3
                eig_layer_K(j, i, k, 3) = (abs(prob_K_layer3_mat(k, 1)))^2 + (abs(prob_K_layer3_mat(k, 2)))^2;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kp valley 的计算
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Vm_Kp, D] = eig(HKp_m_ham);
            eig_enes_Kp(j, i, 1:2) = diag(D);

            [Vb_Kp, D] = eig(HKp_b_ham);
            
            eig_enes_Kp(j, i, 3:6) = diag(D);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 计算颜色
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            eig_clc_Kp(j, i, 1:2, :) = [red_RGB; red_RGB];
            eig_clc_Kp(j, i, 3:6, :) = [blue_RGB; blue_RGB; blue_RGB; blue_RGB];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 计算层分布
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            V_Kp = [Vm_Kp, zeros(2,4);zeros(4,2), Vb_Kp];
            % layer 1 sites
            % prob_mKp_layer1_mat = transpose(Vm_Kp) * [1, 0; 0, 1] * 1 / sqrt(2);
            % prob_bKp_layer1_mat = transpose(Vb_Kp) * [1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            prob_Kp_layer1_mat = transpose(V_Kp) * [1, 0; 0, 1; 1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            
            % layer 2 sites
            % prob_bKp_layer2_mat = transpose(Vb_Kp) * [0, 0; 0, 0; 1, 0; 0, 1];
            prob_Kp_layer2_mat = transpose(V_Kp) * [0, 0; 0, 0; 0, 0; 0, 0; 1, 0; 0, 1];
            
            % layer 3 sites
            % prob_mKp_layer3_mat = - transpose(Vm_Kp) * [1, 0; 0, 1] * 1 / sqrt(2);
            % prob_bKp_layer3_mat = transpose(Vb_Kp) * [1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            prob_Kp_layer3_mat = transpose(V_Kp) * [-1, 0; 0, -1; 1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            
            % 存入矩阵中
            for k = 1:6
                % layer 1
                eig_layer_Kp(j, i, k, 1) = (abs(prob_Kp_layer1_mat(k, 1)))^2 + (abs(prob_Kp_layer1_mat(k, 2)))^2;
                % layer 2
                eig_layer_Kp(j, i, k, 2) = (abs(prob_Kp_layer2_mat(k, 1)))^2 + (abs(prob_Kp_layer2_mat(k, 2)))^2;
                % layer 3
                eig_layer_Kp(j, i, k, 3) = (abs(prob_Kp_layer3_mat(k, 1)))^2 + (abs(prob_Kp_layer3_mat(k, 2)))^2;
            end
        else
            Delta1_mat = zeros(2, 4);
            Delta1_mat(1, 1) = Delta1;
            Delta1_mat(2, 2) = Delta1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % K valley 的计算
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            HK_ham = [HK_m_ham, Delta1_mat; Delta1_mat', HK_b_ham];
            HKp_ham = [HKp_m_ham, Delta1_mat; Delta1_mat', HKp_b_ham];

            [V_K, D] = eig(HK_ham);
            eig_enes_K(j, i, :) = diag(D);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 计算颜色
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:6
                red_weight = sum(abs(V_K(1:2, k)).^2);
                blue_weight = sum(abs(V_K(3:6, k)).^2);
                eig_clc_K(j, i, k, :) = [red_RGB * red_weight + blue_RGB * blue_weight];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 计算层分布
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % layer 1 sites
            prob_K_layer1_mat = transpose(V_K) * [1, 0; 0, 1; 1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            
            % layer 2 sites
            prob_K_layer2_mat = transpose(V_K) * [0, 0; 0, 0; 0, 0; 0, 0; 1, 0; 0, 1];
            
            % layer 3 sites
            prob_K_layer3_mat = transpose(V_K) * [-1, 0; 0, -1; 1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            
            % 存入矩阵中
            for k = 1:6
                % layer 1
                eig_layer_K(j, i, k, 1) = (abs(prob_K_layer1_mat(k, 1)))^2 + (abs(prob_K_layer1_mat(k, 2)))^2;
                % layer 2
                eig_layer_K(j, i, k, 2) = (abs(prob_K_layer2_mat(k, 1)))^2 + (abs(prob_K_layer2_mat(k, 2)))^2;
                % layer 3
                eig_layer_K(j, i, k, 3) = (abs(prob_K_layer3_mat(k, 1)))^2 + (abs(prob_K_layer3_mat(k, 2)))^2;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kp valley 的计算
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [V_Kp, D] = eig(HKp_ham);
            eig_enes_Kp(j, i, :) = diag(D);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 计算颜色
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:6
                red_weight = sum(abs(V_Kp(1:2, k)).^2);
                blue_weight = sum(abs(V_Kp(3:6, k)).^2);
                eig_clc_Kp(j, i, k, :) = [red_RGB * red_weight + blue_RGB * blue_weight];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 计算层分布
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % layer 1 sites
            prob_Kp_layer1_mat = transpose(V_Kp) * [1, 0; 0, 1; 1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            
            % layer 2 sites
            prob_Kp_layer2_mat = transpose(V_Kp) * [0, 0; 0, 0; 0, 0; 0, 0; 1, 0; 0, 1];
            
            % layer 3 sites
            prob_Kp_layer3_mat = transpose(V_Kp) * [-1, 0; 0, -1; 1, 0; 0, 1; 0, 0; 0, 0] * 1 / sqrt(2);
            
            % 存入矩阵中
            for k = 1:6
                % layer 1
                eig_layer_Kp(j, i, k, 1) = (abs(prob_Kp_layer1_mat(k, 1)))^2 + (abs(prob_Kp_layer1_mat(k, 2)))^2;
                % layer 2
                eig_layer_Kp(j, i, k, 2) = (abs(prob_Kp_layer2_mat(k, 1)))^2 + (abs(prob_Kp_layer2_mat(k, 2)))^2;
                % layer 3
                eig_layer_Kp(j, i, k, 3) = (abs(prob_Kp_layer3_mat(k, 1)))^2 + (abs(prob_Kp_layer3_mat(k, 2)))^2;
            end
        end

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 三维图 (仅包含双层)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hole side
eig_index = 3;
figure
set(gca,'FontName','Times New Roman','FontSize',7);
[akx_mesh, aky_mesh] = meshgrid(akx_list, aky_list);
surf(akx_mesh, aky_mesh, eig_enes_K(:,:,eig_index), squeeze(eig_layer_K(:,:,eig_index,:)), 'FaceAlpha',1, 'EdgeColor', 'interp');
hold on
% electron side
eig_index = 4;
set(gca,'FontName','Times New Roman','FontSize',7);
[akx_mesh, aky_mesh] = meshgrid(akx_list, aky_list);
surf(akx_mesh, aky_mesh, eig_enes_K(:,:,eig_index), squeeze(eig_layer_K(:,:,eig_index,:)), 'FaceAlpha',1, 'EdgeColor', 'interp');
shading interp % 去掉曲面上的网格
camlight % 添加光照，使得曲面更立体
lighting phong % 哑光效果，去除表面“纹路”
zlim([-30, 30]);
view(3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 手动设置渲染器
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exportsetupdlg; % 设置 Render 和 dpi 即可
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 背景透明化       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pause(6)
% set(gcf, 'color', 'none') % 背景透明化 
% set(gca, 'color', 'none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 背景透明化       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hole side
eig_index = 3;
[hole_peak_value_list, hole_peak_location_list] = find_lcm_matrix(squeeze(eig_enes_K(:,:,eig_index)));
hole_Nx = 4;
hole_Ny = 4;
hole_ene_delta = 1; % 1meV的能量间隔
hole_ene_start = -1.5;
hole_ene_level_list = zeros(1, hole_Nx * hole_Ny);
for i = 1:hole_Nx * hole_Ny
    hole_ene_level_list(i) = hole_ene_start - (i - 1) * hole_ene_delta;
end
figure
set(gcf, 'unit', 'inch', 'position', [10, 5, 10.00, 8.00]) % figure
for i = 1:hole_Nx * hole_Ny
    subplot(hole_Nx, hole_Ny, i)
    % contour(akx_list, aky_list, eig_enes_K(:,:,i), "ShowText",true, 'LevelList', ene_level_list(i))
    contour(akx_list, aky_list, eig_enes_K(:,:,eig_index), 'LevelList', hole_ene_level_list(i))
end

% electron side
eig_index = 4;
[ele_peak_value_list, ele_peak_location_list] = find_lcm_matrix(squeeze(-eig_enes_K(:,:,eig_index)));
ele_peak_value_list = - ele_peak_value_list;
ele_Nx = 4;
ele_Ny = 4;
ele_ene_delta = 1; % 1meV的能量间隔
ele_ene_start = 5.5;
ele_ene_level_list = zeros(1, hole_Nx * hole_Ny);
for i = 1:ele_Nx * ele_Ny
    ele_ene_level_list(i) = ele_ene_start + (i - 1) * ele_ene_delta;
end
figure
set(gcf, 'unit', 'inch', 'position', [10, 5, 10.00, 8.00]) % figure
for i = 1:ele_Nx * ele_Ny
    subplot(ele_Nx, ele_Ny, i)
    % contour(akx_list, aky_list, eig_enes_K(:,:,i), "ShowText",true, 'LevelList', ene_level_list(i))
    contour(akx_list, aky_list, eig_enes_K(:,:,eig_index), 'LevelList', ele_ene_level_list(i))
end