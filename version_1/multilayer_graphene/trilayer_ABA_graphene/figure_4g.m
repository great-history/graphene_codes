%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 运行前请先运行qhe_LL_plot.m 以获得跃迁参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\code_test")

gamma0 = hopping_params_before(1) * 1000; % 转化为meV单位
gamma1 = hopping_params_before(2) * 1000;
gamma2 = hopping_params_before(3) * 1000;
gamma3 = hopping_params_before(4) * 1000;  
gamma4 = hopping_params_before(5) * 1000;
gamma5 = hopping_params_before(6) * 1000;
delta = hopping_params_before(7) * 1000;
Delta2 = hopping_params_before(8) * 1000;
Delta1 = 0; % 0 / 25 / 50 / 150 / 180 / 250

% k_start / k_end 取为 [-0.1,-0.1], [0.1,0.1] 和 取为 [-1,-1], [1,1]
tic
Nx = 401;
Ny = 401;
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

%% MATLAB模式化作图
figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　画布(figure)和坐标轴(axes)的尺寸设置 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf, 'unit', 'inch', 'position', [10,10, 2.5, 3]) % figure
set(gca, 'Position', [0.18, 0.15, 0.80, 0.75]) % axess
set(gca, 'LabelFontSizeMultiplier', 1) % 这两句话要加上
set(gca, 'TitleFontSizeMultiplier', 1) % 这两句话要加上
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　坐标轴的方向                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca, 'XDir', 'normal')
set(gca, 'YDir', 'normal')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　坐标轴标签的命                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('ak', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial')
ylabel('E(meV)', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial')
set(gca, 'XColor', 'k')
set(gca, 'YColor', 'k')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　坐标轴的范围                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis([-0.04 0.04 -10 20])
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　坐标轴上的标度                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'xtick', [-0.04 -0.02 0 0.02 0.04]);
set(gca,'xticklabel', [-0.04 -0.02 0 0.02 0.04]);
xtickformat('%,.2f');
set(gca,'ytick', [-10, 0, 10, 20]);
set(gca,'yticklabel', [-10, 0, 10, 20]);
set(gca, 'TickLength', [0.01, 0.03])  % 第一个是Minor Tick，第二个是Major Tick，类似于Origin
set(gca, 'Box', 'on') % 框的四边
border_width = 1.5;
set(gca, 'LineWidth', border_width)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　绘制曲线                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
line_width = 0.5;
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
                plot(akx_list, eig_enes_K(center_y, :, i),'Color', red_RGB, 'LineWidth', line_width)
                hold on
            end

            % 双层
            for i = 3:6
                plot(akx_list, eig_enes_K(center_y, :, i),'Color', blue_RGB, 'LineWidth', line_width)
                hold on
            end
        else % Delta1不等于0
            % 单层
            for i = 1:6
                for j = 1:(Nx - 1)
                   plot(akx_list(j:j+1), eig_enes_K(center_y, j:j+1, i),'Color', 1 / 2 * (eig_clc_K(center_y, j, i, :) + eig_clc_K(center_y, j, i, :)), 'LineWidth', line_width)
                   hold on;
                end
            end
        end
    case ["y", "K"]
        if Delta1 == 0
             % 单层
            for i = 1:2
                plot(akx_list, eig_enes_K(:, center_x, i),'Color', red_RGB, 'LineWidth', line_width)
                hold on
            end

            % 双层
            for i = 3:6
                plot(akx_list, eig_enes_K(:, center_x, i),'Color', blue_RGB, 'LineWidth', line_width)
                hold on
            end
        else % Delta1不等于0
            % 单层
            for i = 1:6
                for j = 1:(Nx - 1)
                   plot(aky_list(j:j+1), eig_enes_K(j:j+1, center_x, i),'Color', 1 / 2 * (eig_clc_K(j, center_x, i, :) + eig_clc_K(j, center_x, i, :)), 'LineWidth', line_width)
                   hold on;
                end
            end
        end
    case ["x", "Kp"]
        if Delta1 == 0
             % 单层
            for i = 1:2
                plot(akx_list, eig_enes_Kp(center_y, :, i),'Color', red_RGB, 'LineWidth', line_width)
                hold on
            end

            % 双层
            for i = 3:6
                plot(akx_list, eig_enes_Kp(center_y, :, i),'Color', blue_RGB, 'LineWidth', line_width)
                hold on
            end
        else % Delta1不等于0
            % 单层
            for i = 1:6
                for j = 1:(Nx - 1)
                   plot(akx_list(j:j+1), eig_enes_Kp(center_y, j:j+1, i),'Color', 1 / 2 * (eig_clc_Kp(center_y, j, i, :) + eig_clc_Kp(center_y, j, i, :)), 'LineWidth', line_width)
                   hold on;
                end
            end
        end
    case ["y", "Kp"]
        if Delta1 == 0
             % 单层
            for i = 1:2
                plot(akx_list, eig_enes_Kp(:, center_x, i),'Color', red_RGB, 'LineWidth', line_width)
                hold on
            end

            % 双层
            for i = 3:6
                plot(akx_list, eig_enes_Kp(:, center_y, i),'Color', blue_RGB, 'LineWidth', line_width)
                hold on
            end
        else % Delta1不等于0
            % 单层
            for i = 1:6
                for j = 1:(Nx - 1)
                   plot(aky_list(j:j+1), eig_enes_Kp(j:j+1, center_y, i),'Color', 1 / 2 * (eig_clc_Kp(j, center_y, i, :) + eig_clc_Kp(j, center_y, i, :)), 'LineWidth', line_width)
                   hold on;
                end
            end
        end
    otherwise
        disp('No such type !!!')
end

%title('D = 0 V/nm', 'FontSize', 9, 'FontWeight', 'bold')
% title('At 1.0 Gpa', 'FontSize', 9, 'FontWeight', 'Normal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　设置图例(Legend)                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lg_current = legend([hm(1), hb(1)], {'Mono', 'Bi'}, 'Box', 'off', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7, 'Location', 'east');
% lg_current.ItemTokenSize = [10, 10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　将坐标轴的框置于顶层                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca, 'Layer', 'top')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　手动设置渲染器                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exportsetupdlg; % 设置 Render 和 dpi 即可
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%　背景透明化                          %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pause(6)
% set(gcf, 'color', 'none') % 背景透明化 
% set(gca, 'color', 'none')