classdef trilayer_ABA_class
    % 为常量属性赋值(一些基本的物理常数)
    properties (Constant)
        h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
        one_electron_coulomb = 1.602176634 * 10^(-19); % 单位是C
        epsilon_0 = 8.85 * 10^(-12) * 10^(-9); % 单位是F/nm
        % epsilon_bn = 6.6;
        d_interlayer = 0.335; % 单位是nm
        a_intralayer = 0.246; % 单位是nm
    end
    
    % 定义一些可变的参数，比如 hopping_params
    properties
        % hopping parameters
        gamma0
        gamma1
        gamma2  
        gamma3
        gamma4  
        gamma5  
        delta  
        Delta2
        
        % Delta1相关参数
        Delta1 % 当前的Delta1大小
        Delta1_start
        Delta1_end
        Delta1_steps
        Delta1_list
        eigvals_LL_K_ED
        eigvals_LL_Kp_ED
        eig_info_HK_ED_select_cell
        eig_info_HKp_ED_select_cell
        dos_ED_mat
        density_ED_mat
        dos_nD_mat
        flag_ED
        
        % B_field相关参数
        B_field % 当前的磁场强度
        B_start
        B_end
        B_steps
        B_fields_list
        eigvals_LL_K_EB
        eigvals_LL_Kp_EB
        eig_info_HK_EB_select_cell
        eig_info_HKp_EB_select_cell
        % eigvals_LL_K_b
        % eigvals_LL_Kp_b
        % eigvals_LL_K_m
        % eigvals_LL_Kp_m
        eig_info_HK_b_EB_select_cell
        eig_info_HKp_b_EB_select_cell
        eig_info_HK_m_EB_select_cell
        eig_info_HKp_m_EB_select_cell
        dos_EB_mat
        density_EB_mat
        dos_nB_mat
        flag_EB
        
        % energy window
        ene_ub
        ene_lb
        energy_cnp_list % 电中性点处的能量
        flag_cnp
        eigval_LL_K_b0_list % 存放bilayer-like branch的LL_K_0
        eigval_LL_Kp_b0_list % 存放bilayer-like branch的LL_Kp_0
        
        % 
        dist_index_array  % 分别用1，2，3代表intralayer, nearest interlayer, next-nearest interlayer
        T_mat             % 得到 sublattice basis 与 (m, b) basis 之间的变换矩阵
        
        % LL parameters
        LL_index_cutoff
        dims_m
        dims_b
        dims
        
        % 作图
        energy_list % 能量范围
        density_list % 载流子浓度范围
    end
    
    % 定义一些函数，比如 能带计算 // 朗道能级计算 // 朗道扇形图 // nD图
    methods (Access = 'public')
        % 初始化 model = trilayer_ABA_class(output_value_array(1,:));
        function model = trilayer_ABA_class(hopping_params)
            % 添加路径
            addpath('D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\utilities\');
            addpath('D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\trilayer_ABA_hamiltonian\');
            addpath('D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\trilayer_ABA_fitting_procedure\');
            addpath('D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\plot_funcs\');
            addpath('D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\interaction_effects\');
            
            model.flag_EB = false;
            model.flag_ED = false;
            model.flag_cnp = false;
            
            % 计算dist_array
            dist_index_array = zeros(6);
            for alpha = 1:6
                layer_alpha = mod(alpha, 2) + fix(alpha / 2);
                for beta = 1:6
                    layer_beta = mod(beta, 2) + fix(beta / 2);
                    if abs(layer_beta - layer_alpha) == 1
                        dist_index_array(alpha, beta) = 2;
                    elseif abs(layer_beta - layer_alpha) == 2
                        dist_index_array(alpha, beta) = 3;
                    else % 对应 alpha == beta
                        dist_index_array(alpha, beta) = 1;
                    end
                end
            end
            model.dist_index_array = dist_index_array;
            
            % 得到 sublattice basis 与 (m, b) basis 之间的变换矩阵
            T_mat = zeros(6); % 从{|alpha>}变换到{|phi_k>}的矩阵
            T_mat(1,1) = 1 / sqrt(2);
            T_mat(1,3) = 1 / sqrt(2);
            T_mat(2,2) = 1 / sqrt(2);
            T_mat(2,4) = 1 / sqrt(2);
            T_mat(3,5) = 1;
            T_mat(4,6) = 1;
            T_mat(5,1) = - 1 / sqrt(2);
            T_mat(5,3) = 1 / sqrt(2);
            T_mat(6,2) = - 1 / sqrt(2);
            T_mat(6,4) = 1 / sqrt(2);
            model.T_mat = T_mat;
            
            % Delta1是外场的参数
            model.gamma0 = hopping_params(1);
            model.gamma1 = hopping_params(2);
            model.gamma2 = hopping_params(3);
            model.gamma3 = hopping_params(4);  
            model.gamma4 = hopping_params(5);
            model.gamma5 = hopping_params(6);

            model.delta = hopping_params(7);
            model.Delta2 = hopping_params(8);
            % model.Delta1 = 0.0; 
        end
        
        % 朗道能级计算 // 朗道扇形图 // nD图 // nB图
        function model = trilayer_ABA_LLs_EB(model, Delta1, B_start, B_end, B_steps, LL_index_cutoff, ene_ub, ene_lb, ene_eps)
            % 选取能量位于[ene_lb - ene_eps, ene_ub + ene_eps]区间内的所有朗道能级
            % 计算朗道能级作为磁场和能量的函数
            if ene_eps < 0
                ene_eps = - ene_eps;
            end
            model.Delta1 = Delta1;
            model.B_start = B_start;
            model.B_end = B_end;
            model.B_steps = B_steps;
            model.B_fields_list = linspace(B_start, B_end, B_steps);
            model.energy_cnp_list = zeros(1, B_steps);
            
            model.ene_ub = ene_ub;
            model.ene_lb = ene_lb;
            
            % 计算矩阵的维数
            model.LL_index_cutoff = LL_index_cutoff;
            model.dims_m = 2 * LL_index_cutoff + 1;
            model.dims_b = 4 * LL_index_cutoff;
            model.dims = model.dims_b + model.dims_m;
            
            % 将gamma_i转换为v_i，同时某些gamma_i需要扩大sqrt(2)倍!!!!
            gamma1_new = model.gamma1 * sqrt(2);
            v0 = sqrt(3) / 2 * model.gamma0 * (model.a_intralayer * 10^(-9)) / model.h_bar;
            v3 = sqrt(3) / 2 * sqrt(2) * model.gamma3 * (model.a_intralayer * 10^(-9)) / model.h_bar; % 不要忘了乘以sqrt(2) !!!
            v4 = sqrt(3) / 2 * sqrt(2) * model.gamma4 * (model.a_intralayer * 10^(-9)) / model.h_bar; % 不要忘了乘以sqrt(2) !!!
            
            % 计算LL fan diagram
            tic
            if Delta1 == 0.0
                [model.eigvals_LL_K_EB, model.eigvals_LL_Kp_EB, model.eig_info_HK_b_EB_select_cell, model.eig_info_HKp_b_EB_select_cell, ...
                    model.eig_info_HK_m_EB_select_cell, model.eig_info_HKp_m_EB_select_cell] = ...
                            trilayer_ABA_LL_solver_EB_without_Delta1_with_ene_win(v0, v3, v4, gamma1_new, model.delta, model.gamma2, model.gamma5, model.Delta2, ...
                                                                                  model.B_fields_list, model.B_steps, model.LL_index_cutoff, model.dims_m, model.dims_b, model.ene_ub + ene_eps, model.ene_lb - ene_eps);
            else
                [model.eigvals_LL_K_EB, model.eigvals_LL_Kp_EB, model.eig_info_HK_EB_select_cell, model.eig_info_HKp_EB_select_cell] = ...
                            trilayer_ABA_LL_solver_EB_with_Delta1_with_ene_win(v0, v3, v4, gamma1_new, model.delta, model.gamma2, model.gamma5, model.Delta1, model.Delta2, ...
                                                                               model.B_fields_list, model.B_steps, model.LL_index_cutoff, model.dims_m, model.dims_b, model.ene_ub + ene_eps, ene_lb - ene_eps);
            end
            toc
            
            model.flag_ED = false;
            model.flag_EB = true;
            model.flag_cnp = false;
        end
        
        function model = trilayer_ABA_LLs_ED(model, B_field, Delta1_start, Delta1_end, Delta1_steps, LL_index_cutoff, ene_ub, ene_lb, ene_eps)
            % 计算朗道能级作为磁场和能量的函数
            if ene_eps < 0
                ene_eps = - ene_eps;
            end
            
            model.B_field = B_field;
            
            model.Delta1_steps = Delta1_steps;
            model.Delta1_start = Delta1_start;
            model.Delta1_end = Delta1_end;
            model.Delta1_list = linspace(Delta1_start, Delta1_end, Delta1_steps); % 以eV为单位
            model.energy_cnp_list = zeros(1, model.Delta1_steps);
            
            model.ene_ub = ene_ub;
            model.ene_lb = ene_lb;
            
            % 计算矩阵的维数
            model.LL_index_cutoff = LL_index_cutoff;
            model.dims_m = 2 * LL_index_cutoff + 1;
            model.dims_b = 4 * LL_index_cutoff;
            model.dims = model.dims_b + model.dims_m;
            
            % 将gamma_i转换为v_i，同时某些gamma_i需要扩大sqrt(2)倍!!!!
            gamma1_new = model.gamma1 * sqrt(2);
            v0 = sqrt(3) / 2 * model.gamma0 * (model.a_intralayer * 10^(-9)) / model.h_bar;
            v3 = sqrt(3) / 2 * sqrt(2) * model.gamma3 * (model.a_intralayer * 10^(-9)) / model.h_bar; % 不要忘了乘以sqrt(2) !!!
            v4 = sqrt(3) / 2 * sqrt(2) * model.gamma4 * (model.a_intralayer * 10^(-9)) / model.h_bar; % 不要忘了乘以sqrt(2) !!!
            
            % mag_length = 25.66 / sqrt(B_field); % 以nm为单位
            % x0 = model.h_bar * v0 / mag_length * 10^9; % 以eV为单位
            % x3 = model.h_bar * v3 / mag_length * 10^9; % 以eV为单位
            % x4 = model.h_bar * v4 / mag_length * 10^9; % 以eV为单位
            
            % 计算LL fan diagram as E & D
            tic
            [model.eigvals_LL_K_ED, model.eigvals_LL_Kp_ED, model.eig_info_HK_ED_select_cell, model.eig_info_HKp_ED_select_cell] = ...
                        trilayer_ABA_LL_solver_ED_fixed_B_with_ene_win(v0, v3, v4, gamma1_new, model.delta, model.gamma2, model.gamma5, model.Delta2, ...
                                                                       model.Delta1_list, Delta1_steps, B_field, LL_index_cutoff, model.dims_m, model.dims_b, ene_ub + ene_eps, ene_lb - ene_eps);
            toc
            model.flag_ED = true;
            model.flag_EB = false;
            model.flag_cnp = false;
        end
        
        function varargout = trilayer_ABA_LLs_EB_cnp_shift(model, cnp_index1, cnp_index2)
            % 退出机制
            if (~model.flag_EB)
                disp("还没进行LL EB的计算");
                return % 此时还没找出cnp的指标
            end
            
            % cnp_index 是电中性点上下的两个指标
            % get energy_cnp_list
            model.energy_cnp_list = zeros(1, model.B_steps);
            for B_index = model.B_steps:-1:1
                eigvals_LL_sort = sort([model.eigvals_LL_K_EB(B_index, :), model.eigvals_LL_Kp_EB(B_index, :)], 'ascend');
                model.energy_cnp_list(B_index) = (eigvals_LL_sort(cnp_index1) + eigvals_LL_sort(cnp_index2)) / 2;
            end

            % energy shifting with regard to cnp
            if model.Delta1 == 0.0
                for B_index = 1:model.B_steps
                    model.eig_info_HK_m_EB_select_cell{B_index, 5} = model.eig_info_HK_m_EB_select_cell{B_index, 2} - model.energy_cnp_list(B_index);
                    model.eig_info_HKp_m_EB_select_cell{B_index, 5} = model.eig_info_HKp_m_EB_select_cell{B_index, 2} - model.energy_cnp_list(B_index);
                    model.eig_info_HK_b_EB_select_cell{B_index, 5} = model.eig_info_HK_b_EB_select_cell{B_index, 2} - model.energy_cnp_list(B_index);
                    model.eig_info_HKp_b_EB_select_cell{B_index, 5} = model.eig_info_HKp_b_EB_select_cell{B_index, 2} - model.energy_cnp_list(B_index);
                end
            else
                for B_index = 1:model.B_steps
                    model.eig_info_HK_EB_select_cell{B_index, 5} = model.eig_info_HK_EB_select_cell{B_index, 2} - model.energy_cnp_list(B_index);
                    model.eig_info_HKp_EB_select_cell{B_index, 5} = model.eig_info_HKp_EB_select_cell{B_index, 2} - model.energy_cnp_list(B_index);
                end
            end
            
            model.flag_cnp = true;
            varargout{1} = model;
        end
        
        function varargout = trilayer_ABA_LLs_ED_cnp_shift(model, cnp_index1, cnp_index2)
            % 退出机制
            if (~model.flag_ED)
                disp("还没进行LL ED的计算");
                return % 此时还没找出cnp的指标
            end
            
            % get energy_cnp_list
            model.energy_cnp_list = zeros(1, model.Delta1_steps);
            for D_index = model.Delta1_steps:-1:1
                eigvals_LL_sort = sort([model.eigvals_LL_K_ED(D_index, :), model.eigvals_LL_Kp_ED(D_index, :)], 'ascend');
                model.energy_cnp_list(D_index) = (eigvals_LL_sort(cnp_index1) + eigvals_LL_sort(cnp_index2)) / 2;
            end

            % energy shifting with regard to cnp
            for D_index = 1:model.Delta1_steps
                model.eig_info_HK_ED_select_cell{D_index, 5} = model.eig_info_HK_ED_select_cell{D_index, 2} - model.energy_cnp_list(D_index);
                model.eig_info_HKp_ED_select_cell{D_index, 5} = model.eig_info_HKp_ED_select_cell{D_index, 2} - model.energy_cnp_list(D_index);
            end
            
            model.flag_cnp = true;
            varargout{1} = model;
        end
        
        function varargout = trilayer_ABA_LLs_nB(model, ene_broadening, steps, density_steps) % varargout可以是[model, ... ...]或者nothing
            % 计算态密度(DOS)作为n和B的函数  % ene_brodening LL的能量展宽  % steps作为一个调控能量间隔的参数
            if ~model.flag_EB
                disp("还没进行LL EB的计算");
                return
            end
            
            Ene_steps = round((model.ene_ub - model.ene_lb) / ene_broadening) * 2 * steps + 1;
            energy_list_new = linspace(model.ene_lb, model.ene_ub, Ene_steps);
            delta_ene = abs(energy_list_new(1) - energy_list_new(2));
            
            if model.Delta1 == 0.0
                dos_HK_m_EB_mat  = get_dos_asfo_EB(model.eig_info_HK_m_EB_select_cell,  energy_list_new, Ene_steps, model.B_fields_list, model.B_steps, ene_broadening);
                dos_HKp_m_EB_mat = get_dos_asfo_EB(model.eig_info_HKp_m_EB_select_cell, energy_list_new, Ene_steps, model.B_fields_list, model.B_steps, ene_broadening);
                dos_HK_b_EB_mat  = get_dos_asfo_EB(model.eig_info_HK_b_EB_select_cell,  energy_list_new, Ene_steps, model.B_fields_list, model.B_steps, ene_broadening);
                dos_HKp_b_EB_mat = get_dos_asfo_EB(model.eig_info_HKp_b_EB_select_cell, energy_list_new, Ene_steps, model.B_fields_list, model.B_steps, ene_broadening);

                model.dos_EB_mat = dos_HK_m_EB_mat + dos_HKp_m_EB_mat + dos_HK_b_EB_mat + dos_HKp_b_EB_mat;
                
                [model.density_EB_mat, density_max, density_min] = get_n_asfo_EB(model.dos_EB_mat, energy_list_new, Ene_steps, model.B_steps, delta_ene);
                model.density_list = linspace(density_min, density_max, density_steps);
                model.dos_nB_mat = get_dos_asfo_nB(model.dos_EB_mat, model.density_EB_mat, model.density_list, model.B_steps, density_steps);
                model.energy_list = energy_list_new;
                
                varargout{1} = model;
                varargout{2} = dos_HK_m_EB_mat;
                varargout{3} = dos_HKp_m_EB_mat;
                varargout{4} = dos_HK_b_EB_mat;
                varargout{5} = dos_HKp_b_EB_mat;
            else
                dos_HK_EB_mat  = get_dos_asfo_EB(model.eig_info_HK_EB_select_cell,  energy_list_new, Ene_steps, model.B_fields_list, model.B_steps, ene_broadening);
                dos_HKp_EB_mat = get_dos_asfo_EB(model.eig_info_HKp_EB_select_cell, energy_list_new, Ene_steps, model.B_fields_list, model.B_steps, ene_broadening);

                model.dos_EB_mat = dos_HK_EB_mat + dos_HKp_EB_mat;
                
                [model.density_EB_mat, density_max, density_min] = get_n_asfo_EB(model.dos_EB_mat, energy_list_new, Ene_steps, model.B_steps, delta_ene);
                model.density_list = linspace(density_min, density_max, density_steps);
                model.dos_nB_mat = get_dos_asfo_nB(model.dos_EB_mat, model.density_EB_mat, model.density_list, model.B_steps, density_steps);
                model.energy_list = energy_list_new;
                
                varargout{1} = model;
                varargout{2} = dos_HK_EB_mat;
                varargout{3} = dos_HKp_EB_mat;
            end

        end
        
        function varargout = trilayer_ABA_LLs_nD(model, ene_broadening, steps, density_steps) % varargout可以是[model, ... ...]或者nothing
            % 计算态密度(DOS)作为n和B的函数  % ene_brodening LL的能量展宽  % steps作为一个调控能量间隔的参数
            if ~model.flag_ED
                disp("还没进行LL ED的计算");
                return
            end
            
            Ene_steps = round((model.ene_ub - model.ene_lb) / ene_broadening) * 2 * steps + 1;
            energy_list_new = linspace(model.ene_lb, model.ene_ub, Ene_steps);
            delta_ene = abs(energy_list_new(1) - energy_list_new(2));
            
            % 计算态密度(DOS) as a function of E & D
            dos_HK_ED_mat  = get_dos_asfo_ED(model.eig_info_HK_ED_select_cell,  energy_list_new, Ene_steps, model.Delta1_list, model.Delta1_steps, model.B_field, ene_broadening);
            dos_HKp_ED_mat = get_dos_asfo_ED(model.eig_info_HKp_ED_select_cell, energy_list_new, Ene_steps, model.Delta1_list, model.Delta1_steps, model.B_field, ene_broadening);
            model.dos_ED_mat = dos_HK_ED_mat + dos_HKp_ED_mat;

            % 计算载流子浓度(n) as a function of E & D
            [model.density_ED_mat, density_max, density_min] = get_n_asfo_ED(model.dos_ED_mat, energy_list_new, Ene_steps, model.Delta1_steps, delta_ene);

            % 计算态密度(DOS) as a function of n & D
            model.density_list = linspace(density_min, density_max, density_steps);
            model.dos_nD_mat = get_dos_asfo_nD(model.dos_ED_mat, model.density_ED_mat, model.density_list, model.Delta1_steps, density_steps);
            model.energy_list = energy_list_new;
            
            % 输出
            varargout{1} = model;
            varargout{2} = dos_HK_ED_mat;
            varargout{3} = dos_HKp_ED_mat;
            
        end
        
        % 作图
        function current_fig = representative_imgs_plot(model)
            current_fig = figure;
            if model.flag_EB
                subplot(1,3,1)
                imagesc(model.B_fields_list, model.energy_list, model.dos_EB_mat');
                set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
                colorbar
                subplot(1,3,2)
                imagesc(model.B_fields_list, model.energy_list, model.density_EB_mat');
                set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
                colorbar
                subplot(1,3,3)
                % imagesc(B_fields_list, density_list, dos_nB_mat');
                imagesc(model.density_list, model.B_fields_list, model.dos_nB_mat);
                set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
                colorbar
            elseif model.flag_ED
                subplot(1,3,1)
                imagesc(model.Delta1_list, model.energy_list, model.dos_ED_mat');
                set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
                colorbar
                subplot(1,3,2)
                imagesc(model.Delta1_list, model.energy_list, model.density_ED_mat');
                set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
                colorbar
                subplot(1,3,3)
                imagesc(model.density_list, model.Delta1_list, model.dos_nD_mat);
                set(gca, 'YDir', 'normal') % imagsec默认y轴是上小下大
                colorbar
            end
        end
        
        function current_fig = representative_lines_plot(model, line_info_cell, ene_lb_new, ene_ub_new, save_path)
            % line_info_cell 包含了有关曲线的 线型 / 颜色 / 粗细 / ......
            if model.flag_EB
                if model.Delta1 == 0.0
                    eigvals_LL_cell = cell(4, 4); % 分为 m_K / m_Kp / b_K / b_Kp 四种
                    eigvals_LL_cell{1,1} = model.eigvals_LL_K_EB(:, 1:model.dims_m); % m_K
                    eigvals_LL_cell{1,2} = model.dims_m; % 维数
                    eigvals_LL_cell{1,3} = line_info_cell{1, 1}; % line color 'b--'
                    eigvals_LL_cell{1,4} = line_info_cell{1, 2}; % line width

                    eigvals_LL_cell{2,1} = model.eigvals_LL_Kp_EB(:, 1:model.dims_m); % m_Kp
                    eigvals_LL_cell{2,2} = model.dims_m;
                    eigvals_LL_cell{2,3} = line_info_cell{2, 1};
                    eigvals_LL_cell{2,4} = line_info_cell{2, 2};
                    
                    eigvals_LL_cell{3,1} = model.eigvals_LL_K_EB(:, (1 + model.dims_m):end);
                    eigvals_LL_cell{3,2} = model.dims_b; % 维数
                    eigvals_LL_cell{3,3} = line_info_cell{3, 1}; % line color 'b--'
                    eigvals_LL_cell{3,4} = line_info_cell{3, 2}; % line width

                    eigvals_LL_cell{4,1} = model.eigvals_LL_Kp_EB(:, (1 + model.dims_m):end);
                    eigvals_LL_cell{4,2} = model.dims_b;
                    eigvals_LL_cell{4,3} = line_info_cell{4, 1};
                    eigvals_LL_cell{4,4} = line_info_cell{4, 2};
                else
                    eigvals_LL_cell = cell(2, 4);
                    eigvals_LL_cell{1,1} = model.eigvals_LL_K_EB;
                    eigvals_LL_cell{1,2} = model.dims; % 维数
                    eigvals_LL_cell{1,3} = line_info_cell{1, 1}; % line color 'b--'
                    eigvals_LL_cell{1,4} = line_info_cell{1, 2}; % line width

                    eigvals_LL_cell{2,1} = model.eigvals_LL_Kp_EB;
                    eigvals_LL_cell{2,2} = model.dims;
                    eigvals_LL_cell{2,3} = line_info_cell{2, 1};
                    eigvals_LL_cell{2,4} = line_info_cell{2, 2}; 
                end
                
                current_fig = plot_LLs(model.B_start, model.B_end, ene_lb_new, ene_ub_new, model.B_fields_list, eigvals_LL_cell, save_path);
                
            elseif model.flag_ED
                eigvals_LL_cell_Delta1 = cell(2, 4);
                eigvals_LL_cell_Delta1{1,1} = model.eigvals_LL_K_ED;
                eigvals_LL_cell_Delta1{1,2} = model.dims; % 维数
                eigvals_LL_cell_Delta1{1,3} = line_info_cell{1, 1}; % line color
                eigvals_LL_cell_Delta1{1,4} = line_info_cell{1, 2}; % line width

                eigvals_LL_cell_Delta1{2,1} = model.eigvals_LL_Kp_ED;
                eigvals_LL_cell_Delta1{2,2} = model.dims;
                eigvals_LL_cell_Delta1{2,3} = line_info_cell{2, 1};
                eigvals_LL_cell_Delta1{2,4} = line_info_cell{2, 2};

                current_fig = plot_LLs(model.Delta1_start, model.Delta1_end, ene_lb_new, ene_ub_new, model.Delta1_list, eigvals_LL_cell_Delta1, save_path);
                
            end
        end
        
        function current_fig = select_LL_line_plot(current_fig)
            % type : 
            figure(current_fig);
            
            
        end
        
        function index = find_closet_field_index(model, field)
            if model.flag_EB
                diff_list = abs(model.B_fields_list - field);
                [~, index] = min(diff_list);
            elseif model.flag_ED
                diff_list = abs(model.Delta1_list - field);
                [~, index] = min(diff_list);
            else
                index= nan;
            end
        end
        
        % 提取出zeroth LL的信息
        function zeroth_LL_info_cell = get_zeroth_LL_info_cell(model, B_index1, B_index2)
            % 用来存放zeroth LL信息的cell
            zeroth_LL_info_cell = cell(6, 4); % 第一个指标存放名称(LL_K_m0 // LL_Kp_m0 // LL_K_b0 // LL_Kp_b0 // LL_K_b1 // LL_Kp_b1)
                                              % 第二个指标是 本征值(eigvals)
                                              % 第三个指标是 本征态各个component的分量(在基组phi_k下)
                                              % 第四个指标是 本征态各个component的分量(在基组alpha下)
            
            B_steps_new = B_index2 - B_index1 + 1;
            % LL_K_m0
            LL_K_m0_eigvals_list = zeros(1, B_steps_new);
            LL_K_m0_eigvec_phi_component_array  = zeros(model.LL_index_cutoff + 1, 2, B_steps_new);
            LL_K_m0_eigvec_alpha_component_array  = zeros(model.LL_index_cutoff + 1, 6, B_steps_new);
            % LL_Kp_m0
            LL_Kp_m0_eigvals_list = zeros(1, B_steps_new);
            LL_Kp_m0_eigvec_phi_component_array = zeros(model.LL_index_cutoff + 1, 2, B_steps_new);
            LL_Kp_m0_eigvec_alpha_component_array = zeros(model.LL_index_cutoff + 1, 6, B_steps_new);
            % LL_K_b0
            LL_K_b0_eigvals_list  =  zeros(1, B_steps_new);
            LL_K_b0_eigvec_phi_component_array   = zeros(model.LL_index_cutoff + 1, 4, B_steps_new);
            LL_K_b0_eigvec_alpha_component_array   = zeros(model.LL_index_cutoff + 1, 6, B_steps_new);
            % LL_Kp_b0
            LL_Kp_b0_eigvals_list =  zeros(1, B_steps_new);
            LL_Kp_b0_eigvec_phi_component_array = zeros(model.LL_index_cutoff + 1, 4, B_steps_new);
            LL_Kp_b0_eigvec_alpha_component_array = zeros(model.LL_index_cutoff + 1, 6, B_steps_new);
            % LL_K_b1
            LL_K_b1_eigvals_list  =  zeros(1, B_steps_new);
            LL_K_b1_eigvec_phi_component_array = zeros(model.LL_index_cutoff + 1, 4, B_steps_new);
            LL_K_b1_eigvec_alpha_component_array = zeros(model.LL_index_cutoff + 1, 6, B_steps_new);
            % LL_Kp_b1
            LL_Kp_b1_eigvals_list =  zeros(1, B_steps_new);
            LL_Kp_b1_eigvec_phi_component_array = zeros(model.LL_index_cutoff + 1, 4, B_steps_new);
            LL_Kp_b1_eigvec_alpha_component_array = zeros(model.LL_index_cutoff + 1, 6, B_steps_new);
            
            count = 0;
            for B_index = B_index1:B_index2
                count = count + 1;
                LL_K_m0_index = model.eig_info_HK_m_EB_select_cell{B_index, 4} == 0;
                LL_K_m0_eigvals_list(count) = model.eig_info_HK_m_EB_select_cell{B_index, 2}(LL_K_m0_index);
                eigvec_K_m0 = model.eig_info_HK_m_EB_select_cell{B_index, 1}(:, LL_K_m0_index);
                LL_K_m0_eigvec_phi_component_array(:, :, count) = get_LL_m_components_each_phi_k(eigvec_K_m0, model.LL_index_cutoff, + 1);
                LL_K_m0_eigvec_alpha_component_array(:, :, count) = get_LL_m_components_each_alpha(model.T_mat, LL_K_m0_eigvec_phi_component_array(:, :, count), model.LL_index_cutoff);
                
                LL_Kp_m0_index = model.eig_info_HKp_m_EB_select_cell{B_index, 4} == 0;
                LL_Kp_m0_eigvals_list(count) = model.eig_info_HKp_m_EB_select_cell{B_index, 2}(LL_Kp_m0_index);
                eigvec_Kp_m0 = model.eig_info_HKp_m_EB_select_cell{B_index, 1}(:, LL_Kp_m0_index);
                LL_Kp_m0_eigvec_phi_component_array(:, :, count) = get_LL_m_components_each_phi_k(eigvec_Kp_m0, model.LL_index_cutoff, - 1);
                LL_Kp_m0_eigvec_alpha_component_array(:, :, count) = get_LL_m_components_each_alpha(model.T_mat, LL_Kp_m0_eigvec_phi_component_array(:, :, count), model.LL_index_cutoff);

                LL_K_b0_index = model.eig_info_HK_b_EB_select_cell{B_index, 4} == 0;
                LL_K_b0_eigvals_list(count) = model.eig_info_HK_b_EB_select_cell{B_index, 2}(LL_K_b0_index);
                eigvec_K_b0 = model.eig_info_HK_b_EB_select_cell{B_index, 1}(:, LL_K_b0_index);
                LL_K_b0_eigvec_phi_component_array (:, :, count) = get_LL_b_components_each_phi_k(eigvec_K_b0, model.LL_index_cutoff, + 1);
                LL_K_b0_eigvec_alpha_component_array(:, :, count) = get_LL_b_components_each_alpha(model.T_mat, LL_K_b0_eigvec_phi_component_array(:, :, count), model.LL_index_cutoff);
                
                LL_Kp_b0_index = model.eig_info_HKp_b_EB_select_cell{B_index, 4} == 0;
                LL_Kp_b0_eigvals_list(count) = model.eig_info_HKp_b_EB_select_cell{B_index, 2}(LL_Kp_b0_index);
                eigvec_Kp_b0 = model.eig_info_HKp_b_EB_select_cell{B_index, 1}(:, LL_Kp_b0_index);
                LL_Kp_b0_eigvec_phi_component_array(:, :, count) = get_LL_b_components_each_phi_k(eigvec_Kp_b0, model.LL_index_cutoff, - 1);
                LL_Kp_b0_eigvec_alpha_component_array(:, :, count) = get_LL_b_components_each_alpha(model.T_mat, LL_Kp_b0_eigvec_phi_component_array(:, :, count), model.LL_index_cutoff);
                
                LL_K_b1_index = model.eig_info_HK_b_EB_select_cell{B_index, 4} == 1;
                LL_K_b1_eigvals_list(count) = model.eig_info_HK_b_EB_select_cell{B_index, 2}(LL_K_b1_index);
                eigvec_K_b1 = model.eig_info_HK_b_EB_select_cell{B_index, 1}(:, LL_K_b1_index);
                LL_K_b1_eigvec_phi_component_array(:, :, count) = get_LL_b_components_each_phi_k(eigvec_K_b1, model.LL_index_cutoff, + 1);
                LL_K_b1_eigvec_alpha_component_array(:, :, count) = get_LL_b_components_each_alpha(model.T_mat, LL_K_b1_eigvec_phi_component_array(:, :, count), model.LL_index_cutoff);
                
                LL_Kp_b1_index = model.eig_info_HKp_b_EB_select_cell{B_index, 4} == 1;
                LL_Kp_b1_eigvals_list(count) = model.eig_info_HKp_b_EB_select_cell{B_index, 2}(LL_Kp_b1_index);
                eigvec_Kp_b1 = model.eig_info_HKp_b_EB_select_cell{B_index, 1}(:, LL_Kp_b1_index);
                LL_Kp_b1_eigvec_phi_component_array(:, :, count) = get_LL_b_components_each_phi_k(eigvec_Kp_b1, model.LL_index_cutoff, - 1);
                LL_Kp_b1_eigvec_alpha_component_array(:, :, count) = get_LL_b_components_each_alpha(model.T_mat, LL_Kp_b1_eigvec_phi_component_array(:, :, count), model.LL_index_cutoff);
            end
            
            zeroth_LL_info_cell{1,1} = "LL_K_m0";
            zeroth_LL_info_cell{1,2} = LL_K_m0_eigvals_list;
            zeroth_LL_info_cell{1,3} = LL_K_m0_eigvec_phi_component_array;
            zeroth_LL_info_cell{1,4} = LL_K_m0_eigvec_alpha_component_array;
            
            zeroth_LL_info_cell{2,1} = "LL_Kp_m0";
            zeroth_LL_info_cell{2,2} = LL_Kp_m0_eigvals_list;
            zeroth_LL_info_cell{2,3} = LL_Kp_m0_eigvec_phi_component_array;
            zeroth_LL_info_cell{2,4} = LL_Kp_m0_eigvec_alpha_component_array;
            
            zeroth_LL_info_cell{3,1} = "LL_K_b0";
            zeroth_LL_info_cell{3,2} = LL_K_b0_eigvals_list;
            zeroth_LL_info_cell{3,3} = LL_K_b0_eigvec_phi_component_array ;
            zeroth_LL_info_cell{3,4} = LL_K_b0_eigvec_alpha_component_array;
            
            zeroth_LL_info_cell{4,1} = "LL_Kp_b0";
            zeroth_LL_info_cell{4,2} = LL_Kp_b0_eigvals_list;
            zeroth_LL_info_cell{4,3} = LL_Kp_b0_eigvec_phi_component_array;
            zeroth_LL_info_cell{4,4} = LL_Kp_b0_eigvec_alpha_component_array;
            
            zeroth_LL_info_cell{5,1} = "LL_K_b1";
            zeroth_LL_info_cell{5,2} = LL_K_b1_eigvals_list;
            zeroth_LL_info_cell{5,3} = LL_K_b1_eigvec_phi_component_array;
            zeroth_LL_info_cell{5,4} = LL_K_b1_eigvec_alpha_component_array;
            
            zeroth_LL_info_cell{6,1} = "LL_Kp_b1";
            zeroth_LL_info_cell{6,2} = LL_Kp_b1_eigvals_list;
            zeroth_LL_info_cell{6,3} = LL_Kp_b1_eigvec_phi_component_array;
            zeroth_LL_info_cell{6,4} = LL_Kp_b1_eigvec_alpha_component_array;
            
        end
        
        % 提取出Lowest LL的信息
        function lowest_LL_info_cell = get_lowest_LL_info_cell(model, B_index1, B_index2, LL_K_m_index_bound, LL_Kp_m_index_bound, LL_K_b_index_bound, LL_Kp_b_index_bound)
            % 用来存放lowest LL信息的cell
            B_steps_new = B_index2 - B_index1 + 1;
            num_LL_K_m  = LL_K_m_index_bound(2) - LL_K_m_index_bound(1) + 1;
            num_LL_Kp_m = LL_Kp_m_index_bound(2) - LL_Kp_m_index_bound(1) + 1;
            num_LL_K_b  = LL_K_b_index_bound(2) - LL_K_b_index_bound(1) + 1;
            num_LL_Kp_b = LL_Kp_b_index_bound(2) - LL_Kp_b_index_bound(1) + 1;
            num_LL = num_LL_K_m + num_LL_Kp_m + num_LL_K_b + num_LL_Kp_b;
            lowest_LL_info_cell = cell(num_LL, 4); % 第一个指标存放名称(LL_K_m0 // LL_Kp_m0 // LL_K_b0 // LL_Kp_b0 // LL_K_b1 // LL_Kp_b1)
                                                   % 第二个指标是 本征值(eigvals)
                                                   % 第三个指标是 本征态各个component的分量(在基组phi_k下)
                                                   % 第四个指标是 本征态各个component的分量(在基组alpha下)
            % LL_K_m_x
            for ii = 1:num_LL_K_m
                LL_K_mx_index = LL_K_m_index_bound(1) + ii - 1;
                LL_K_mx_eigvals_list = zeros(1, B_steps_new);
                LL_K_mx_eigvec_phi_component_array  = zeros(model.LL_index_cutoff + 1, 2, B_steps_new);
                LL_K_mx_eigvec_alpha_component_array  = zeros(model.LL_index_cutoff + 1, 6, B_steps_new);
                
                count = 0;
                for B_index = B_index1:B_index2
                    count = count + 1;
                    index = model.eig_info_HK_m_EB_select_cell{B_index, 4} == LL_K_mx_index;
                    LL_K_mx_eigvals_list(count) = model.eig_info_HK_m_EB_select_cell{B_index, 2}(index);
                    eigvec_K_mx = model.eig_info_HK_m_EB_select_cell{B_index, 1}(:, index);
                    LL_K_mx_eigvec_phi_component_array(:, :, count) = get_LL_m_components_each_phi_k(eigvec_K_mx, model.LL_index_cutoff, + 1);
                    LL_K_mx_eigvec_alpha_component_array(:, :, count) = get_LL_m_components_each_alpha(model.T_mat, LL_K_mx_eigvec_phi_component_array(:, :, count), model.LL_index_cutoff);
                end
                
                lowest_LL_info_cell{ii,1} = ['LL_K_m', num2str(LL_K_mx_index)];
                lowest_LL_info_cell{ii,2} = LL_K_mx_eigvals_list;
                lowest_LL_info_cell{ii,3} = LL_K_mx_eigvec_phi_component_array;
                lowest_LL_info_cell{ii,4} = LL_K_mx_eigvec_alpha_component_array;
            end
            
            % LL_Kp_m_x
            for ii = 1:num_LL_Kp_m
                LL_Kp_mx_index = LL_Kp_m_index_bound(1) + ii - 1;
                LL_Kp_mx_eigvals_list = zeros(1, B_steps_new);
                LL_Kp_mx_eigvec_phi_component_array  = zeros(model.LL_index_cutoff + 1, 2, B_steps_new);
                LL_Kp_mx_eigvec_alpha_component_array  = zeros(model.LL_index_cutoff + 1, 6, B_steps_new);
                
                count = 0;
                for B_index = B_index1:B_index2
                    count = count + 1;
                    index = model.eig_info_HKp_m_EB_select_cell{B_index, 4} == LL_Kp_mx_index;
                    LL_Kp_mx_eigvals_list(count) = model.eig_info_HKp_m_EB_select_cell{B_index, 2}(index);
                    eigvec_Kp_mx = model.eig_info_HKp_m_EB_select_cell{B_index, 1}(:, index);
                    LL_Kp_mx_eigvec_phi_component_array(:, :, count) = get_LL_m_components_each_phi_k(eigvec_Kp_mx, model.LL_index_cutoff, - 1);
                    LL_Kp_mx_eigvec_alpha_component_array(:, :, count) = get_LL_m_components_each_alpha(model.T_mat, LL_Kp_mx_eigvec_phi_component_array(:, :, count), model.LL_index_cutoff);
                end
                
                lowest_LL_info_cell{ii + num_LL_K_m,1} = ['LL_Kp_m', num2str(LL_Kp_mx_index)];
                lowest_LL_info_cell{ii + num_LL_K_m,2} = LL_Kp_mx_eigvals_list;
                lowest_LL_info_cell{ii + num_LL_K_m,3} = LL_Kp_mx_eigvec_phi_component_array;
                lowest_LL_info_cell{ii + num_LL_K_m,4} = LL_Kp_mx_eigvec_alpha_component_array;
            end
            
            % LL_K_b_x
            for ii = 1:num_LL_K_b
                LL_K_bx_index = LL_K_b_index_bound(1) + ii - 1;
                LL_K_bx_eigvals_list = zeros(1, B_steps_new);
                LL_K_bx_eigvec_phi_component_array  = zeros(model.LL_index_cutoff + 1, 4, B_steps_new);
                LL_K_bx_eigvec_alpha_component_array  = zeros(model.LL_index_cutoff + 1, 6, B_steps_new);
                
                count = 0;
                for B_index = B_index1:B_index2
                    count = count + 1;
                    index = model.eig_info_HK_b_EB_select_cell{B_index, 4} == LL_K_bx_index;
                    LL_K_bx_eigvals_list(count) = model.eig_info_HK_b_EB_select_cell{B_index, 2}(index);
                    eigvec_K_bx = model.eig_info_HK_b_EB_select_cell{B_index, 1}(:, index);
                    LL_K_bx_eigvec_phi_component_array(:, :, count) = get_LL_b_components_each_phi_k(eigvec_K_bx, model.LL_index_cutoff, + 1);
                    LL_K_bx_eigvec_alpha_component_array(:, :, count) = get_LL_b_components_each_alpha(model.T_mat, LL_K_bx_eigvec_phi_component_array(:, :, count), model.LL_index_cutoff);
                end
                
                lowest_LL_info_cell{ii + num_LL_K_m + num_LL_Kp_m,1} = ['LL_K_b', num2str(LL_K_bx_index)];
                lowest_LL_info_cell{ii + num_LL_K_m + num_LL_Kp_m,2} = LL_K_bx_eigvals_list;
                lowest_LL_info_cell{ii + num_LL_K_m + num_LL_Kp_m,3} = LL_K_bx_eigvec_phi_component_array;
                lowest_LL_info_cell{ii + num_LL_K_m + num_LL_Kp_m,4} = LL_K_bx_eigvec_alpha_component_array;
            end
            
            % LL_Kp_b_x
            for ii = 1:num_LL_Kp_b
                LL_Kp_bx_index = LL_Kp_b_index_bound(1) + ii - 1;
                LL_Kp_bx_eigvals_list = zeros(1, B_steps_new);
                LL_Kp_bx_eigvec_phi_component_array  = zeros(model.LL_index_cutoff + 1, 4, B_steps_new);
                LL_Kp_bx_eigvec_alpha_component_array  = zeros(model.LL_index_cutoff + 1, 6, B_steps_new);
                
                count = 0;
                for B_index = B_index1:B_index2
                    count = count + 1;
                    index = model.eig_info_HKp_b_EB_select_cell{B_index, 4} == LL_Kp_bx_index;
                    LL_Kp_bx_eigvals_list(count) = model.eig_info_HKp_b_EB_select_cell{B_index, 2}(index);
                    eigvec_Kp_bx = model.eig_info_HKp_b_EB_select_cell{B_index, 1}(:, index);
                    LL_Kp_bx_eigvec_phi_component_array(:, :, count) = get_LL_b_components_each_phi_k(eigvec_Kp_bx, model.LL_index_cutoff, - 1);
                    LL_Kp_bx_eigvec_alpha_component_array(:, :, count) = get_LL_b_components_each_alpha(model.T_mat, LL_Kp_bx_eigvec_phi_component_array(:, :, count), model.LL_index_cutoff);
                end
                
                lowest_LL_info_cell{ii + num_LL_K_m + num_LL_Kp_m + num_LL_K_b,1} = ['LL_Kp_b', num2str(LL_Kp_bx_index)];
                lowest_LL_info_cell{ii + num_LL_K_m + num_LL_Kp_m + num_LL_K_b,2} = LL_Kp_bx_eigvals_list;
                lowest_LL_info_cell{ii + num_LL_K_m + num_LL_Kp_m + num_LL_K_b,3} = LL_Kp_bx_eigvec_phi_component_array;
                lowest_LL_info_cell{ii + num_LL_K_m + num_LL_Kp_m + num_LL_K_b,4} = LL_Kp_bx_eigvec_alpha_component_array;
            end
        end
        
        % 计算相互作用涉及的能量尺度(energy scale)
        function [E_zeeman, E_F, E_exchange, E_H] = get_energy_scales(model, epsilon_bn, B_field_now, d_interlayer_now)
            % epsilon_bn        % BN的介电常数(无量纲化)
            % E_zeeman          % 计算公式 ：0.116 * B_field / 1000; % 单位是eV
            % E_F               % 计算公式 ：one_electron_coulomb / (4 * pi * epsilon_0 * mag_length); % 单位是 eV fock interaction strength
            % E_exchange        % 计算公式 ：E_F / epsilon_bn; % 由于hBN存在介电常数（这里是in-plane dielectric constant），这里取为6.6，它可以有效减小交换相互作用的强度
            % E_H               % 计算公式 ：d_interlayer / (2 * mag_length) * E_F; % 单位是 eV hartree interaction strength, 注意在ABC trilayer那里是(2 * d / l_b)
            mag_length = 25.6 / sqrt(B_field_now); % 以nm为单位
            E_zeeman = 0.116 * B_field_now / 1000; % 单位是eV
            E_F = model.one_electron_coulomb / (4 * pi * model.epsilon_0 * mag_length); % 单位是 eV fock interaction strength
            E_exchange = E_F / epsilon_bn; % 由于hBN存在介电常数（这里是in-plane dielectric constant），这里取为6.6，它可以有效减小交换相互作用的强度
            E_H = d_interlayer_now / (2 * mag_length) * E_F; % 单位是 eV hartree interaction strength, 注意在ABC trilayer那里是(2 * d / l_b)
        end
    end
    
    % 静态方法
    methods (Static)
        
    end
    
    % 定义一些只能内部访问的函数
    methods (Access = 'private') % Access by class members only 
        
    end
    
end % classdef