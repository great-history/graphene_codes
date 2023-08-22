% todo:Kp的寻找low_ene
%% 一些参数
onsite_a = 0.015; % 0.015
onsite_b = 0.0;

vf = 2.1354; % eV为单位

%% parameters from PRB 2,033150(2020) Topological flat bands and correlated states in twisted monolayer-bilayer graphene
% gamma0 = 2.61;
% gamma1 = 0.361;
% gamma3 = 0.283;
% gamma4 = 0.140;
% delta_dimer = 0.050;
% Delta_asymm = 0.0; % 0.0 / 0.005 / 0.020 / -0.005 / -0.020
% flag_chiral = true;  % true: AB-AB / false: AB-BA
% 
% w_aa = 0.055;
% w_ab = 0.11;
% theta = 1.08 / 180 * pi; % 角度制要转化为弧度制 0.9 / 1.05 / 1.1 / 1.33 / 1.5 / 2.0 

%% parameters from Topological flat bands in twisted trilayer graphene
% gamma0_list = [2.464, 2.464, 2.464, 2.610, 2.610, 2.610];
% gamma1_list = [0.4, 0.4, 0.4, 0.360, 0.360, 0.360];
% gamma3_list = [0.2, 0.320, 0.320, 0.283, 0.283, 0.283];
% gamma4_list = [0.0, 0.138, 0.044, 0.138, 0.138, 0.138];
% delta_dimer_list = [0.0, 0.0, 0.050, 0.015, 0.015, 0.015];
% Delta_asymm_list = [0.0, 0.0, 0.0, 0.0, 0.02, -0.02];
% angle_list = [1.12, 1.12, 1.12, 1.12, 1.12, 1.12];
% 
% opt = 3;
% gamma0 = gamma0_list(opt);
% gamma1 = gamma1_list(opt);
% gamma3 = gamma3_list(opt);
% gamma4 = gamma4_list(opt);
% delta_dimer = delta_dimer_list(opt);
% Delta_asymm = Delta_asymm_list(opt);
% flag_chiral = true;  % true: AB-AB / false: AB-BA
% 
% w_aa = 0.0797;
% w_ab = 0.0975;
% theta = angle_list(opt) / 180 * pi; % 角度制要转化为弧度制 0.9 / 1.05 / 1.1 / 1.33 / 1.5 / 2.0 

%% parameters from Competing ground states and abundant orbital magnetism in twisted monolayer-bilayer graphene
gamma0 = 2.610;
gamma1 = 0.361;
gamma3 = 0.283;
gamma4 = 0.140;
delta_dimer = 0;
Delta_asymm = -0.04;
flag_chiral = true;  % true: AB-AB / false: AB-BA

w_aa = 0.110 * 0.5;
w_ab = 0.110;
theta = 1.24 / 180 * pi; % 角度制要转化为弧度制 0.9 / 1.05 / 1.1 / 1.33 / 1.5 / 2.0 

%% 第一步：作截断  trunction of k vecs
q_trunc = 3; % 
[bm_sup_mn, bm_sup_vecs] = get_q_couple(theta, q_trunc);
% 我们可以得到维数：(不包括自旋和谷的自由度),因为TBG在正常情况下都是自旋和谷都简并的
num_q_couple = size(bm_sup_vecs, 2);
dims = 3 * 2 * num_q_couple;

%% 第二步：高对称线
ak_norm = 8 * pi * abs(sin(theta / 2)) / 3; 
k_vertice_list = [[0; 0],[sqrt(3)/2; 1/2], [0; 1/2], [0; 1]];
np_side_list = [101, 101, 101];
[ak_len_array, akx_array, aky_array] = brillouin_k_line(k_vertice_list, np_side_list, ak_norm);

%% 第三步：构造哈密顿量 @ each k
eig_vals_K = zeros(num_k, dims);
eig_vals_Kp = zeros(num_k, dims);

% 我只关心能量最低的几条能带，对它们进行连续演化
low_ene_bound = 0.1;

tic
for i = 1:num_k
% for i = 1:2
    akx = akx_array(i);
    aky = aky_array(i);
    [H_tmbg_K, H_tmbg_Kp] = construct_tmbg_continuum_model(gamma0, gamma1, gamma3, gamma4, delta_dimer, Delta_asymm, w_aa, w_ab, ...
                                                            akx, aky, theta, num_q_couple, bm_sup_vecs, bm_sup_mn);
    
    % % valley K
    [eigvecs_K, eigvals] = eig(H_tmbg_K);
    eig_vals_K_temp = diag(eigvals);
    eig_vals_K(i, :) = eig_vals_K_temp;
    
    % 对低能能带进行排列
    if i == 1 % 第一个k点
        [~, new_order] = sort(abs(eig_vals_K_temp));
        count = 0;
        for jj = 1:dims
            index = new_order(jj);
            if abs(eig_vals_K_temp(index)) < low_ene_bound
                count = count + 1;
            end
        end
        
        % 数低能能带的数目
        num_K_low_band = count;
        if num_K_low_band == 0
            % 如果低能能带的数目等于0则无法进行排序
            continue
        end
        
        % 定义存储low_ene的几个变量
        eig_vals_K_low_ene = zeros(num_k, num_K_low_band);
        eig_vecs_K_low_ene = zeros(num_k, dims, num_K_low_band);  % 知道本征态可以用来计算Berry curvature / orbital magnetism / Chern number
        eig_vecs_K_low_ene_last = zeros(1, dims, num_K_low_band);
        eig_vals_K_low_ene_last = zeros(1, num_K_low_band);
        
        % 更新low_ene_last
        eig_vecs_K_low_ene_last(1,:,:) = eigvecs_K(:, new_order(1:num_K_low_band));
        eig_vals_K_low_ene_last(1,:) = eig_vals_K_temp(new_order(1:num_K_low_band));
        
        % 存入low_ene
        eig_vals_K_low_ene(i, :) = eig_vals_K_low_ene_last;
        eig_vecs_K_low_ene(i, :, :) = eig_vecs_K_low_ene_last;
    else
        if num_K_low_band == 0
            continue
        end
        
        % 寻找low_ene
        eps = sqrt(0.5); % 取为sqrt(0.5)的原因是当存在双重简并时，肯定有条能带的内积是大于sqrt(0.5)，而另一条小于sqrt(0.5)，但是对于多重简并就不成立了

        order_now = [];
        for jj = 1:num_K_low_band
            vec_last = eig_vecs_K_low_ene_last(1, :, jj);
            val_last = eig_vals_K_low_ene_last(1,jj);

            degeneracy = 0;
            for ii = 1:dims
                vec_now = eigvecs_K(:, ii);
                val_now = eig_vals_K_temp(ii);

                if abs(val_now - val_last) < 0.01
                    inner_dot = abs(dot(vec_now, vec_last));

                    if inner_dot >= eps
                        if ismember(ii, order_now)
                            continue
                        else
                            order_now = [order_now, ii];
                        end
                    end
                end
            end
            % disp(degeneracy)
        end

        % 更新low_ene_last
        for jj = 1:num_K_low_band  
            % 有时候order_now中的数目可能比num_K_low_band多，但我们只要保证能量最低的那些态是连续的就行
            eig_vecs_K_low_ene_last(1,:,jj) = eigvecs_K(:, order_now(jj));
            eig_vals_K_low_ene_last(1,jj) = eig_vals_K_temp(order_now(jj));
        end

        % 存入low_ene
        eig_vals_K_low_ene(i, :) = eig_vals_K_low_ene_last;
        eig_vecs_K_low_ene(i, :, :) = eig_vecs_K_low_ene_last;
    end
    
    % % valley Kp
    [eigvecs_Kp, eigvals] = eig(H_tmbg_Kp);
    eig_vals_Kp_temp = diag(eigvals);
    eig_vals_Kp(i, :) = eig_vals_Kp_temp;
    
    % 对低能能带进行排列
    if i == 1 % 第一个k点
        [~, new_order] = sort(abs(eig_vals_Kp_temp));
        count = 0;
        for jj = 1:dims
            index = new_order(jj);
            if abs(eig_vals_Kp_temp(index)) < low_ene_bound
                count = count + 1;
            end
        end
        
        % 数低能能带的数目
        num_Kp_low_band = count;
        if num_Kp_low_band == 0
            % 如果低能能带的数目等于0则无法进行排序
            continue
        end
        
        % 定义存储low_ene的几个变量
        eig_vals_Kp_low_ene = zeros(num_k, num_Kp_low_band);
        eig_vecs_Kp_low_ene = zeros(num_k, dims, num_Kp_low_band);  % 知道本征态可以用来计算Berry curvature / orbital magnetism / Chern number
        eig_vecs_Kp_low_ene_last = zeros(1, dims, num_Kp_low_band);
        eig_vals_Kp_low_ene_last = zeros(1, num_Kp_low_band);
        
        % 更新low_ene_last
        eig_vecs_Kp_low_ene_last(1,:,:) = eigvecs_Kp(:, new_order(1:num_Kp_low_band));
        eig_vals_Kp_low_ene_last(1,:) = eig_vals_Kp_temp(new_order(1:num_Kp_low_band));
        
        % 存入low_ene
        eig_vals_Kp_low_ene(i, :) = eig_vals_Kp_low_ene_last;
        eig_vecs_Kp_low_ene(i, :, :) = eig_vecs_Kp_low_ene_last;
    else
        if num_Kp_low_band == 0
            continue
        end
        
        % 寻找low_ene
%         eps = sqrt(0.5); %
% 
%         order_now = [];
%         for jj = 1:num_Kp_low_band
%             vec_last = eig_vecs_Kp_low_ene_last(1, :, jj);
%             val_last = eig_vals_Kp_low_ene_last(1,jj);
% 
%             degeneracy = 0;
%             for ii = 1:dims
%                 vec_now = eigvecs_Kp(:, ii);
%                 val_now = eig_vals_Kp_temp(ii);
% 
%                 if abs(val_now - val_last) < 0.01
%                     inner_dot = abs(dot(vec_now, vec_last));
% 
%                     if inner_dot >= eps
%                         if ismember(ii, order_now)
%                             continue
%                         else
%                             order_now = [order_now, ii];
%                         end
%                     end
%                 end
%             end
%             % disp(degeneracy)
%         end
% 
%         % 更新low_ene_last
%         for jj = 1:num_Kp_low_band  
%             % 有时候order_now中的数目可能比num_Kp_low_band多，但我们只要保证能量最低的那些态是连续的就行
%             eig_vecs_Kp_low_ene_last(1,:,jj) = eigvecs_Kp(:, order_now(jj));
%             eig_vals_Kp_low_ene_last(1,jj) = eig_vals_Kp_temp(order_now(jj));
%         end
% 
%         % 存入low_ene
%         eig_vals_Kp_low_ene(i, :) = eig_vals_Kp_low_ene_last;
%         eig_vecs_Kp_low_ene(i, :, :) = eig_vecs_Kp_low_ene_last;
    end
end
toc

% ene_shift = eig_vals_K_low_ene(1); % 将低能的K简并点平移到零能处\
ene_shift = 0.0;

figure
plot(ak_len_array, eig_vals_K - ene_shift, 'b-','LineWidth', 2)
hold on;
plot(ak_len_array, eig_vals_Kp - ene_shift, 'r--','LineWidth', 1.5)
grid on
xlim([min(ak_len_array), max(ak_len_array)])
%% parameters from Topological flat bands in twisted trilayer graphene
% ylim_upper_list = [0.02, 0.05, 0.04, 0.04, 0.03, 0.05];
% ylim_lower_list = [-0.02, -0.05, -0.04, -0.02, -0.045, -0.02];
% ylim([ylim_lower_list(opt), ylim_upper_list(opt)]) % 30meV

%% parameters from Competing ground states and abundant orbital magnetism in twisted monolayer-bilayer graphene
ylim([-0.125, 0.175])
%% 
% figure
% plot(ak_len_array, eig_vals_K_low_ene - ene_shift, 'b-','LineWidth', 2)
% hold on
% plot(ak_len_array, eig_vals_Kp_low_ene - ene_shift, 'r--','LineWidth', 1.5)
% ylim([-0.004, 0.004]) % 30meV
% grid on
% xlim([min(ak_len_array), max(ak_len_array)])
% ylim([min([min(eig_vals_K_low_ene), min(eig_vals_Kp_low_ene)]) - 0.01, max([max(eig_vals_K_low_ene), max(eig_vals_Kp_low_ene)]) + 0.01]) % 30meV