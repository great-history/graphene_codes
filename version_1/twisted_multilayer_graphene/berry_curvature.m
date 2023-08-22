addpath("brillouin_zone\")
%% 一些参数
onsite_a = 0.015; % 0.015
onsite_b = 0.0;

vf = 2.1354; % eV为单位
gamma0 = vf * 2 / sqrt(3);
gamma1 = 0.4;
gamma3 = 0.32;
gamma4 = 0.044;
delta_dimer = 0.050;
Delta_asymm = 0.0; % 0.0 / 0.005 / 0.020 / -0.005 / -0.020
flag_chiral = true; % true: AB-AB / false: AB-BA

w_aa = 0.0797;
w_ab = 0.0975;
theta = 2.0 / 180 * pi; % 角度制要转化为弧度制 0.9 / 1.05 / 1.1 / 1.33 / 1.5 / 2.0 

%% 第一步：作截断  trunction of k vecs
q_trunc = 3; % 
[bm_sup_mn, bm_sup_vecs] = get_q_couple(theta, q_trunc);
% 我们可以得到维数：(不包括自旋和谷的自由度),因为TBG在正常情况下都是自旋和谷都简并的
num_q_couple = size(bm_sup_vecs, 2);
dims = 4 * 2 * num_q_couple;

%% 首先对第一布里渊区进行撒点
np_side = 5;
[bz_bound, bz_inner, k_K_list] = brillouin_sampling(np_side);

%% 对每个k点进行计算
num_k_bound = size(bz_bound, 2);
num_k_inner = size(bz_inner, 2);
num_k = num_k_bound + num_k_inner;

%% 第三步：构造哈密顿量 @ each k
eig_vals_K = zeros(num_k, dims);
eig_vals_Kp = zeros(num_k, dims);

% 我只关心能量最低的几条能带，对它们进行连续演化
low_ene_bound = 0.2;

tic
for i = 1:num_k
    if i <= num_k_inner % 布里渊区内部的k点
        akx = bz_inner(1,i);
        aky = bz_inner(2,i);
    else % 布里渊区边界上的k点
        akx = bz_bound(1, i - num_k_inner);
        aky = bz_bound(2, i - num_k_inner);
    end
    [H_tdbg_K, H_tdbg_Kp] = construct_tdbg_continuum_model(gamma0, gamma1, gamma3, gamma4, delta_dimer, Delta_asymm, flag_chiral, w_aa, w_ab, ...
                                                            akx, aky, theta, num_q_couple, bm_sup_vecs, bm_sup_mn);
    
    % % valley K
    [eigvecs_K, eigvals] = eig(H_tdbg_K);
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
        
        
    end
    
    % % valley Kp
    [eigvecs_Kp, eigvals] = eig(H_tdbg_Kp);
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
    end
end
toc


%% 将所有k点的贝里曲率进行求和，最后还有乘dkx * dky，得到陈数

