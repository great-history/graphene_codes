%% 一些参数
onsite_a = 0.015; % 0.015 / 0.0
onsite_b = 0.0;
% vf = 2.1354; % eV为单位
gamma0 = 2.4657;
w_aa = 0.0797;
w_ab = 0.0975;
theta = 1.05 / 180 * pi; % 角度制要转化为弧度制

%% 第一步：作截断  trunction of k vecs
q_trunc = 3; % 
[bm_sup_mn, bm_sup_vecs] = get_q_couple(theta, q_trunc);
% 我们可以得到维数：(不包括自旋和谷的自由度),因为TBG在正常情况下都是自旋和谷都简并的
num_q_couple = size(bm_sup_vecs, 2);
dims = 2 * 2 * num_q_couple;

%% 第二步：Brillouin Zone (G1 / G2)
ak_norm = 8 * pi * abs(sin(theta / 2)) / 3; 
G1_vec = ak_norm * [sqrt(3) / 2; - 3 / 2];
G2_vec = ak_norm * [sqrt(3) / 2; 3 / 2];

% 得到k点的数目
num_k1 = 51;
num_k2 = 101;
num_k = num_k1 * num_k2;

% 得到akx, aky_array
akx_array = zeros(num_k1, num_k2);
aky_array = zeros(num_k1, num_k2);
for jj = 1:num_k2
    for ii = 1:num_k1
        akx_array(ii, jj) = (ii - 1) / (num_k1 - 1) * G1_vec(1) + (jj - 1) / (num_k2 - 1) * G2_vec(1);
        aky_array(ii, jj) = (ii - 1) / (num_k1 - 1) * G1_vec(2) + (jj - 1) / (num_k2 - 1) * G2_vec(2);
    end
end

%% 计算能带结构
wilson_loop_array = zeros(num_k1, 1);

k1_index = 1;
eig_vals_K = zeros(num_k2, dims);

% 我只关心能量最低的几条能带，对它们进行连续演化
low_ene_bound = 0.33;

% 定义存储low_ene的几个变量
eig_vals_K_low_ene_cell = cell(num_k2, 1); % 存放低能的本征值
eig_vecs_K_low_ene_cell = cell(num_k2, 1); % 存放低能的本征态 ： 知道本征态可以用来计算Berry curvature / orbital magnetism / Chern number
num_band_list = zeros(num_k2, 1);

%% 每个k对应的哈密顿量进行对角化并找到指定能量范围内的所有态
tic
for k2_index = 1:num_k2
    %% 对(akx, aky)的哈密顿量进行对角化
    akx = akx_array(k1_index, k2_index);
    aky = aky_array(k1_index, k2_index);
    [H_tbg_K, H_tbg_Kp] = construct_tbg_continuum_model(onsite_a, onsite_b, gamma0, w_aa, w_ab, akx, aky, theta, num_q_couple, bm_sup_vecs, bm_sup_mn);
    [eig_vals_K, eig_vals_K_low_ene_temp, eig_vecs_K_low_ene_temp, num_K_low_band] = calc_ham_and_get_low_ene(H_tbg_K, eig_vals_K, low_ene_bound, dims, k2_index);
    
    %% 存入low_ene
    eig_vals_K_low_ene_cell{k2_index} = eig_vals_K_low_ene_temp;
    eig_vecs_K_low_ene_cell{k2_index} = eig_vecs_K_low_ene_temp;
    num_band_list(k2_index) = length(eig_vals_K_low_ene_temp);
end
toc

%% 确定准简并子空间
tic
ene_eps = 5e-4; 
[subspace_index_cell, subspace_range_cell] = get_quasi_deg_subspace(eig_vals_K_low_ene_cell, num_band_list, num_k2, ene_eps);
toc

%% 确定彼此有能隙的子空间
tic
overlap_eps = 0.1; % 如果重叠小于overlap_eps的话就被认为是有重叠的，一般overlap_eps取为0.1 / 0.01，如果取1e-3 / 1e-4则会把所有子空间都Mix在一起
[eig_vals_K_low_ene_gssp_cache, eig_vecs_K_low_ene_gssp_cache, dim_gssp_list, num_gssp] = get_gssp(eig_vecs_K_low_ene_cell, eig_vals_K_low_ene_cell, num_band_list, ...
                                                                                                   subspace_index_cell, subspace_range_cell, num_k2, dims, overlap_eps);
toc

%% 清楚没必要的变量
% clear eig_vals_K_low_ene_cell
% clear eig_vals_K_low_ene_temp
% clear eig_vecs_K_low_ene_cell
% clear subspace_index_cell
% clear subspace_range_cell

%% 在每个子空间中找出能带（两个判据：能量连续 + 群速度连续）
% 首先确定能量绝对值最低的两条能带(平带)在哪个gssp中
eig_vals_K_low_ene_gssp_temp = eig_vals_K_low_ene_gssp_cache{1}(1,:);
eig_vals_min = abs(sum(eig_vals_K_low_ene_gssp_temp) / length(eig_vals_K_low_ene_gssp_temp));
flat_band_gssp_index = 1;

for gssp_index = 2:num_gssp
    eig_vals_K_low_ene_gssp_temp = eig_vals_K_low_ene_gssp_cache{gssp_index}(1, :);
    eig_vals_temp = sum(eig_vals_K_low_ene_gssp_temp) / length(eig_vals_K_low_ene_gssp_temp);
    eig_vals_temp = abs(eig_vals_temp);
    
    if eig_vals_temp < eig_vals_min
        eig_vals_min = eig_vals_temp;
        flat_band_gssp_index = gssp_index;
    end
end

%% 计算每个k2下的两条平带对应的wilson loop
phase_accum = 0;

vec_last = eig_vecs_K_low_ene_gssp_cache{flat_band_gssp_index}(1, :);
vec_last = transpose(vec_last);
count_phase_change = 0;

for k_index = 1:num_k2
    if (k_index + 1) <= num_k2
        vec_temp = eig_vecs_K_low_ene_gssp_cache{flat_band_gssp_index}(k_index + 1, :);
        vec_temp = transpose(vec_temp);
    else
        vec_temp = eig_vecs_K_low_ene_gssp_cache{flat_band_gssp_index}(1, :);
        vec_temp = transpose(vec_temp);
    end
    
    % U(1) link
    u1_link = dot(vec_last, vec_temp);
    u1_link = u1_link / norm(u1_link);
    phase = angle(u1_link);
    if abs(pi - abs(phase)) < 0.1 % 与pi或-pi太接近，这时候容易多出一个2*pi的误差
        count_phase_change = count_phase_change + 1;
        
        vec_temp = - vec_temp; % 取负号就是加一个pi的相位
        u1_link = dot(vec_last, vec_temp);
        u1_link = u1_link / norm(u1_link);
        phase = angle(u1_link);
    end
    
    % 累加
    phase_accum = phase_accum + phase;
    
    vec_last = vec_temp;
end

%% 作图
figure
plot(eig_vals_K, 'b-','LineWidth', 2)
ylim([-0.004, 0.004]) % 30meV
grid on
xlim([1, num_k2])
ylim([-0.02, 0.02]) % 30meV % ylim([-0.2, 0.2]) % 30meV

% flat band
figure
for gssp_index = 1:num_gssp
    if gssp_index == flat_band_gssp_index
        plot(eig_vals_K_low_ene_gssp_cache{gssp_index}, 'r-','LineWidth', 2)
    else
        plot(eig_vals_K_low_ene_gssp_cache{gssp_index}, '-','LineWidth', 2)
    end
    hold on
end

grid on
xlim([1, num_k2])
ylim([-0.02, 0.02]) % ylim([-low_ene_bound, low_ene_bound]) % 30meV
