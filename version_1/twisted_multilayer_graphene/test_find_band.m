addpath("brillouin_zone\")
addpath("helper\")
addpath("multilayer_graphene\")
addpath("velocity_operator\")
addpath("utilities\")

%% 一些参数
onsite_a = 0.010; % 0.015 / 0.0
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

%% 第二步：高对称线
% % 从K到Gamma再到M再到K' ：k_vertice_list = [[0; 0],[sqrt(3)/2; 1/2], [0; 1/2], [0; 1]]; np_side_list = [101, 101, 101];
% % 从K到K'到Gamma再到K ： k_vertice_list = [[0; 0], [0; 3]]; np_side_list = [301];
% % 从K到K'再到K ： 

% ak_norm = 8 * pi * abs(sin(theta / 2)) / 3; 
% k_vertice_list = [[0; 0], [0; 1], [sqrt(3)/2; 3/2]]; 
% np_side_list = [101, 101]; % 
% k_direction_list = [[0; 1], [sqrt(3)/2; 1/2]];

ak_norm = 8 * pi * abs(sin(theta / 2)) / 3; 
k_vertice_list = [[0; 0], [0; 1], [sqrt(3)/2; 3/2], [sqrt(3); 1], [sqrt(3); 0], [sqrt(3)/2; -1/2], [0;0]]; 
np_side_list = [101, 101, 101, 101, 101, 101]; % 
k_direction_list = [[0; 1], [sqrt(3)/2; 1/2]];

% k_vertice_list = [[0; 0],[sqrt(3)/2; 1/2], [0; 1/2], [0; 1]]; 
% np_side_list = [101, 101, 101];
% k_direction_list = [[sqrt(3)/2; 1/2], [-1; 0], [0; 1]];

[ak_len_array, akx_array, aky_array] = brillouin_k_line(k_vertice_list, np_side_list, ak_norm);
% 得到k点的数目
akx_array = akx_array(1:end-1); % 去掉最后一个k点
aky_array = aky_array(1:end-1); % 去掉最后一个k点
ak_len_array = ak_len_array(1:end-1); % 去掉最后一个k点
num_k = size(akx_array, 2);


%% 第三步：构造哈密顿量 @ each k
eig_vals_K = zeros(num_k, dims);

% 我只关心能量最低的几条能带，对它们进行连续演化
low_ene_bound = 0.31;

% 定义存储low_ene的几个变量
eig_vals_K_low_ene_cell = cell(num_k, 1); % 存放低能的本征值
eig_vecs_K_low_ene_cell = cell(num_k, 1); % 存放低能的本征态 ： 知道本征态可以用来计算Berry curvature / orbital magnetism / Chern number
num_band_list = zeros(num_k, 1);

%% 每个k对应的哈密顿量进行对角化并找到指定能量范围内的所有态
tic
for i = 1:num_k
    %% 对(akx, aky)的哈密顿量进行对角化
    akx = akx_array(i);
    aky = aky_array(i);
    [H_tbg_K, H_tbg_Kp] = construct_tbg_continuum_model(onsite_a, onsite_b, gamma0, w_aa, w_ab, akx, aky, theta, num_q_couple, bm_sup_vecs, bm_sup_mn);
    [eig_vals_K, eig_vals_K_low_ene_temp, eig_vecs_K_low_ene_temp, num_K_low_band] = calc_ham_and_get_low_ene(H_tbg_K, eig_vals_K, low_ene_bound, dims, i);
    
    %% 存入low_ene
    eig_vals_K_low_ene_cell{i} = eig_vals_K_low_ene_temp;
    eig_vecs_K_low_ene_cell{i} = eig_vecs_K_low_ene_temp;
    num_band_list(i) = length(eig_vals_K_low_ene_temp);
end
toc

%% 确定准简并子空间
tic
ene_eps = 5e-4; 
[subspace_index_cell, subspace_range_cell] = get_quasi_deg_subspace(eig_vals_K_low_ene_cell, num_band_list, num_k, ene_eps);
toc

%% 确定彼此有能隙的子空间
tic
overlap_eps = 0.01; % 如果重叠小于overlap_eps的话就被认为是有重叠的，一般overlap_eps取为0.1 / 0.01，如果取1e-3 / 1e-4则会把所有子空间都Mix在一起
[eig_vals_K_low_ene_gssp_cache, eig_vecs_K_low_ene_gssp_cache, dim_gssp_list, num_gssp] = get_gssp(eig_vecs_K_low_ene_cell, eig_vals_K_low_ene_cell, num_band_list, ...
                                                                                                            subspace_index_cell, subspace_range_cell, num_k, dims, overlap_eps);
toc
%% 清楚没必要的变量
clear eig_vals_K_low_ene_cell
clear eig_vals_K_low_ene_temp
clear eig_vecs_K_low_ene_cell
clear subspace_index_cell
clear subspace_range_cell
clear i

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

% 找出这些gssp中的所有简并点
% 根据群速度连续来找出能带(只能在同一方向上，如果改变了方向则不适用)

%% 计算陈数(Chern number)
Chern_number = 0;
phase_array = zeros(num_k, 1);

vec_last = eig_vecs_K_low_ene_gssp_cache{flat_band_gssp_index}(1, :);
vec_last = transpose(vec_last);
count_phase_change = 0;
for k_index = 1:num_k
    if (k_index + 1) <= num_k
        vec_temp = eig_vecs_K_low_ene_gssp_cache{flat_band_gssp_index}(k_index + 1, :);
        vec_temp = transpose(vec_temp);
    else
        vec_temp = eig_vecs_K_low_ene_gssp_cache{flat_band_gssp_index}(1, :);
        vec_temp = transpose(vec_temp);
    end
    
    % U(1) link
    u1_link = dot(vec_last, vec_temp);
    u1_link = u1_link / norm(u1_link);
    phase = imag(log(u1_link));
    if abs(pi - abs(phase)) < 0.1 % 与pi或-pi太接近，这时候容易多出一个2*pi的误差
        count_phase_change = count_phase_change + 1;
        
        vec_temp = - vec_temp; % 取负号就是加一个pi的相位
        u1_link = dot(vec_last, vec_temp);
        u1_link = u1_link / norm(u1_link);
        phase = imag(log(u1_link));
    end
    
    phase_array(k_index) = phase;
    Chern_number = Chern_number + phase;
    
    vec_last = vec_temp;
end

%% 作图
figure
plot(ak_len_array, eig_vals_K, 'b-','LineWidth', 2)
ylim([-0.004, 0.004]) % 30meV
grid on
xlim([min(ak_len_array), max(ak_len_array)])
ylim([-0.02, 0.02]) % 30meV % ylim([-0.2, 0.2]) % 30meV

% flat band 
figure
for gssp_index = 1:num_gssp
    if gssp_index == flat_band_gssp_index
        plot(ak_len_array, eig_vals_K_low_ene_gssp_cache{gssp_index}, 'r-','LineWidth', 2)
    else
        plot(ak_len_array, eig_vals_K_low_ene_gssp_cache{gssp_index}, '-','LineWidth', 2)
    end
    hold on
end

grid on
xlim([min(ak_len_array), max(ak_len_array)])
ylim([-0.02, 0.02]) % ylim([-low_ene_bound, low_ene_bound]) % 30meV