%% 计算在mBZ上的tbg能带(三维图)
addpath("brillouin_zone\")
addpath("helper\")
addpath("multilayer_graphene\")
addpath("velocity_operator\")

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

%% 第二步：整个第一布里渊区撒点
np_side = 101;
ak_norm = 8 * pi * abs(sin(theta / 2)) / 3; 
[bz_bound, bz_inner, k_K_list] = brillouin_sampling_twisted_multilayer_graphene(np_side, ak_norm);

%% 根据x坐标对bz_inner中的点进行排序
bz_point_cell = cell(2 * np_side - 1, 1);
x_coord_list = ak_norm * linspace(-sqrt(3)/2, sqrt(3)/2, 2 * np_side - 1);
x_diff = abs(x_coord_list(1) - x_coord_list(2));

num_k_inner = size(bz_inner, 2);
num_k_bound = size(bz_bound, 2);
num_k = num_k_bound + num_k_inner;

for i = 1:(2 * np_side - 1)
    x_now = x_coord_list(i);
    bz_point_list = [];
    % 对内部的点进行筛选
    for j = 1:num_k_inner
        if abs(bz_inner(1, j) - x_now) < x_diff / 100
            bz_point_list = [bz_point_list, bz_inner(:,j)];
        end
    end
    
    % 对外部的点进行筛选
    for j = 1:num_k_bound
        if abs(bz_bound(1, j) - x_now) < x_diff / 100
            bz_point_list = [bz_point_list, bz_bound(:,j)];
        end
    end
    
    % 将bz_point_list按照x坐标进行排序
    [~, new_order] = sort(bz_point_list(2,:));
    bz_point_list = bz_point_list(:, new_order);
    
    bz_point_cell{i} = bz_point_list;
end

% 对于转角体系还要平移一下
num_lines = 2 * np_side - 1;
for line_index = 1:num_lines
    bz_point_cell{line_index} = bz_point_cell{line_index} + ak_norm * [sqrt(3)/2;1/2];
end

%% 不放心的话，可以画一下
% scatter(bz_bound(1,:), bz_bound(2,:),50,'r','p','filled');
% hold on
% scatter(bz_inner(1,1:end-1), bz_inner(2,1:end-1),30,'b','s');
% hold on
% scatter(bz_inner(1,end), bz_inner(2,end),30,'k','+');
% hold on

%% 第三步：计算之前的准备工作
num_lines = 2 * np_side - 1;
eig_vals_K_cell = cell(num_lines, 1);  % 每个cell存储的是某条沿着y方向的k_line上的能带(保存有所有的能量但是并没有将低能能带找出来)
low_ene_bound = 0.50; % 我只关心能量在[-low_ene_bound, low_ene_bound]范围内的几条能带，对它们进行连续演化
eig_vals_K_low_ene_gssp_cache_cell = cell(num_lines, 1); % 记录每条Line上选出来的低能能带的能量本征值
eig_vecs_K_low_ene_gssp_cache_cell = cell(num_lines, 1); % 记录每条Line上选出来的低能能带的能量本征态
dim_gssp_list_cell = cell(num_lines, 1); % 记录每条Line上存在能隙的子空间个数
num_gssp_list = zeros(num_lines, 1);

poolobj = gcp('nocreate');
if isempty(poolobj)
    disp('启动并行运算，核心数：8');
    % Perform a basic check by entering this code, where "local" is one kind of cluster profile.
    parpool('local', 8);
else
    disp(['并行运o算已启动，核心数：' num2str(poolobj.NumWorkers)]);
end

% 对每条Line进行计算
tic
parfor line_index = 1:num_lines
% for line_index = 1:num_lines
    num_k = size(bz_point_cell{line_index}, 2);
    akx_array = bz_point_cell{line_index}(1,:);
    aky_array = bz_point_cell{line_index}(2,:);
    
    eig_vals_K = zeros(num_k, dims);
    % 定义存储low_ene的几个变量
    eig_vals_K_low_ene_cell = cell(num_k, 1); % 存放低能的本征值
    eig_vecs_K_low_ene_cell = cell(num_k, 1); % 存放低能的本征态 ： 知道本征态可以用来计算Berry curvature / orbital magnetism / Chern number
    num_band_list = zeros(num_k, 1);
    
    %% 第四步：构造哈密顿量 @ each k : 每个k对应的哈密顿量进行对角化并找到指定能量范围内的所有态
    % tic
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
    % toc
    
    %% 确定准简并子空间
    % tic
    ene_eps = 5e-4; 
    [subspace_index_cell, subspace_range_cell] = get_quasi_deg_subspace(eig_vals_K_low_ene_cell, num_band_list, num_k, ene_eps);
    % toc
    
    %% 确定彼此有能隙的子空间
    % tic
    overlap_eps = 0.01; % 如果重叠小于overlap_eps的话就被认为是有重叠的，一般overlap_eps取为0.1 / 0.01，如果取1e-3 / 1e-4则会把所有子空间都Mix在一起
    [eig_vals_K_low_ene_gssp_cache, eig_vecs_K_low_ene_gssp_cache, dim_gssp_list, num_gssp] = get_gssp(eig_vecs_K_low_ene_cell, eig_vals_K_low_ene_cell, num_band_list, ...
                                                                                                                subspace_index_cell, subspace_range_cell, num_k, dims, overlap_eps);
    % toc
    
    %% 存到cell中
    eig_vals_K_cell{line_index} = eig_vals_K;
    eig_vals_K_low_ene_gssp_cache_cell{line_index} = eig_vals_K_low_ene_gssp_cache;
    eig_vecs_K_low_ene_gssp_cache_cell{line_index} = eig_vecs_K_low_ene_gssp_cache; % 记录每条Line上选出来的低能能带的能量本征态
    dim_gssp_list_cell{line_index} = dim_gssp_list; % 记录每条Line上存在能隙的子空间个数
    num_gssp_list(line_index) = num_gssp;
    
%     % 首先确定能量绝对值最低的两条能带(平带)在哪个gssp中
%     eig_vals_K_low_ene_gssp_temp = eig_vals_K_low_ene_gssp_cache{1}(1,:);
%     eig_vals_min = abs(sum(eig_vals_K_low_ene_gssp_temp) / length(eig_vals_K_low_ene_gssp_temp));
%     flat_band_gssp_index = 1;
% 
%     for gssp_index = 2:num_gssp
%         eig_vals_K_low_ene_gssp_temp = eig_vals_K_low_ene_gssp_cache{gssp_index}(1, :);
%         eig_vals_temp = sum(eig_vals_K_low_ene_gssp_temp) / length(eig_vals_K_low_ene_gssp_temp);
%         eig_vals_temp = abs(eig_vals_temp);
% 
%         if eig_vals_temp < eig_vals_min
%             eig_vals_min = eig_vals_temp;
%             flat_band_gssp_index = gssp_index;
%         end
%     end
    
    
end
toc

%% 如果不放心，可以做一下图看看
% figure
% plot(ak_len_array, eig_vals_K, 'b-','LineWidth', 2)
% ylim([-0.004, 0.004]) % 30meV
% grid on
% xlim([min(ak_len_array), max(ak_len_array)])
% ylim([-0.02, 0.02]) % 30meV % ylim([-0.2, 0.2]) % 30meV

% flat band
line_index_select_list = [1, 11, 21, 31, 41, 51];
for i = 1:6
    line_index = line_index_select_list(i);
    figure
    for gssp_index = 1:num_gssp_list(line_index)
        plot(bz_point_cell{line_index}(2,:), eig_vals_K_low_ene_gssp_cache_cell{line_index}{gssp_index}, '-','LineWidth', 2)
        hold on
    end

    grid on
    xlim([min(bz_point_cell{line_index}(2,:)), max(bz_point_cell{line_index}(2,:))])
    ylim([-0.02, 0.02]) % ylim([-low_ene_bound, low_ene_bound]) % 30meV 
end

%% 定出整个布里渊区上的能带图
for line_index = 1:num_lines
    
end