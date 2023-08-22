addpath("brillouin_zone\")
addpath("helper\")
addpath("multilayer_graphene\")
addpath("velocity_operator\")

%% 一些参数
onsite_a = 0.008; % 0.015 / 0.0
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

ak_norm = 8 * pi * abs(sin(theta / 2)) / 3; 
% k_vertice_list = [[0; 0], [0; 1], [sqrt(3)/2; 3/2]]; 
% np_side_list = [101, 101]; % 
k_vertice_list = [[0; 0],[sqrt(3)/2; 1/2], [0; 1/2], [0; 1]]; 
k_unit_list = [[sqrt(3)/2; 1/2], [-1; 0], [0; 1]];
np_side_list = [101, 101, 101];
[ak_len_array, akx_array, aky_array] = brillouin_k_line(k_vertice_list, np_side_list, ak_norm);
% 得到k点的数目
num_k = size(akx_array, 2);

%% 第三步：构造哈密顿量 @ each k
eig_vals_K = zeros(num_k, dims);

% 定义存储low_ene的几个变量
eig_vals_K_low_ene_cell_cell = {}; % 存放低能的本征值
eig_vecs_K_low_ene_cell_cell = {}; % 存放低能的本征态 ： 知道本征态可以用来计算Berry curvature / orbital magnetism / Chern number

eig_vals_K_low_ene_temp_cell = {}; % 存放当前k点的低能的本征值
eig_vecs_K_low_ene_temp_cell = {}; % 存放当前k点的低能的本征态
eig_vals_K_low_ene_last_cell = {}; % 存放上一个k点的低能的本征值
eig_vecs_K_low_ene_last_cell = {}; % 存放上一个k点的低能的本征态

num_band_list = [];
num_subspace_list = [];
quasi_deg_cell_cell = {}; % 近似简并：quasi degenerate
quasi_deg_temp_cell = {}; % 把简并的能带指标放到一起
quasi_deg_last_cell = {}; % 把简并的能带指标放到一起

% overlap_low_ene_temp = {};
overlap_low_ene_cell_cell = {}; % 存在temp与last之间的overlap

% 我只关心能量最低的几条能带，对它们进行连续演化
low_ene_bound = 0.31;
ene_eps = 1e-4; 

%% 开始计算i = 1
i = 1;
% 对低能能带进行排列
if i == 1 % 如果是第一个k点
    akx = akx_array(i);
    aky = aky_array(i);
    [H_tbg_K, H_tbg_Kp] = construct_tbg_continuum_model(onsite_a, onsite_b, gamma0, w_aa, w_ab, akx, aky, theta, num_q_couple, bm_sup_vecs, bm_sup_mn);
    [eig_vals_K, eig_vals_K_low_ene_temp, eig_vecs_K_low_ene_temp, num_K_low_band] = calc_tbg_ham(H_tbg_K, eig_vals_K, low_ene_bound, dims, i);
    
    % 数低能能带的数目
    if ~(num_K_low_band == 0)
        % 更新low_ene_last
        quasi_deg_last_cell{1}(1) = 1;
        eig_vals_K_low_ene_last_cell{1} = eig_vals_K_low_ene_temp(1);
        eig_vecs_K_low_ene_last_cell{1} = eig_vecs_K_low_ene_temp(:, 1);
        
        num_subspace = 1;
        for kk = 2:num_K_low_band
            val_temp = eig_vals_K_low_ene_temp(kk);
            vec_temp = eig_vecs_K_low_ene_temp(:, kk);
            
            if abs(val_temp - eig_vals_K_low_ene_last_cell{num_subspace}(end)) < ene_eps % 发生了近似简并
                quasi_deg_last_cell{num_subspace}(end + 1) = kk; % 存入能带指标kk
                eig_vals_K_low_ene_last_cell{num_subspace}(end + 1) = val_temp; % 存入近似简并的本征值
                eig_vecs_K_low_ene_last_cell{num_subspace} = [eig_vecs_K_low_ene_last_cell{num_subspace}, vec_temp]; % 存入近似简并的本征态
            else % 没有发生简并
                num_subspace = num_subspace + 1;
                quasi_deg_last_cell{num_subspace} = kk; % 存入能带指标kk
                eig_vals_K_low_ene_last_cell{num_subspace} = val_temp; % 存入新的本征值
                eig_vecs_K_low_ene_last_cell{num_subspace} = vec_temp; % 存入新的本征态q
            end
        end
        
        % 存入low_ene
        num_band_list = [num_band_list, num_K_low_band];
        num_subspace_list = [num_subspace_list, num_subspace];
        quasi_deg_cell_cell{i} = quasi_deg_last_cell;
        
        eig_vals_K_low_ene_cell_cell{i} = eig_vals_K_low_ene_last_cell;
        eig_vecs_K_low_ene_cell_cell{i} = eig_vecs_K_low_ene_last_cell;
    end
end

%% 开始计算第二到第num_k的k点
tic
for i = 2:num_k
    %% 对(akx, aky)的哈密顿量进行对角化
    akx = akx_array(i);
    aky = aky_array(i);
    [H_tbg_K, H_tbg_Kp] = construct_tbg_continuum_model(onsite_a, onsite_b, gamma0, w_aa, w_ab, akx, aky, theta, num_q_couple, bm_sup_vecs, bm_sup_mn);
    [eig_vals_K, eig_vals_K_low_ene_temp, eig_vecs_K_low_ene_temp, num_K_low_band] = calc_tbg_ham(H_tbg_K, eig_vals_K, low_ene_bound, dims, i);
    
    %% 计算overlap matrix
    overlap_low_ene = calc_overlap_matrix(eig_vals_K_low_ene_temp, eig_vecs_K_low_ene_temp, ...
                                          eig_vals_K_low_ene_last_cell, eig_vecs_K_low_ene_last_cell, ...
                                          num_K_low_band, num_subspace);
    
    %% 找出那些有overlap的态中的近似简并态
    innner_dot_threshold = 0.1; % 判断是否存在overlap的判据
    
    num_subspace = 0;
    num_K_low_band_temp = 0;
    
    quasi_deg_temp_cell = {};
    eig_vals_K_low_ene_temp_cell = {};
    eig_vecs_K_low_ene_temp_cell = {};
    
    overlap_low_ene_temp = [];
    % 首先寻找第一个态
    for jj = 1:num_K_low_band
        inner_dot_sum = sqrt(sum(overlap_low_ene(jj, :))); % 把所有的和加起来
        if inner_dot_sum > innner_dot_threshold % band_index对应的态与Last之间存在overlap
            num_subspace = num_subspace + 1;
            num_K_low_band_temp = num_K_low_band_temp + 1;
            overlap_low_ene_temp = [overlap_low_ene_temp; overlap_low_ene(jj, :)];
            
            quasi_deg_temp_cell{num_subspace} = num_K_low_band_temp;
            eig_vals_K_low_ene_temp_cell{num_subspace} = eig_vals_K_low_ene_temp(jj);
            eig_vecs_K_low_ene_temp_cell{num_subspace} = eig_vecs_K_low_ene_temp(:, jj);
            
            break
        end
    end
    
    % 然后寻找其它态
    % fprintf("第一个本征态的指标为%d\n", jj)
    for band_index = (jj + 1):num_K_low_band  % 从能量绝对值较低的态到较高的态
        vec_temp = eig_vecs_K_low_ene_temp(:, band_index);
        val_temp = eig_vals_K_low_ene_temp(band_index);
        
        inner_dot_sum = sqrt(sum(overlap_low_ene(band_index, :))); % 把所有的和加起来
        
        if inner_dot_sum > innner_dot_threshold % band_index对应的态与Last之间存在overlap
            num_K_low_band_temp = num_K_low_band_temp + 1;
            overlap_low_ene_temp = [overlap_low_ene_temp; overlap_low_ene(band_index, :)];
            
            if abs(val_temp - eig_vals_K_low_ene_temp_cell{num_subspace}(end)) < ene_eps % 发生了近似简并
                quasi_deg_temp_cell{num_subspace}(end + 1) = num_K_low_band_temp; % 存入能带指标
                eig_vals_K_low_ene_temp_cell{num_subspace}(end + 1) = val_temp; % 存入近似简并的本征值
                eig_vecs_K_low_ene_temp_cell{num_subspace} = [eig_vecs_K_low_ene_temp_cell{num_subspace}, vec_temp]; % 存入近似简并的本征态
            else % 没有发生简并
                num_subspace = num_subspace + 1;
                
                quasi_deg_temp_cell{num_subspace} = num_K_low_band_temp; % 存入能带指标
                eig_vals_K_low_ene_temp_cell{num_subspace} = val_temp; % 存入新的本征值
                eig_vecs_K_low_ene_temp_cell{num_subspace} = vec_temp; % 存入新的本征态q
            end
            
            % fprintf("找到了能带%d\n", band_index)
        end
    end
    num_K_low_band = num_K_low_band_temp;
    % fprintf("\n\n")
    
    %% 重新确定low_ene_bound
    % low_ene_bound = abs(eig_vals_K_low_ene_temp_cell{end}(end));
    % fprintf("当前low ene bound变为%f\n", low_ene_bound)
    
    %% overlap_low_ene_temp的行指标也变为subspace_index而不是band_index
    for subspace_index = length(quasi_deg_temp_cell):-1:1
        subspace = quasi_deg_temp_cell{subspace_index};
        if length(subspace) >= 2 % 存在简并
            overlap_low_ene_temp(subspace(1), :) = sum(overlap_low_ene_temp(subspace(1):subspace(end), :));
            overlap_low_ene_temp(subspace(1)+1:subspace(end), :) = [];
        end
    end
    
    %% 存入low_ene
    num_band_list = [num_band_list, num_K_low_band];
    num_subspace_list = [num_subspace_list, num_subspace];
    quasi_deg_cell_cell{i} = quasi_deg_temp_cell;
    
    overlap_low_ene_cell_cell{i} = sparse(overlap_low_ene_temp); % 由于0元素较多，用稀疏矩阵存储比较好
    eig_vals_K_low_ene_cell_cell{i} = eig_vals_K_low_ene_temp_cell;
    eig_vecs_K_low_ene_cell_cell{i} = eig_vecs_K_low_ene_temp_cell;
    
    % 将temp改成last
    quasi_deg_last_cell = quasi_deg_temp_cell;
    eig_vals_K_low_ene_last_cell = eig_vals_K_low_ene_temp_cell;
    eig_vecs_K_low_ene_last_cell = eig_vecs_K_low_ene_temp_cell;
end
toc

figure
plot(ak_len_array, eig_vals_K, 'b-','LineWidth', 2)
ylim([-0.004, 0.004]) % 30meV
grid on
xlim([min(ak_len_array), max(ak_len_array)])
ylim([-0.2, 0.2]) % 30meV
ylim([-0.02, 0.02]) % 30meV

%% 确定能带的演化
% 首先确定最少能带数
num_K_low_band = min(num_band_list);

% 然后确定最少能带数时的最多子空间数

% 首先进行subspace separation


% 然后进行band separation

% 使用群速度来进行能带区分，但这种方法处理Dirac Point时会失效，但处理一维的能带是没问题的


%% 对角化转角石墨烯哈密顿量的函数
function [eig_vals_K, eig_vals_K_low_ene_temp, eig_vecs_K_low_ene_temp, num_K_low_band] = calc_tbg_ham(H_tbg_K, eig_vals_K, low_ene_bound, dims, i)
    %% 对(akx, aky)的哈密顿量进行对角化
    % valley K
    [eig_vecs_K_temp, eigvals] = eig(H_tbg_K);
    eig_vals_K_temp = diag(eigvals);
    eig_vals_K(i, :) = eig_vals_K_temp;
    % 按照绝对值重新进行排序, 可以减少计算次数
    [~, new_order] = sort(abs(eig_vals_K_temp));
    eig_vals_K_temp = eig_vals_K_temp(new_order);
    eig_vecs_K_temp = eig_vecs_K_temp(:, new_order);
    
    %% 得到低能能带数num_K_low_band
    num_K_low_band = 0;
    for jj = 1:dims
        % if abs(eig_vals_K_temp(jj)) < low_ene_bound + 0.01  % 适用于low_ene_bound会自动调整的情况
        if abs(eig_vals_K_temp(jj)) < low_ene_bound  % 适用于low_ene_bound固定的情况
            num_K_low_band = num_K_low_band + 1;
        else
            break
        end
    end
    
    if ~(num_K_low_band == 0)
        eig_vals_K_low_ene_temp = eig_vals_K_temp(1:num_K_low_band);
        eig_vecs_K_low_ene_temp = eig_vecs_K_temp(:, 1:num_K_low_band);
        
        [~, new_order] = sort(eig_vals_K_low_ene_temp);
        eig_vals_K_low_ene_temp = eig_vals_K_low_ene_temp(new_order); % 不按绝对值从小到大排
        eig_vecs_K_low_ene_temp = eig_vecs_K_low_ene_temp(:, new_order);
    else
        eig_vals_K_low_ene_temp = [];
        eig_vecs_K_low_ene_temp = [];
    end
    
end

%% 计算overlap_matrix
function overlap_low_ene = calc_overlap_matrix(eig_vals_K_low_ene_temp, eig_vecs_K_low_ene_temp, eig_vals_K_low_ene_last_cell, ...
                                               eig_vecs_K_low_ene_last_cell, num_K_low_band, num_subspace)
    overlap_low_ene = zeros(num_K_low_band, num_subspace);
    for band_index = 1:num_K_low_band  % 从能量绝对值较低的态到较高的态
        val_temp = eig_vals_K_low_ene_temp(band_index);
        vec_temp = eig_vecs_K_low_ene_temp(:, band_index);
        
        % 通过计算overlap建立联系
        for deg_index = 1:num_subspace
            vecs_last = eig_vecs_K_low_ene_last_cell{deg_index}; % 可能存在简并的情况
            vals_last = eig_vals_K_low_ene_last_cell{deg_index}; % 可能存在简并的情况

            if abs(val_temp - vals_last(end)) < 0.01
                % 计算内积
                inner_dot = 0.0;
                for jj = 1:size(vecs_last, 2)
                    inner_dot = inner_dot + (abs(dot(vec_temp, vecs_last(:, jj))))^2;  % 加abs是因为有可能是复数
                end
                
                overlap_low_ene(band_index, deg_index) = inner_dot;
            end
        end
    end
end

function get_quasi_deg_subspaces()
    
end

%% 交换overlap_matrix中的行
function exchange_rows_by_subspaces(num_subspace, num_K_low_band)
    num_find = 0;
    for subspace_index = 1:num_subspace
        band_start = num_find + 1;
        
        for band_index = band_start:num_K_low_band
            if ~(overlap_low_ene(band_index, subspace_index) == 0)
                num_find = num_find + 1;
                if ~(band_index == num_find)
                    % 交换overlap matrix中的行
                    vec_overlap_temp = overlap_low_ene(band_index, :);
                    overlap_low_ene(band_index, :) = overlap_low_ene(num_find, :);
                    overlap_low_ene(num_find, :) = vec_overlap_temp;
    
                    % 交换eig_vecs_K_low_ene_temp中的行
                    vec_eig_temp = eig_vecs_K_low_ene_temp(:, band_index);
                    eig_vecs_K_low_ene_temp(:, band_index) = eig_vecs_K_low_ene_temp(:, num_find);
                    eig_vecs_K_low_ene_temp(:, num_find) = vec_eig_temp;
    
                    % 交换eig_vals_K_low_ene_temp中的行
                    val_eig_temp = eig_vals_K_low_ene_temp(band_index);
                    eig_vals_K_low_ene_temp(band_index) = eig_vals_K_low_ene_temp(num_find);
                    eig_vals_K_low_ene_temp(num_find) = val_eig_temp;
                end
            end
        end
    end
end

%% 合并两个order
function order_new = join_two_order(order1, order2)
    order_new = zeros(size(order1));
    for i = 1:length(order2)
        order_new(i) = order1(order2(i));
    end
end