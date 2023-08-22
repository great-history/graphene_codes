%% 添加路径
addpath('.\bilayer_hamitonians\')
addpath('.\utils\')
%% 参数设置
gamma0 = 2.8; % 完全pi键
gamma1 = - 0.4; % 完全sigma键
gamma3 = 0.3;
gamma4 = 0.04;
delta_dimer = 0.022;
U = 0.015; % 0.05

cc_bond1 = [sqrt(3)/2;1/2];
cc_bond2 = [-sqrt(3)/2;1/2];
cc_bond3 = [0;-1];

%% 对布里渊区进行撒点
np_side = 100;
ak_norm = 4 * pi / (sqrt(3) * 3);
[bz_bound, bz_inner, k_K_list] = brillouin_sampling_multilayer_graphene(np_side, ak_norm);

%% 根据y坐标对bz_inner中的点进行排序
bz_point_cell = cell(2 * np_side - 1, 1);
y_coord_list = ak_norm * linspace(sqrt(3)/2, -sqrt(3)/2, 2 * np_side - 1);
y_diff = (y_coord_list(1) - y_coord_list(2));

num_k_inner = size(bz_inner, 2);
num_k_bound = size(bz_bound, 2);
num_k = num_k_bound + num_k_inner;

for i = 1:(2 * np_side - 1)
    y_now = y_coord_list(i);
    bz_point_list = [];
    % 对内部的点进行筛选
    for j = 1:num_k_inner
        if abs(bz_inner(2, j) - y_now) < y_diff / 100
            bz_point_list = [bz_point_list, bz_inner(:,j)];
        end
    end
    
    % 对外部的点进行筛选
    for j = 1:num_k_bound
        if abs(bz_bound(2, j) - y_now) < y_diff / 100
            bz_point_list = [bz_point_list, bz_bound(:,j)];
        end
    end
    
    % 将bz_point_list按照x坐标进行排序
    [~, new_order] = sort(bz_point_list(1,:));
    bz_point_list = bz_point_list(:, new_order);
    
    bz_point_cell{i} = bz_point_list;
end

%% 计算哈密顿量
eig_val_cell = cell(2 * np_side - 1, 1);
eig_vec_cell = cell(2 * np_side - 1, 1);

tic
for i = 1:(2 * np_side - 1)
    num_k_temp = length(bz_point_cell{i});
    k_point_list = bz_point_cell{i};
    eig_val_array_temp = zeros(num_k_temp, 4);
    eig_vec_array_temp = zeros(num_k_temp, 4, 4);
    
    for k_index = 1:num_k_temp
        akx = k_point_list(1, k_index);
        aky = k_point_list(2, k_index);
        
        H_ham = construct_bilayer_tight_binding_model(gamma0, gamma1, gamma3, gamma4, delta_dimer, akx, aky);
        H_ham = construct_multilayer_graphene_onsite(U, 0, akx, aky, 2, H_ham, "tight_binding");
        % hem = helper_check_hermite(H_ham, 1e-8);
        % if hem == 0
        %   disp("monolayer Ham is not hermitian")
        % end
        % 对角化
        % valley K
        [eig_vecs_temp, eig_vals_temp] = eig(H_ham);
        eig_vals_temp = diag(eig_vals_temp);

        % 按照绝对值重新进行排序, 可以减少计算次数
        [~, new_order] = sort(eig_vals_temp);
        
        eig_val_array_temp(k_index, :) = eig_vals_temp(new_order);
        eig_vec_array_temp(k_index, :, :) = eig_vecs_temp(:, new_order);
    end
    eig_val_cell{i} = eig_val_array_temp;
    eig_vec_cell{i} = eig_vec_array_temp;
end
toc

%% 作图
% for i = 1:(2 * np_side - 1)
%     figure
%     plot(eig_val_cell{i}, 'b-','LineWidth', 2)
% end
% 
figure
plot(eig_val_cell{1}, 'b-','LineWidth', 2)
hold on
plot(eig_val_cell{end}, 'r--','LineWidth', 2)