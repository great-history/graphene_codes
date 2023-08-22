addpath("./bilayer_hamitonians/")
%% 参数设置
gamma0 = 3.16; % 完全pi键
gamma1 = - 0.381; % 完全sigma键
gamma3 = 0.38;
gamma4 = 0.14;
delta_dimer = 0.022;
U = 0.0; % 0.05

cc_bond1 = [sqrt(3)/2;1/2];
cc_bond2 = [-sqrt(3)/2;1/2];
cc_bond3 = [0;-1];

%% 对布里渊区进行撒点
np_side = 50;
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

%% 进行插值得到meshgrid(适用于在整个布里渊区上画三维能带图)
% 首先将所有k_points和eig_vals合并
k_point_list = [];
ene_list = [];
for i = 1:(2 * np_side - 1)
    k_point_list = [k_point_list, bz_point_cell{i}];
    ene_list = [ene_list, transpose(eig_val_cell{i})];
end

% 我们得到第一布里渊区的最近邻六个第二布里渊区

recip_vec_list = ak_norm * [[0;0], [3/2; sqrt(3)/2], [0; sqrt(3)], [-3/2; sqrt(3)/2], [-3/2; -sqrt(3)/2], [0; - sqrt(3)], [3/2; -sqrt(3)/2]];
for g_index = 2:length(recip_vec_list)
    k_point_list = [k_point_list, k_point_list(:, 1:num_k) + recip_vec_list(:, g_index)];
    ene_list = [ene_list, ene_list(:, 1:num_k)];
end

% 然后用griddata进行插值
[kx_q,ky_q] = meshgrid(-3/2:1/(10*np_side):3/2, -sqrt(3):1/(10*np_side):sqrt(3));
kx_q = ak_norm * kx_q;
ky_q = ak_norm * ky_q;
ene_q = cell(4,1);

for band_index = 1:4
    ene_q{band_index} = griddata(k_point_list(1,:), k_point_list(2,:), ene_list(band_index,:), kx_q, ky_q);
end

%% 作图
% figure
% for i = 1:(2 * np_side - 1)
%     plot(eig_val_cell{i}, 'b-','LineWidth', 2)
%     hold on
%     plot(eig_val_cell{i}, 'r--','LineWidth', 2)
% end

% figure
% plot(eig_val_cell{1}, 'b-','LineWidth', 2)
% hold on
% plot(eig_val_cell{end}, 'r--','LineWidth', 2)

% 画三维图
% for band_index = 1:4
%     figure
%     mesh(kx_q,ky_q,ene_q{band_index})
% end

% 画2d pcolormap

for band_index = 1:4
    figure
    pcolor(kx_q, ky_q, ene_q{band_index});
    hold on
    grid off
    axis([min(min(kx_q)) max(max(kx_q)) min(min(ky_q)) max(max(ky_q))])
    box on
    shading interp
    
    % 画第一布里渊区的边界
    plot(k_K_list(1,:), k_K_list(2,:), 'k--');
    plot([k_K_list(1,end), k_K_list(1,1)], [k_K_list(2,end), k_K_list(2,1)], 'k--');
end