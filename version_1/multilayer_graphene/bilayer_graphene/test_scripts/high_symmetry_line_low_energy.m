%% 参数设置
gamma0 = 3.16; % 完全pi键
gamma1 = - 0.381; % 完全sigma键
gamma3 = 0.38;
gamma4 = 0.14;
delta_dimer = 0.022;
U = 0.05;

cc_bond1 = [sqrt(3)/2;1/2];
cc_bond2 = [-sqrt(3)/2;1/2];
cc_bond3 = [0;-1];

%% 对布里渊区进行撒点
np_side = 201;
ak_norm = 0.5;
akx_array = ak_norm * linspace(-1,1,np_side);
aky_array = ak_norm * linspace(-1,1,np_side);

%% 计算哈密顿量
eig_val_K_matrix = zeros(np_side, np_side, 4);
eig_vec_K_cell = cell(np_side, 1);
eig_val_Kp_matrix = zeros(np_side, np_side, 4, 4);
eig_vec_Kp_cell = cell(np_side, 1);

tic
for ky_index = 1:np_side
    num_k_temp = np_side;
    aky = aky_array(ky_index);
    eig_val_K_array_temp = zeros(num_k_temp, 4);
    eig_vec_K_array_temp = zeros(num_k_temp, 4, 4);
    eig_val_Kp_array_temp = zeros(num_k_temp, 4);
    eig_vec_Kp_array_temp = zeros(num_k_temp, 4, 4);
    
    for kx_index = 1:num_k_temp
        akx = akx_array(kx_index);
        
        [HK_ham, HKp_ham] = construct_bilayer_continuum_model(gamma0, gamma1, gamma3, gamma4, delta_dimer, akx, aky);
        
        hem = helper_check_hermite(HK_ham, 1e-8);
        if hem == 0
          disp("monolayer Ham is not hermitian")
        end
        
        % 对角化
        % valley K
        [eig_vecs_K_temp, eig_vals_K_temp] = eig(HK_ham);
        eig_vals_K_temp = diag(eig_vals_K_temp);

        % 按照绝对值重新进行排序, 可以减少计算次数
        [~, new_order] = sort(eig_vals_K_temp);
        
        eig_val_K_array_temp(kx_index, :) = eig_vals_K_temp(new_order);
        eig_vec_K_array_temp(kx_index, :, :) = eig_vecs_K_temp(:, new_order);
    end
    
    eig_val_K_matrix(ky_index, :, :) = eig_val_K_array_temp;
    eig_vec_K_cell{ky_index} = eig_vec_K_array_temp;
end
toc

%% 作图
% for i = 1:(2 * np_side - 1)
%     figure
%     plot(eig_val_cell{i}, 'b-','LineWidth', 2)
% end
% 
figure
plot(akx_array, eig_val_K_matrix(:, 1, 1), 'b-','LineWidth', 2)
hold on
plot(akx_array, eig_val_K_matrix(:, end, 1), 'r--','LineWidth', 2)
xlim([min(akx_array), max(akx_array)])

% 画三维图
for band_index = 1:4
    figure
    s = surf(akx_array, aky_array, eig_val_K_matrix(:, :, band_index),'FaceAlpha',0.5);
    s.EdgeColor = 'none';
end