a = 0.246; % nm
block_length = sqrt(3) * a;

onsite_params = [0.0, 0.0, 0.0, 0.0];

gamma0 = 3.165;
gamma1 = 0.391;
gamma3 = 0.31515;
gamma4 = 0.04424;
gamma_params = [gamma0, gamma1, gamma3, gamma4];

B_field = 150;
% 加入磁场之后，我们可以找出一个面积S, 其对应的磁通\Phi = BS与量子磁通\Phi_0 = BS的比值\Phi / \Phi_0
length_ratio = (a / 25.66)^2;

y_shifts = [0.0, 1/3];
BS_BS0 = 2 * pi * length_ratio * B_field * sqrt(3) / (8 * pi);

N_rows = 220; % 有11行共有5.5个原胞的厚度
dims = 2 * N_rows;

[H0, H1, H_1, y_coords_layer1, y_coords_layer2] = construct_bilayer_GNR_zigzag(gamma_params, onsite_params, N_rows, y_shifts, BS_BS0);
% y_coords = [y_coords_layer1; y_coords_layer2];
y_coords = [y_coords_layer1; y_coords_layer1];
Im = 1j;

% 生成k格点(一维)
k_points = 1000;
% aks = linspace(0, 2*pi, k_points);
aks = linspace(-pi, pi, k_points);

% order_flag = false;
order_flag = false;
[eig_enes, guide_centers, vel_array, vel_dirs] = ribbon_band_solver(H0, H1, H_1, y_coords, aks, eps, order_flag);

% eig_enes = zeros(dims, k_points);
% 
% guide_centers = zeros(dims, k_points);
% vel_dirs = zeros(dims, k_points);
% vel_array = zeros(dims, k_points);
% 
% tic
% idx = 1;
% Hk = H0 + H1 * exp(Im * aks(idx)) + H_1 * exp(-Im * aks(idx));
% [eigvecs_last, eigvals_last] = eig(Hk);
% eigvals_last = diag(eigvals_last);
% [~, degeneracy, eig_num, first_idxs, last_idxs] = helper_get_degeneracy(eigvals_last, dims);
% 
% Vel_k = Im * ( H1 * exp(Im * aks(idx)) - H_1 * exp(-Im * aks(idx)) );
% [eigvecs_last, vel_dirs(:, idx), vel_array(:, idx)] = get_bloch_velocity(Vel_k, eigvecs_last, degeneracy, eig_num, first_idxs, last_idxs, dims);
% 
% eig_enes(:,idx) = eigvals_last;
% 
% for num = 1:dims
%     guide_center = 0;
%     guide_centers(num, idx) = dot(y_coords, abs(eigvecs_last(:, num)).^2);
% end
% 
% eps = 0.2 * gamma0;
% for idx = 2:k_points
%     Hk = H0 + H1 * exp(Im * aks(idx)) + H_1 * exp(-Im * aks(idx));
%     [eigvecs_now, eigvals_now] = eig(Hk);
%     eigvals_now = diag(eigvals_now);
%     [~, degeneracy, eig_num, first_idxs, last_idxs] = helper_get_degeneracy(eigvals_now, dims); % 得到简并度
%     
%     % bloch velocity 
%     Vel_k = Im * (H1 * exp(Im * aks(idx)) - H_1 * exp(-Im * aks(idx)));
%     [eigvecs_now, vel_dirs(:, idx), vel_array(:, idx)] = get_bloch_velocity(Vel_k, eigvecs_now, degeneracy, eig_num, first_idxs, last_idxs, dims);
%     
% %     delta_k = aks(idx) - aks(idx - 1);
% %     eigvals_predict = eigvals_now - delta_k * vel_array(:, idx);
% %     new_order = get_new_order(eigvecs_last, eigvals_last, eigvecs_now, eigvals_predict, degeneracy, first_idxs, last_idxs, dims, eig_num, eps);
% %     
% %     % 交换位置
% %     eigvecs_now(:,:) = eigvecs_now(:,new_order);
% %     eigvals_now(:) = eigvals_now(new_order);
% %     vel_dirs(:, idx) = vel_dirs(new_order, idx);
% %     vel_array(:, idx) = vel_array(new_order, idx);
%     
%     % 计算guide center
%     for num = 1:dims
%         guide_center = 0;
%         guide_centers(num, idx) = dot(y_coords, abs(eigvecs_now(:, num)).^2);
%     end
%     
%     eig_enes(:,idx) = eigvals_now;
%     
%     eigvecs_last = eigvecs_now;
%     eigvals_last = eigvals_now;
% end
% toc

figure
axis([aks(1) aks(end) -0.6 0.6])
hold on
for i = 1:dims
plot(aks, eig_enes(i,:))
end
hold on

figure
axis([aks(1) aks(end) -0.6 0.6])
sizeMarker = linspace(y_coords(end), y_coords(1), dims);    % 比0大，值越大标记越大
colorMarker = guide_centers;   % 颜色渐变
for i = 1:dims
%     plot(aks, eig_enes(i,:))
    patch([aks NaN],[eig_enes(i,:) NaN],[guide_centers(i,:) NaN],'Marker','o','MarkerSize',4,'EdgeColor','interp','MarkerFaceColor','flat')
    hold on
end

function new_order = get_new_order(eigvecs_last, eigvals_last, eigvecs_now, eigvals_now, degeneracy, first_idxs, last_idxs, dims, eig_num, eps)
    new_order = zeros(dims, 1);
    is_find = zeros(1,dims);
    for num = 1:eig_num
        if degeneracy(num) == 1  % 简并度为1的态
            first = first_idxs(num);
            vec_now = eigvecs_now(:,first);
            val_now = eigvals_now(first);
            % 开始地毯式搜索
            target = 0;
            max_id = 0.0;
            for jj = 1:dims
                if is_find(jj) == 0 && abs(val_now - eigvals_last(jj)) < eps
                    inner_dot = abs(dot(vec_now, eigvecs_last(:,jj)));
                    if inner_dot > max_id  % 理想上应该是逼近于1的
                        max_id = inner_dot;
                        target = jj;
                    end
                end
            end
            new_order(target) = first;  % 注意不是new_order(first) = target !!!
            is_find(target) = 1;
            
        else  % 多重简并的态
            first = first_idxs(num);
            last = last_idxs(num);
            
            inner_dots = zeros(dims,1);
            for jj = 1:dims
                if is_find(jj) == 0
                    inner_dot = 0.0;
                    for kk = first:last
                        inner_dot = inner_dot + abs(dot(eigvecs_last(:,jj), eigvecs_now(:,kk)));
                    end
                    inner_dots(jj) = inner_dot;
                else
                    inner_dots(jj) = 0.0;
                end
            end
            
            [~, idxs] = sort(inner_dots);
            new_order(idxs(end-degeneracy(num)+1:end)) = first:last;
            is_find(idxs(end-degeneracy(num)+1:end)) = 1;
        end
    end
end