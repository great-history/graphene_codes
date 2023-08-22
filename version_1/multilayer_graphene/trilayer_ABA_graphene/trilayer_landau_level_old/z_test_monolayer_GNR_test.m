%% add magnetic field on zigzag
% a = 0.246; % nm
% block_length = sqrt(3) * a;
% gamma0 = 3.1; % eV
% 
% t = 1;
% B_field = 250;
% 
% % 加入磁场之后，我们可以找出一个面积S, 其对应的磁通\Phi = BS与量子磁通\Phi_0 = BS的比值\Phi / \Phi_0
% length_ratio = (a / 25.66)^2;
% % BS_BS0 = 2 * pi * length_ratio * B_field * sqrt(3) / (8 * pi);  % B_field = 10T时，如果宽度只有100个原胞的话，朗道能级是看不出来的
% % 如果要看到朗道能级，要么就是宽度要足够大，要么就是磁场足够高
% 
% N_rows = 220;  % 有11行共有5.5个原胞的厚度
% dims = N_rows;
% N_block = floor(N_rows / 4);
% half_block = mod(N_rows, 4);
% % 生成y坐标(沿着ribbon宽度方向)
% y_coords = zeros(dims, 1);
% % for n = 1:N_block
% %     y_coords(4*(n-1) + 1) = (2 * n - 13 / 6);
% %     y_coords(4*(n-1) + 2) = (2 * n - 11 / 6);
% %     y_coords(4*(n-1) + 3) = (2 * n - 7 / 6);
% %     y_coords(4*(n-1) + 4) = (2 * n - 5 / 6);
% % end
% y_shift = 0;
% for n = 1:N_block
%     y_coords(4*(n-1) + 1) = (2 * n - 5 / 3) + y_shift;
%     y_coords(4*(n-1) + 2) = (2 * n - 4 / 3) + y_shift;
%     y_coords(4*(n-1) + 3) = (2 * n - 2 / 3) + y_shift;
%     y_coords(4*(n-1) + 4) = (2 * n - 1 / 3) + y_shift;
% end
% 
% BS_BS0 = 2 * pi * length_ratio * B_field * sqrt(3) / (8 * pi);
% % BS_BS0 = 1 / 90;
% 
% % construct the tight-binding model under a magnetic field
% H0 = zeros(dims, dims);
% H1 = zeros(dims, dims);
% mass = 0;
% 
% Im = 1j;
% for n = 1:N_block
% %     phi1 = BS_BS0 * (2 * n - 1 / 2);  % 12
% %     phi2 = BS_BS0 * (2 * n - 3 / 2);  % 34 或 43
%     phi1 = BS_BS0 * (y_coords(4*(n-1) + 4) + y_coords(4*(n-1) + 3)) / 2;
%     phi2 = BS_BS0 * (y_coords(4*(n-1) + 2) + y_coords(4*(n-1) + 1)) / 2;
%     
%     H0(4*(n-1)+4, 4*(n-1)+4) = mass;
%     H0(4*(n-1)+3, 4*(n-1)+3) = -mass;
%     H0(4*(n-1)+2, 4*(n-1)+2) = mass;
%     H0(4*(n-1)+1, 4*(n-1)+1) = -mass;
%     
%     H0(4*(n-1)+3, 4*(n-1)+4) = t * exp(-Im * phi1);
%     H0(4*(n-1)+4, 4*(n-1)+3) = conj(H0(4*(n-1)+3, 4*(n-1)+4));
%     
%     H0(4*(n-1)+2, 4*(n-1)+3) = t;
%     H0(4*(n-1)+3, 4*(n-1)+2) = conj(H0(4*(n-1)+2, 4*(n-1)+3));
%     
%     H0(4*(n-1)+1, 4*(n-1)+2) = t * exp(Im * phi2);
%     H0(4*(n-1)+2, 4*(n-1)+1) = conj(H0(4*(n-1)+1, 4*(n-1)+2));
%     
%     if ~(n == N_block)
%         H0(4*n + 1, 4*(n-1) + 4) = t;
%         H0(4*(n-1) + 4, 4*n + 1) = conj(H0(4*n + 1, 4*(n-1)+4));
%     end
%     
%     % H1不是厄密的矩阵
%     H1(4*(n-1)+4, 4*(n-1)+3) = t * exp(-Im * phi1); 
%     H1(4*(n-1)+1, 4*(n-1)+2) = t * exp(-Im * phi2);
%     
%     % H_{-1}是H1的厄密，故不用创建
% end
% H_1 = H1'; % H_{-1}是H1的厄密
% --------------------------------------------------------------------------------------------------------------------------
a = 0.246; % nm
block_length = sqrt(3) * a;
gamma0 = 3.1; % eV
mass = 0;

t = 1;
B_field = 150;

% 加入磁场之后，我们可以找出一个面积S, 其对应的磁通\Phi = BS与量子磁通\Phi_0 = BS的比值\Phi / \Phi_0
length_ratio = (a / 25.66)^2;
N_rows = 220; % 有11行共有5.5个原胞的厚度
dims = N_rows;
Im = 1j;
y_shift = 0;
BS_BS0 = 2 * pi * length_ratio * B_field * sqrt(3) / (8 * pi);

[H0, H1, H_1, y_coords] = construct_monolayer_GNR_zigzag(t, [0.0,0.0], N_rows, y_shift, BS_BS0);
% hem = helper_check_hermite(H0, 1e-8);


% 生成k格点(一维)
k_points = 1000;
aks = linspace(0, 2*pi, k_points);
% aks = linspace(-pi, pi, k_points);

eps = 0.1 * gamma0;

% order_flag = false;
order_flag = true;
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
% eps = 0.1 * t;
% % for idx = 2:k_points
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
%     % new_order = helper_get_new_order(eigvecs_last, eigvals_last, eigvecs_now, eigvals_now, degeneracy, first_idxs, last_idxs, dims, eig_num, eps);
%     delta_k = aks(idx) - aks(idx - 1);
%     eigvals_predict = eigvals_now - delta_k * vel_array(:, idx);
%     
%     new_order = helper_get_new_order(eigvecs_last, eigvals_last, eigvecs_now, eigvals_predict, degeneracy, first_idxs, last_idxs, dims, eig_num, eps);
%     % 交换位置
%     eigvecs_now(:,:) = eigvecs_now(:,new_order);
%     eigvals_now(:) = eigvals_now(new_order);
%     vel_dirs(:, idx) = vel_dirs(new_order, idx);
%     vel_array(:, idx) = vel_array(new_order, idx);
%     eig_enes(:,idx) = eigvals_now;
%     
%     % 计算guide center
%     for num = 1:dims
%         guide_center = 0;
%         guide_centers(num, idx) = dot(y_coords, abs(eigvecs_now(:, num)).^2);
%     end
%     
%     eigvecs_last = eigvecs_now;
%     eigvals_last = eigvals_now;
% end
% toc

% % % perform a shift along the kx:通过找零朗道能级的中心
% % k_start_index = 410;
% % k_end_index = 1000;
% % k_mid_index = ceil((k_start_index + k_end_index) / 2);
% % k_mid_index_ori = ceil(k_points / 2);
% % 
% % if k_mid_index > k_mid_index_ori
% %     delta_k_index = k_mid_index - k_mid_index_ori;
% %     select = [(delta_k_index+1):k_points, 1:delta_k_index];
% %     % re-order k points:记得要减2*pi
% %     aks1 = aks(1:delta_k_index) + 2 * pi;
% %     aks = [aks((delta_k_index+1):k_points), aks1];
% % else
% %     delta_k_index = k_mid_index_ori - k_mid_index;
% %     select = [(k_points - delta_k_index + 1):k_points, 1:(k_points - delta_k_index)];
% %     % re-order k points:记得要加2*pi
% %     aks1 = aks((k_points - delta_k_index + 1):k_points) - 2 * pi;
% %     aks = [aks1, aks(1:(k_points - delta_k_index))];
% % end
% 
% % re-order
% % eig_enes = eig_enes(:,select);
% % vel_array = vel_array(:, select);
% % vel_dirs = vel_dirs(:, select);
% % guide_centers = guide_centers(:, select);

figure
axis([aks(1) aks(end) -0.6 0.6])
hold on
for i = 1:dims
plot(aks, eig_enes(i,:))
hold on
end

figure
axis([aks(1) aks(end) -0.6 0.6])
hold on
sizeMarker = linspace(y_coords(end), y_coords(1), dims);    % 比0大，值越大标记越大
colorMarker = guide_centers;   % 颜色渐变
for i = 1:dims
%     plot(aks, eig_enes(i,:))
    patch([aks NaN],[eig_enes(i,:) NaN],[guide_centers(i,:) NaN],'Marker','o','MarkerSize',4,'EdgeColor','interp','MarkerFaceColor','flat')
    hold on
end

% 
% 
% subplot(1,2,1)
% scatter(x, y, sizeMarker, colorMarker, 'o', 'filled')
% subplot(1,2,2)
% patch([x NaN],[y NaN],[colorMarker NaN],'Marker','o','EdgeColor','interp','MarkerFaceColor','flat')

% new_order = zeros(dims, 1);
% is_find = zeros(1,dims);
% delta_k = aks(2) - aks(1);
% inner_dots = zeros(dims, 1);
% for num = 1:eig_num
%     if degeneracy(num) == 1  % 简并度为1的态
%         first = first_idxs(num);
%         vec_now = eigvecs_now(:,first);
%         val_now = eigvals_now(first);
%         val_predict = val_now - delta_k * vel_array(num,2);
%         % 开始地毯式搜索
%         target = 0;
%         max_id = 0.0;
%         for jj = 1:dims
%             if is_find(jj) == 0
%                 inner_dot = abs(dot(vec_now, eigvecs_last(:,jj)));
%                 if inner_dot > max_id  % 理想上应该是逼近于1的
%                     max_id = inner_dot;
%                     target = jj;
%                 end
%             end
%         end
%         target
%         new_order(target) = first;  % 注意不是new_order(first) = target !!!
%         is_find(target) = 1;
%     end
% end