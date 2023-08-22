a = 0.246; % nm
block_length = sqrt(3) * a;

gamma0 = 3.1;
gamma1 = 0.390;
gamma2 = -0.028;
gamma3 = 0.315;  
gamma4 = 0.041;
gamma5 = 0.050;

delta = 0.046;
Delta1 = 0.01;
Delta2 = 0.0018;

% gamma0 = 3100;
% gamma1 = 390;
% gamma2 = -28;
% gamma3 = 315;  
% gamma4 = 41;
% gamma5 = 50;
% 
% delta = 46;
% Delta1 = 0;
% Delta2 = 1.8;

gamma_params = [gamma0, gamma1, gamma2, gamma3, gamma4, gamma5];
diagonal_terms = [Delta1, Delta2, delta];
% diagonal_terms = [0.0, 0.0, 0.0, 0.0];

B_field = 150;
% 加入磁场之后，我们可以找出一个面积S, 其对应的磁通\Phi = BS与量子磁通\Phi_0 = BS的比值\Phi / \Phi_0
length_ratio = (a / 25.66)^2;
BS_BS0 = 2 * pi * length_ratio * B_field * sqrt(3) / (8 * pi);

N_rows = 220; % 有11行共有5.5个原胞的厚度
dims = 3 * N_rows;

type = "symmetry";
y_shifts = [0.0, 1 / 3, 0.0];
[H0, H1, H_1, y_coords_layer1, y_coords_layer2, y_coords_layer3] = construct_trilayer_GNR_zigzag(type, gamma_params, diagonal_terms, N_rows, y_shifts, BS_BS0);
helper_check_hermite(H0, 1e-8)
Im = 1j;

% 生成k格点(一维)
k_points = 500;
aks = linspace(0, 2*pi, k_points);
% % aks = linspace(-pi, pi, k_points);
if Delta1 == 0
    H0_m = H0(1:N_rows, 1:N_rows);
    H0_b = H0((N_rows + 1):3*N_rows, (N_rows + 1):3*N_rows);
    H1_m = H1(1:N_rows, 1:N_rows);
    H_1_m = H1_m';
    H1_b = H1((N_rows + 1):3*N_rows, (N_rows + 1):3*N_rows);
    H_1_b = H1_b';
    
    tic
    order_flag = true;
    [eig_enes_m, guide_centers_m, vel_array_m, vel_dirs_m] = ribbon_band_solver(H0_m, H1_m, H_1_m, y_coords, aks, eps, order_flag);
    order_flag = true;
    [eig_enes_b, guide_centers_b, vel_array_b, vel_dirs_b] = ribbon_band_solver(H0_b, H1_b, H_1_b, y_coords, aks, eps, order_flag);
    toc
else
    tic
    order_flag = true;
    [eig_enes, guide_centers, vel_array, vel_dirs] = ribbon_band_solver(H0, H1, H_1, y_coords, aks, eps, order_flag);
    toc
    
    figure
    axis([aks(1) aks(end) -0.6 0.6])
    hold on
    for i = 1:dims
    plot(aks, eig_enes(i,:))
    hold on
    end
end


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
% eps = 0.1 * gamma0;
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
%     eig_enes(:,idx) = eigvals_now;
%     
%     eigvecs_last = eigvecs_now;
%     eigvals_last = eigvals_now;
% end
% toc
