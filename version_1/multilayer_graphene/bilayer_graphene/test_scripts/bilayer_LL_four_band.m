%% 添加路径
addpath(".\bilayer_hamitonians\")

%% test
v = 10^6; % 来自gamma0，单位是m/s
gamma1 = 0.40;
% gamma3 = 0.9; % 0.0 / 0.1 * gamma0
v3 = 0.1 * v; % 0.0 * v / 0.1 * v
v4 = 0;
delta = 0;
u = 0.04; % 0.0 / 0.04 / 0.065 / 0.07 / 0.08
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位

B_start = 0.0;
B_end = 6;
B_steps = 601;
B_fields_list = linspace(B_start, B_end, B_steps);

LL_index_cutoff = 80;
dims = 4 * LL_index_cutoff;
eigvals_LL_K = zeros(B_steps, dims);
eigvals_LL_Kp = zeros(B_steps, dims);

tic
for B_index = 1:B_steps
    %% parameters set up
    B_field = B_fields_list(B_index);
    mag_length = 25.66 / sqrt(B_field); % 以nm为单位
    x0 = h_bar * v / mag_length * 10^9;
    x3 = h_bar * v3 / mag_length * 10^9;
    x4 = h_bar * v4 / mag_length * 10^9;
    
    %% construct Hamiltonian @ valley K
    % [Ham_LL_Kp, Ham_LL_K] = construct_bilayer_LL_four_bands_test(x0, gamma1, x3, u, LL_index_cutoff, dims);
    Ham_LL_K = construct_bilayer_LL_four_bands(x0, x3, x4, gamma1, delta, u, +1, LL_index_cutoff, dims);
    Ham_LL_Kp = construct_bilayer_LL_four_bands(x0, x3, x4, gamma1, delta, u, -1, LL_index_cutoff, dims);
    
    % helper_check_hermite(Ham_LL_K, 1e-8);
    
    %% 对哈密顿量进行对角化
    % call the eig sovler
    [eigvec_HK_now, eigval_HK] = eig(Ham_LL_K);
    eigval_HK_diag_now = diag(eigval_HK);
    
    [eigvec_HKp_now, eigval_HKp] = eig(Ham_LL_Kp);
    eigval_HKp_diag_now = diag(eigval_HKp);
    
    % [eigvec_HK_now, eigval_HK_diag_now, eig_num_K] = helper_re_order_states(eigvec_HK_last, eigval_HK_diag_last, ...
    %                                                                         eigvec_HK_now, eigval_HK_diag_now, dims, eps);

    % push into the LLs
    eigvals_LL_K(B_index, :) = eigval_HK_diag_now;
    eigvals_LL_Kp(B_index, :) = eigval_HKp_diag_now;
    
end
toc

%% 进行作图
addpath("..\plot_funcs\")
eigvals_LL_cell = cell(2, 4);
eigvals_LL_cell{1,1} = eigvals_LL_K / gamma1;
eigvals_LL_cell{1,2} = dims; % 维数
eigvals_LL_cell{1,3} = 'b'; % line color
eigvals_LL_cell{1,4} = 0.5; % line width

eigvals_LL_cell{2,1} = eigvals_LL_Kp / gamma1;
eigvals_LL_cell{2,2} = dims;
eigvals_LL_cell{2,3} = 'b--';
eigvals_LL_cell{2,4} = 1.0;

save_path = "";
fig0 = plot_LLs(B_start, B_end, -0.125, 0.125, B_fields_list, eigvals_LL_cell, save_path);

% figure
% % subplot(2,1,1)
% axis([B_start B_end -0.035 0.035])
% hold on
% 
% % plot bilayer-like LL
% for i = 1:dims
%     plot(B_fields_list, eigvals_LL_K(:,i), 'b', 'LineWidth', 0.5)
%     plot(B_fields_list, eigvals_LL_Kp(:,i), 'b--', 'LineWidth', 1.0)
% end

%% 函数
function H_n = get_diagonal_H_n(x, gamma1, u, valley, n)
    if n == 1
        % diagonal elements
        H_n = valley * u / 2 * diag([1, 1, -1, 1]);
        % off-diagonal elements
        H_n(2, 4) = 1j * x * sqrt(2);
        H_n(4, 2) = - 1j * x * sqrt(2);
        H_n(3, 4) = valley * gamma1;
        H_n(4, 3) = valley * gamma1;
    else
        % diagonal elements
        H_n = valley * u / 2 * diag([1, -1, -1, 1]);
        % off-diagonal elements
        H_n(1, 4) = 1j * x * sqrt(2 * n);
        H_n(4, 1) = - 1j * x * sqrt(2 * n);
        H_n(2, 3) = - 1j * x * sqrt(2 * (n - 1));
        H_n(3, 2) = 1j * x * sqrt(2 * (n - 1));
        H_n(3, 4) = valley * gamma1;
        H_n(4, 3) = valley * gamma1;
    end
end

function W_n = get_off_diagonal_W_n(x3, n)
    W_n = zeros(4);
    if n == 2
        W_n(2, 2) = - 1j * x3 * 2;
    else
        W_n(1, 2) = - 1j * x3 * sqrt(2 * n);
    end
end

function Ham_LL_valley = get_bilayer_LL_four_bands(x, gamma1, x3, u, valley, LL_index_cutoff, dims)
    % 构建在谷valley处的LL能级哈密顿量
    Ham_LL_valley = zeros(dims);
    
    % 对 n = 1 block进行单独处理
    n = 1; 
    index_start = 4 * (n - 1) + 1;
    index_end = 4 * n;
    
    % H_n
    H_n = get_diagonal_H_n(x, gamma1, u, valley, n);
    Ham_LL_valley(index_start:index_end, index_start:index_end) = H_n;
    
    % W1
    W_n = get_off_diagonal_W_n(x3, 1);
    Ham_LL_valley(index_start:index_end, (index_start + 8):(index_end + 8)) = W_n;
    Ham_LL_valley((index_start + 8):(index_end + 8), index_start:index_end) = W_n';

    % W2
    W_n = get_off_diagonal_W_n(x3, 2);
    Ham_LL_valley(index_start:index_end, (index_start + 12):(index_end + 12)) = W_n;
    Ham_LL_valley((index_start + 12):(index_end + 12), index_start:index_end) = W_n';
    
    % 对 n = 2 ： (LL_index_cutoff - 3) block进行单独处理
    for n = 2:(LL_index_cutoff - 3)
        index_start = 4 * (n - 1) + 1;
        index_end = 4 * n;
        
        % H_n
        H_n = get_diagonal_H_n(x, gamma1, u, valley, n);
        Ham_LL_valley(index_start:index_end, index_start:index_end) = H_n;

        % Wn (n >= 3)
        W_n = get_off_diagonal_W_n(x3, n + 1);
        Ham_LL_valley(index_start:index_end, (index_start + 12):(index_end + 12)) = W_n;
        Ham_LL_valley((index_start + 12):(index_end + 12), index_start:index_end) = W_n';
    end
    
    % 对 n = (LL_index_cutoff - 2) : LL_index_cutoff block 进行单独处理
    for n = (LL_index_cutoff - 2):LL_index_cutoff
        index_start = 4 * (n - 1) + 1;
        index_end = 4 * n;

        % H_n
        H_n = get_diagonal_H_n(x, gamma1, u, valley, n);
        Ham_LL_valley(index_start:index_end, index_start:index_end) = H_n;
        Ham_LL_valley(index_start:index_end, index_start:index_end) = H_n;
    end
    
end

function [Ham_LL_K, Ham_LL_Kp] = construct_bilayer_LL_four_bands_test(x, gamma1, x3, u, LL_index_cutoff, dims)
    % valley : + 1 for K valley ; -1 for Kp valley
    %% construct Hamiltonian @ valley K and valley Kp
    valley = 1; % 1 for K, -1 for Kp
    Ham_LL_K = get_bilayer_LL_four_bands(x, gamma1, x3, u, valley, LL_index_cutoff, dims);
    
    valley = - 1; % 1 for K, -1 for Kp
    Ham_LL_Kp = get_bilayer_LL_four_bands(x, gamma1, x3, u, valley, LL_index_cutoff, dims);
end