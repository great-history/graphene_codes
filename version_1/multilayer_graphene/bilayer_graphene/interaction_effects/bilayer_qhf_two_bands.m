%% parameters set up
format long
v = 1 * 10^6; % 来自gamma0，单位是m/s
gamma1 = 0.40;
% gamma3 = 0.1 * gamma0; % 0.0 / 0.1 * gamma0
v3 = 0.1 * v; % 0.0 * v / 0.1 * v
v4 = 0;
delta = 0;
u = 0.2 * gamma1; % 0.0 / 0.04 / 0.065 / 0.07 / 0.08
h_bar = 4.1356676969 * 10^(-15) / (2 * pi); % 以eV*s为单位
single_electron = 1.6021892 * 10^(-19); % 以C为单位

B_field = 1.0;
mag_length = 25.66 / sqrt(B_field); % 以nm为单位
x0 = h_bar * v / mag_length * 10^9;
x3 = h_bar * v3 / mag_length * 10^9;
x4 = h_bar * v4 / mag_length * 10^9;

LL_index_cutoff = 12;
dims = 2 * LL_index_cutoff;
eigvals_LL_K = zeros(1, dims);
eigvals_LL_Kp = zeros(1, dims);

%% construct Hamiltonian @ valley K
Ham_LL_K = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, +1, LL_index_cutoff, dims);
Ham_LL_Kp = construct_bilayer_LL_two_bands(x0, x3, gamma1, u, -1, LL_index_cutoff, dims);

% helper_check_hermite(Ham_LL_K, 1e-8);

%% 对哈密顿量进行对角化
% call the eig sovler
[eigvec_HK_now, eigval_HK] = eig(Ham_LL_K);
eigval_HK_diag_now = diag(eigval_HK);

[eigvec_HKp_now, eigval_HKp] = eig(Ham_LL_Kp);
eigval_HKp_diag_now = diag(eigval_HKp);

% push into the LLs
eigvals_LL_K(:) = eigval_HK_diag_now;
eigvals_LL_Kp(:) = eigval_HKp_diag_now;

%% 选取最靠近CNP的6条朗道能级(3条位于导带，3条位于价带) ———— triplets states
% K valley 取价带能量绝对值最小的三个, Kp valley 取导带能量绝对值最小的三个
negative_index_list = find(eigvals_LL_K < 0);
negative_ene_vals = eigvals_LL_K(eigvals_LL_K < 0);
[~, sort_index_list] = sort(negative_ene_vals, 'descend');
eigvals_LLL_K = eigvals_LL_K(negative_index_list(sort_index_list(1:3)));
eigvecs_LLL_K = eigvec_HK_now(:, negative_index_list(sort_index_list(1:3)));

positive_index_list = find(eigvals_LL_Kp > 0);
positive_ene_vals = eigvals_LL_Kp(eigvals_LL_Kp > 0);
[~, sort_index_list] = sort(positive_ene_vals, 'ascend');
eigvals_LLL_Kp = eigvals_LL_Kp(positive_index_list(sort_index_list(1:3)));
eigvecs_LLL_Kp = eigvec_HKp_now(:, positive_index_list(sort_index_list(1:3)));

% check
% for i = 1:3
%     error_max = max(abs(eigvals_LLL_Kp(i) * eigvecs_LLL_Kp(:, i) - Ham_LL_Kp * eigvecs_LLL_Kp(:, i)));
%     disp(error_max)
%     error_max = max(abs(eigvals_LLL_K(i) * eigvecs_LLL_K(:, i) - Ham_LL_K * eigvecs_LLL_K(:, i)));
%     disp(error_max)
% end

%% take the interaction effect into account
% The structure of S matrix element : S^{\lambda \sigma}_{\alpha \beta}
addpath('.\interaction_effects\')
% 这里只考虑Kp的三个LL能级
eigvecs_LLL_Kp_A_component = zeros(LL_index_cutoff + 1, 3);
eigvecs_LLL_Kp_B_component = zeros(LL_index_cutoff + 1, 3);
for ii = 1:3
    [eigvec_A_component, eigvec_B_component] = get_LL_components_each_sublattice_two_bands(eigvecs_LLL_Kp(:, ii), LL_index_cutoff, -1);
    eigvecs_LLL_Kp_A_component(:, ii) = eigvec_A_component;
    eigvecs_LLL_Kp_B_component(:, ii) = eigvec_B_component;
end

% 计算(3*3)*(3*3)的S tensor
d = 0.335; % 双层的层间距离，单位为nm
d_interlayer = d / 25.6 * sqrt(B_field); % 这里计算的是第一层和第二层之间的exchange
tic
disp("开始计算S矩阵")
nm_combination_list = zeros(4, (LL_index_cutoff + 1)^4);
count = 0;
for n1 = 0:LL_index_cutoff
    for n2 = 0:LL_index_cutoff
        for m1 = 0:LL_index_cutoff
            for m2 = 0:LL_index_cutoff
                count = count + 1;
                nm_combination_list(1, count) = n1;
                nm_combination_list(2, count) = n2;
                nm_combination_list(3, count) = m1;
                nm_combination_list(4, count) = m2;
            end
        end
    end
end

% for n1 = 0:LL_index_cutoff
%     for n2 = 0:LL_index_cutoff
%         for m1 = 0:LL_index_cutoff
%             for m2 = 0:LL_index_cutoff
%             end
%         end
%     end
% end

poolobj = gcp('nocreate');
if isempty(poolobj)
    disp('启动并行运算，核心数：6');
    % Perform a basic check by entering this code, where "local" is one kind of cluster profile.
    parpool('local', 6);
else
    disp(['并行运o算已启动，核心数：' num2str(poolobj.NumWorkers)]);
end

S_tensor = zeros(3,3,3,3); % 前两个指标是(\lambda, \sigma)， 后两个指标是(\alpha, \beta)
S_tensor_temp = zeros(count, 3, 3, 3, 3);
parfor ii = 1:count
    n1 = nm_combination_list(1, ii);
    n2 = nm_combination_list(2, ii);
    m1 = nm_combination_list(3, ii);
    m2 = nm_combination_list(4, ii);
    
    X_nm_intralayer = get_exchange_integral(n1, n2, m1, m2); % intralayer
    X_nm_interlayer = get_exchange_integral(n1, n2, m1, m2, d_interlayer); % interlayer
    
    if X_nm_intralayer == 0
        continue
    end
    
    for alpha = 1:3 % 对应n1 产生
        for beta = 1:3 % 对应m2 湮灭
            for lambda = 1:3 % 对应m1 产生
                for sigma = 1:3 % 对应n2 湮灭
                    aa_alpha_sigma = conj(eigvecs_LLL_Kp_A_component(n1 + 1, alpha)) * eigvecs_LLL_Kp_A_component(n2 + 1, sigma);
                    bb_alpha_sigma = conj(eigvecs_LLL_Kp_B_component(n1 + 1, alpha)) * eigvecs_LLL_Kp_B_component(n2 + 1, sigma);
                    aa_lambda_beta = conj(eigvecs_LLL_Kp_A_component(m1 + 1, lambda)) * eigvecs_LLL_Kp_A_component(m2 + 1, beta);
                    bb_lambda_beta = conj(eigvecs_LLL_Kp_B_component(m1 + 1, lambda)) * eigvecs_LLL_Kp_B_component(m2 + 1, beta);

                    current_val = (aa_alpha_sigma * aa_lambda_beta + bb_alpha_sigma * bb_lambda_beta) * X_nm_intralayer + ...
                                  (aa_alpha_sigma * bb_lambda_beta + bb_alpha_sigma * aa_lambda_beta) * X_nm_interlayer;
                    S_tensor_temp(ii, lambda, sigma, alpha, beta) = current_val;
                end
            end
        end
    end
end


for alpha = 1:3 % 对应n1 产生
    for beta = 1:3 % 对应m2 湮灭
        for lambda = 1:3 % 对应m1 产生
            for sigma = 1:3 % 对应n2 湮灭
                
                val = 0.0;
                parfor ii = 1:count
                    val = val + S_tensor_temp(ii, lambda, sigma, alpha, beta);
                end
                
                S_tensor(lambda, sigma, alpha, beta) = val;
            end
        end
    end
end

disp("S矩阵计算结束")
save('S_tensor_B_1T.mat','S_tensor');
save('eigvals_LLL_K_B_1T.mat', 'eigvals_LLL_K')
save('eigvals_LLL_Kp_B_1T.mat', 'eigvals_LLL_Kp')
toc



%% test
% A = rand(2000,2000);
% tic
% A1 = gpuArray(single(A));
% [U,S,V] = svd(A1,'econ');
% A2 = U*S*V';
% A3 = gather(A2);
% toc
% tic
% [U1,S1,V1] = svd(A,'econ');
% A4 = U1*S1*V1';
% toc
% error = norm(A3-A4,'fro');
% tic
% A1 = gpuArray(single(A));
% [vv,dd] = eig(A1);
% toc

% M = rand(2000,2000);            % 生成一个随机矩阵
% tic
% [A1,B1] = eig(M);                    % 求该随机矩阵的特征值和特征向量
% toc
% tic
% M = single(M);                     % 将数据转换为单精度型
% M = gpuArray(M);                % 将数据从CPU中搬到GPU
% [A2,B2] = eig(M);                 % 求特征值和特征向量
% A2 = gather(A2);                 % 将数据从GPU中搬到CPU
% toc

%% 得到sublattice上朗道能级对应的指标
% function [eigvec_A_component, eigvec_B_component] = get_LL_components_each_sublattice_two_bands(eigvec, LL_index_cutoff, valley)
%     % two bands model : (B2, A1) @ K valley & (A1, B2) @ Kp valley
%     eigvec_A_component = zeros(LL_index_cutoff + 1, 1);
%     eigvec_B_component = zeros(LL_index_cutoff + 1, 1);
%     
%     if valley == 1 % for K valley
%         eigvec_B_component(1) = eigvec(1); % |B2, 0>的分量
%         eigvec_B_component(2) = eigvec(2); % |B2, 1>的分量
%         for n = 2:LL_index_cutoff
%             eigvec_B_component(n + 1) = eigvec(2*(n - 1) + 1); % |B2, n>的分量
%             eigvec_A_component(n - 1) = eigvec(2*n); % |A1, n - 2>的分量
%         end
%    
%     else % for Kp valley
%         eigvec_A_component(1) = eigvec(1); % |A1, 0>的分量
%         eigvec_A_component(2) = eigvec(2); % |A1, 1>的分量
%         for n = 2:LL_index_cutoff
%             eigvec_A_component(n + 1) = eigvec(2*(n - 1) + 1); % |A1, n>的分量
%             eigvec_B_component(n - 1) = eigvec(2*n); % |B2, n - 2>的分量
%         end
%         
%     end
% end