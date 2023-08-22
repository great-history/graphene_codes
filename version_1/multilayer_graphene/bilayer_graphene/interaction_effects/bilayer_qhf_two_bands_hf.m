% load('S_tensor_B_1T.mat');
% load('./bilayer_data/LLL_Kp_info_B_1T.mat')
% 
% eigvals_LLL_Kp =LLL_Kp_info_cell{1, 1};
% eigvecs_LLL_Kp = LLL_Kp_info_cell{1, 2};
% eigvecs_LLL_Kp_A_component = LLL_Kp_info_cell{1, 3};
% eigvecs_LLL_Kp_B_component = LLL_Kp_info_cell{1, 4};
% eigvals_LLL_Kp_quasi = LLL_Kp_info_cell{1, 5};
% eigvecs_LLL_Kp_triplet = LLL_Kp_info_cell{1, 6};
% eigvecs_LLL_Kp_A_component_triplet = LLL_Kp_info_cell{1, 7};
% eigvecs_LLL_Kp_B_component_triplet = LLL_Kp_info_cell{1, 8};
% 
% %% energy scale
% B_field = 1;
% gamma1 = 0.40;
% % Zeeman energy
% g_factor = 2;
% mu_b = 5.788 * 10^(-5); % 单位是eV / T
% E_zeeman = g_factor * mu_b * B_field;
% % e-e interaction energy scale
% E_ee = 0.14 * sqrt(B_field) * gamma1;
% filling_factor = 1;
% 
% %% construct the full Hartree-Fock Hamiltonian
% dims = 2 * 3; % spin * triplet
% % initialize the expectation of density matrix
% density_matrix = zeros(dims); % 基组为(|1\up>, |2\up>, |3\up>, |1\down>, |2\down>, |3\down>)
% % initial guess
% density_matrix(4,4) = 1;
% density_matrix(5,5) = 1;
% density_matrix(6,6) = 1;
% 
% % construct the initial full Hartree-Fock Hamiltonian
% H_hf = construct_H_hf(density_matrix, S_tensor, eigvals_LLL_Kp, E_ee, E_zeeman);
% [error, density_matrix] = loss_function(H_hf, density_matrix, filling_factor, dims);
% 
% %% start self-consistent calculation
% tic
% while error >= 1e-8
%     H_hf = construct_H_hf(density_matrix, S_tensor, eigvals_LLL_Kp, E_ee, E_zeeman);
%     [error, density_matrix] = loss_function(H_hf, density_matrix, filling_factor, dims);
% end
% toc
% 
% %% 损失函数
% function [error, density_matrix_temp] = loss_function(H_hf, density_matrix_last, filling_factor, dims)
%     [eigvecs, ~] = eig(H_hf); % MATLAB eig默认是从大到小进行排序
%     
%     density_matrix_temp = zeros(dims);
%     for ii = dims:-1:(dims - filling_factor + 1)
%         density_matrix_temp(ii, ii) = 1;
%     end
%     
%     density_matrix_temp = density_matrix_temp * transpose(eigvecs);
%     density_matrix_temp = conj(eigvecs) * density_matrix_temp;
%     
%     error_matrix = density_matrix_temp - density_matrix_last;
%     b = sum(sum(abs(density_matrix_last).^2));
%     a = sum(sum(abs(error_matrix).^2));
%     error = a / b;
% end
% 
% function H_hf = construct_H_hf(density_matrix, S_tensor, eigvals_LLL_K, E_ee, E_zeeman)
%     H_hf = kron(diag([1,1]), diag(eigvals_LLL_K)) + kron(diag([1, -1]), E_zeeman * diag([1, 1, 1])); % 第一项是近似简并的LL能量,第二项是Zeeman term
%     for alpha = 1:3
%         for beta = 1:3
%             val_up = 0.0;
%             val_down = 0.0;
%             for lambda = 1:3
%                 for sigma = 1:3
%                     val_up = val_up + S_tensor(lambda, sigma, alpha, beta) * density_matrix(lambda, sigma);
%                     val_down = val_down + S_tensor(lambda, sigma, alpha, beta) * density_matrix(lambda + 3, sigma + 3);
%                 end
%             end
% 
%             H_hf(alpha, beta) = H_hf(alpha, beta) - E_ee * val_up; % spin up
%             H_hf(alpha + 3, beta + 3) = H_hf(alpha + 3, beta + 3) - E_ee * val_down; % spin down
%         end
%     end
% end

%% test
S0000 = 0.64;
S0011 = 0.47;
S0212 = 0.29;
S0102 = -0.29;

S_tensor = zeros(3,3,3,3);
for alpha = 1:3
    for beta = 1:3
        for lambda = 1:3
            for sigma = 1:3
                % 非零必须至少有两个指标一样
                if (alpha == beta) && (alpha == lambda) && (alpha == sigma)
                    S_tensor(alpha, beta, lambda, sigma) = S0000;
                    continue % 一定要加continue!!!
                end
                
                if (alpha == beta) && (lambda == sigma)
                    S_tensor(alpha, beta, lambda, sigma) = S0011;
                    continue % 一定要加continue!!!
                end
                
                if (alpha == sigma) && (lambda == beta)
                    S_tensor(alpha, beta, lambda, sigma) = S0011;
                    continue % 一定要加continue!!!
                end
                
                if (beta == sigma)
                    if ~(alpha == lambda) && ~(alpha == beta) && ~(lambda == beta)
                        if beta == 3
                            S_tensor(alpha, beta, lambda, sigma) = S0212;
                        else
                            S_tensor(alpha, beta, lambda, sigma) = S0102;
                        end
                    end
                    
                    continue
                end
                
                if (alpha == lambda)
                    if ~(beta == sigma) && ~(beta == alpha) && ~(sigma == alpha)
                        if alpha == 3
                            S_tensor(alpha, beta, lambda, sigma) = S0212;
                        else
                            S_tensor(alpha, beta, lambda, sigma) = S0102;
                        end
                    end
                    
                    continue
                end
            end
        end
    end
end

%% energy scale
B_field = 1;
% gamma1 = 0.40;
% Zeeman energy
% g_factor = 2;
% mu_b = 5.788 * 10^(-5); % 单位是eV / T
% E_zeeman = g_factor * mu_b * B_field;
E_zeeman = 0.289 * B_field;
% e-e interaction energy scale
% E_ee = 0.14 * sqrt(B_field) * gamma1;
E_ee = 140 * sqrt(B_field) / 6;  % 6为相对介电常数

%% construct the full Hartree-Fock Hamiltonian
addpath('.\interaction_effects\')
% dims = 2 * 3; % spin * triplet

% initialize the expectation of density matrix
density_matrix_up = zeros(3); % 基组为(|1\up>, |2\up>, |3\up>, |1\down>, |2\down>, |3\down>)
density_matrix_down = zeros(3); % 基组为(|1\up>, |2\up>, |3\up>, |1\down>, |2\down>, |3\down>)
% initial guess : 自旋向上都不占据，自旋向下都占据
density_matrix_down(1,1) = 0.6;
density_matrix_down(1,2) = 0.7;
density_matrix_down(2,1) = 0.7;
density_matrix_down(2,2) = 0.5;
density_matrix_down(2,3) = 0.523;
density_matrix_down(3,2) = 0.523;
density_matrix_down(3,3) = 1.5;
% initialize density matrix
[H_hf_up, H_hf_down] = construct_H_hf(density_matrix_up, density_matrix_down, S_tensor, E_ee, E_zeeman);
% 对角化自旋向上块（check:应该没问题）
[eigvecs_up_temp, eigvals_up_temp] = eig(H_hf_up); % MATLAB eig默认是从大到小进行排序
% 对角化自旋向下块
[eigvecs_down_temp, eigvals_down_temp] = eig(H_hf_down); % MATLAB eig默认是从大到小进行排序

filling_factor = 2; % 共有2/6个LLs被填充
density_matrix_up = zeros(3); % 基组为(|1\up>, |2\up>, |3\up>, |1\down>, |2\down>, |3\down>)
density_matrix_down = zeros(3); % 基组为(|1\up>, |2\up>, |3\up>, |1\down>, |2\down>, |3\down>)
% initial guess : 自旋向上都不占据，自旋向下都占据
density_matrix_up(1,1) = 1;
density_matrix_up(2,2) = 0.5;
density_matrix_up(3,3) = 0;
density_matrix_down(1,1) = 1;
density_matrix_down(2,1) = 0.5;
density_matrix_down(1,2) = 0.5;
density_matrix_down(2,2) = 1;
density_matrix_down(3,3) = 0;

density_matrix_up = density_matrix_up * transpose(eigvecs_up_temp);
density_matrix_up = conj(eigvecs_up_temp) * density_matrix_up;
density_matrix_down = density_matrix_down * transpose(eigvecs_down_temp);
density_matrix_down = conj(eigvecs_down_temp) * density_matrix_down;

% construct the initial full Hartree-Fock Hamiltonian
[H_hf_up, H_hf_down] = construct_H_hf(density_matrix_up, density_matrix_down, S_tensor, E_ee, E_zeeman);
[error, density_matrix_up, density_matrix_down] = Hartree_Fock_iteraction(H_hf_up, H_hf_down, density_matrix_up, density_matrix_down, filling_factor);
disp(["误差为", error])
%% start self-consistent calculation
tic
while error >= 1e-8
    [H_hf_up, H_hf_down] = construct_H_hf(density_matrix_up, density_matrix_down, S_tensor, E_ee, E_zeeman);
    [error, density_matrix_up, density_matrix_down] = Hartree_Fock_iteraction(H_hf_up, H_hf_down, density_matrix_up, density_matrix_down, filling_factor);
    % disp(["误差为", error])
end
disp(["误差为", error])
toc