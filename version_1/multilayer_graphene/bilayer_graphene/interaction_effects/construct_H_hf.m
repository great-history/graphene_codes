%% 构造哈密顿量的函数
function [H_hf_up, H_hf_down] = construct_H_hf(density_matrix_up, density_matrix_down, S_tensor, E_ee, E_zeeman)
    % 尽管有相互作用，哈密顿量仍然是关于自旋对角的，因此我们可以分为自旋向上的块和自旋向下的块
    % H_hf = kron(diag([1, -1]), E_zeeman * diag([1, 1, 1])); % 第一项是近似简并的LL能量,第二项是Zeeman term
    H_hf_up = E_zeeman * diag([1, 1, 1]);
    H_hf_down = - E_zeeman * diag([1, 1, 1]);
    for alpha = 1:3
        for beta = 1:3
            val_up = 0.0; % 用来计算自旋向上的矩阵元
            val_down = 0.0; % 用来计算自旋向下的矩阵元
            
            for lambda = 1:3
                for sigma = 1:3
                    % 看thesis发现应该是density_matrix_up(sigma, lambda)？？？
                    % val_up = val_up + S_tensor(alpha, beta, lambda, sigma) * density_matrix_up(lambda, sigma);
                    % val_down = val_down + S_tensor(alpha, beta, lambda, sigma) * density_matrix_down(lambda, sigma);
                    val_up = val_up + S_tensor(alpha, beta, lambda, sigma) * density_matrix_up(sigma, lambda);
                    val_down = val_down + S_tensor(alpha, beta, lambda, sigma) * density_matrix_down(sigma, lambda);
                end
            end

            H_hf_up(alpha, beta) = H_hf_up(alpha, beta) - E_ee * val_up; % spin up
            H_hf_down(alpha, beta) = H_hf_down(alpha, beta) - E_ee * val_down; % spin down
        end
    end
end