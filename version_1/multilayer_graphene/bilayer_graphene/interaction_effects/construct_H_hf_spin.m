function H_hf_spin = construct_H_hf_spin(density_matrix_spin, S_tensor, E_ee, E_zeeman)
    % 尽管有相互作用，哈密顿量仍然是关于自旋对角的，因此我们可以分为自旋向上的块和自旋向下的块
    % H_hf = kron(diag([1, -1]), E_zeeman * diag([1, 1, 1])); % 第一项是近似简并的LL能量,第二项是Zeeman term
    H_hf_spin = E_zeeman * diag([1, 1, 1]);
    for alpha = 1:3
        for beta = 1:3
            val_spin = 0.0;
            
            for lambda = 1:3
                for sigma = 1:3
                    val_spin = val_spin + S_tensor(alpha, beta, lambda, sigma) * density_matrix_spin(sigma, lambda);
                end
            end

            H_hf_spin(alpha, beta) = H_hf_spin(alpha, beta) - E_ee * val_spin;
        end
    end
    
end