function H_hf_spin = construct_ABA_trilayer_H_hf_T2_spin(density_matrix_spin_temp, S_tensor, R_mat, ...
                                                         T2_eigvals_list, E_zeeman, E_H, Delta_mid, E_exchange)
    % S_tensor是用来计算fock term，而R_mat是用来计算hartree term
    % energy scale : T2_eigvals_list // E_zeeman // E_H * Delta_mid / 2 // E_exchange
    
    %% 构造在T2 basis下的Hartree-Fock Hamiltonian
    % eigenvalues
    H_hf_spin = diag(T2_eigvals_list);
    
    % Zeeman term
    H_hf_spin = H_hf_spin - E_zeeman * diag([1,1,1]);
    
    % Hartree-term
    H_hf_spin = H_hf_spin + E_H / 2 * Delta_mid * R_mat;
    
    % fock term
    fock_mat_spin = zeros(3);
    for a = 1:3
        for b = 1:3
            val_spin = 0.0;
            for c = 1:3
                for d = 1:3
                    val_spin = val_spin + density_matrix_spin_temp(c,d) * S_tensor(a,d,c,b);
                end
            end
            fock_mat_spin(a,b) = val_spin;
        end
    end
    fock_mat_spin = - E_exchange * fock_mat_spin;
    
    H_hf_spin = H_hf_spin + fock_mat_spin;
end