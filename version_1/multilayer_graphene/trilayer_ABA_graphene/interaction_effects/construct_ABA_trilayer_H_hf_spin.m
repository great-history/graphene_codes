function H_hf_spin = construct_ABA_trilayer_H_hf_spin(density_matrix_spin_temp, S_tensor, R_mat, ...
                                                      E_LL_list, E_exchange, E_H, Delta_mid, E_zeeman, spin, dim_subspace)
    % 构造ABA trilayer graphene的full Hartree-Fock Hamiltonian in any basis
    % energy scale : T2_eigvals_list // E_zeeman // E_H * Delta_mid / 2 // E_exchange
    
    %% 构造Hartree-Fock Hamiltonian
    % eigenvalues
    H_hf_spin = diag(E_LL_list);
    
    % Zeeman term
    H_hf_spin = H_hf_spin + spin * E_zeeman * diag(ones(1,dim_subspace));
    
    % Hartree-term
    H_hf_spin = H_hf_spin + E_H / 2 * Delta_mid * R_mat;
    
    % fock term
    fock_mat_spin = zeros(dim_subspace);
    for a = 1:dim_subspace
        for b = 1:dim_subspace
            val_spin = 0.0;
            for c = 1:dim_subspace
                for d = 1:dim_subspace
                    val_spin = val_spin + density_matrix_spin_temp(c,d) * S_tensor(a,d,c,b);
                end
            end
            fock_mat_spin(a,b) = val_spin;
        end
    end
    fock_mat_spin = - E_exchange * fock_mat_spin;
    
    H_hf_spin = H_hf_spin + fock_mat_spin;
    
end
