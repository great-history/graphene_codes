function H_ham = construct_bilayer_tight_binding_model(gamma0, gamma1, gamma3, gamma4, delta_dimer, akx, aky)
    %% 构造哈密顿量（紧束缚模型）
    % 不包括除了delta之外的onsite potential
    % basis : A1, B1, A2, B2
    cc_bond1 = [sqrt(3)/2;1/2];
    cc_bond2 = [-sqrt(3)/2;1/2];
    cc_bond3 = [0;-1];
    ak_vec = [akx; aky];
    fk = exp(-1j * dot(ak_vec, cc_bond1)) + exp(-1j * dot(ak_vec, cc_bond2)) + exp(-1j * dot(ak_vec, cc_bond3));
    
    H_ham = zeros(4,4);
    H_ham(2,2) = delta_dimer;
    H_ham(3,3) = delta_dimer;
    
    h0 = gamma0 * fk;
    h3 = gamma3 * fk;
    h4 = gamma4 * fk;
    
    H_ham(1,2) = conj(h0);
    H_ham(1,3) = conj(h4);
    H_ham(1,4) = h3;
    
    H_ham(2,1) = h0;
    H_ham(2,3) = gamma1;
    H_ham(2,4) = conj(h4);
    
    H_ham(3,1) = h4;
    H_ham(3,2) = gamma1;
    H_ham(3,4) = conj(h0);
    
    H_ham(4,1) = conj(h3);
    H_ham(4,2) = h4;
    H_ham(4,3) = h0;
    
    %% 检查厄米性
    % hem = helper_check_hermite(H_ham, 1e-8);
    % if hem == 0
    %   disp("monolayer Ham is not hermitian")
    % end
end

%% 作图
% figure
% plot(ak_len_array, eig_vals, 'b-','LineWidth', 2)
% ylim([-0.004, 0.004]) % 30meV
% grid on
% xlim([min(ak_len_array), max(ak_len_array)])
% ylim([-0.02, 0.02]) % 30meV % ylim([-0.2, 0.2]) % 30meV