function [HK_m_ham, HKp_m_ham, HK_b_ham, HKp_b_ham] = construct_trilayer_ABA_LL_six_bands_continuum_model(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, akx, aky)
    on_site_params_m = [- gamma2 / 2 + Delta2, delta - gamma5 / 2 + Delta2];
    on_site_params_b = [gamma2 / 2 + Delta2, delta + gamma5 / 2 + Delta2, delta - 2 * Delta2, - 2 * Delta2];
    [HK_m_ham, HKp_m_ham] = construct_monolayer_continuum_model(gamma0, on_site_params_m, akx, aky);
    [HK_b_ham, HKp_b_ham] = construct_bilayer_continuum_model(gamma0, sqrt(2) * gamma1, sqrt(2) * gamma3, sqrt(2) * gamma4, on_site_params_b, akx, aky);
    
    hem1 = helper_check_hermite(HK_m_ham, 1e-8);
    hem2 = helper_check_hermite(HKp_m_ham,1e-8);
    hem3 = helper_check_hermite(HK_b_ham,1e-8);
    hem4 = helper_check_hermite(HKp_b_ham,1e-8);
    if hem1 == 0 || hem2 == 0 || hem3 == 0 || hem4 == 0
        disp("the Hamiltonian is not hermitian")
    end
end

function [HK_ham, HKp_ham] = construct_monolayer_continuum_model(gamma0, on_site_params, akx, aky)
    HK_ham = zeros(2,2);
    HKp_ham = zeros(2,2);
    h0 = - sqrt(3) / 2 * gamma0 * (akx - 1j * aky);
    
    HK_ham(1, 2) = conj(h0);
    HK_ham(2, 1) = h0;
    HK_ham = HK_ham + diag(on_site_params);
    
    HKp_ham(1, 2) = - h0;
    HKp_ham(2, 1) = - conj(h0);
    HKp_ham = HKp_ham + diag(on_site_params);
end

function [HK_ham, HKp_ham] = construct_bilayer_continuum_model(gamma0, gamma1, gamma3, gamma4, on_site_params, akx, aky)
    % ext_mat : 描述外界引入的项，包括 电位移场 / Ising SOCs
    %% 构造哈密顿量（紧束缚模型）
    % 不包括除了delta之外的onsite potential
    % basis : A1, B1, A2, B2
    h0 = - sqrt(3) / 2 * gamma0 * (akx - 1j * aky);
    h3 = - sqrt(3) / 2 * gamma3 * (akx - 1j * aky);
    h4 = - sqrt(3) / 2 * gamma4 * (akx - 1j * aky);
    
    %% vally K = (4 * pi / (3 * sqrt(3)), 0)
    HK_ham = diag(on_site_params);
    HK_ham(2,3) = gamma1;
    HK_ham(3,2) = gamma1;
    
    HK_ham(1,2) = conj(h0);
    HK_ham(1,3) = - conj(h4);
    HK_ham(1,4) = h3;
    
    HK_ham(2,1) = h0;
    HK_ham(2,4) = - conj(h4);
    
    HK_ham(3,1) = - h4;
    HK_ham(3,4) = conj(h0);
    
    HK_ham(4,1) = conj(h3);
    HK_ham(4,2) = - h4;
    HK_ham(4,3) = h0;
    
    %% 检查厄米性
    hem = helper_check_hermite(HK_ham, 1e-8);
    if hem == 0
      disp("monolayer Ham is not hermitian")
    end
    
    %% vally Kp = (-4 * pi / (3 * sqrt(3)), 0)
    HKp_ham = diag(on_site_params);
    HKp_ham(2,3) = gamma1;
    HKp_ham(3,2) = gamma1;
    
    HKp_ham(1,2) = - h0;
    HKp_ham(1,3) = h4;
    HKp_ham(1,4) = - conj(h3);
    
    HKp_ham(2,1) = - conj(h0);
    HKp_ham(2,4) = h4;
    
    HKp_ham(3,1) = conj(h4);
    HKp_ham(3,4) = - h0;
    
    HKp_ham(4,1) = - h3;
    HKp_ham(4,2) = conj(h4);
    HKp_ham(4,3) = - conj(h0);
end