function [H_m_K, H_b_K, H_m_Kp, H_b_Kp] = construct_trilayer_effective_Ham_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, akx, aky)
    [H_m_K, H_m_Kp] = construct_monolayer_effective_Ham(Delta2 - gamma2 / 2, Delta2 + delta - gamma5 / 2, gamma0, akx, aky);
    on_sites = [Delta2 + gamma2 / 2, -2*Delta2, -2*Delta2 + delta, Delta2 + delta + gamma5 /2];
    [H_b_K, H_b_Kp] = construct_bilayer_effective_Ham(on_sites, gamma0, sqrt(2) * gamma1, sqrt(2) * gamma3, sqrt(2) * gamma4, akx, aky);

%     hem1 = check_hermite(H_m_K,1e-8);
%     hem2 = check_hermite(H_m_Kp,1e-8);
%     if hem1 == 0 || hem2 == 0
%         disp("the Hamiltonian is not hermitian")
%     end
end

function [H_b_K, H_b_Kp] = construct_bilayer_effective_Ham(on_sites, gamma0, sqrt2_gamma1, sqrt2_gamma3, sqrt2_gamma4, akx, aky)
    H_b_K = zeros(4,4);
    
    H_b_K(1,1) = on_sites(1);
    H_b_K(2,2) = on_sites(2);
    H_b_K(3,3) = on_sites(3);
    H_b_K(4,4) = on_sites(4);
    
    H_b_K(3,4) = sqrt2_gamma1;
    H_b_K(4,3) = sqrt2_gamma1;
    
    Im = 1i;
    f_K_eff = (-akx + Im * aky);
    x0 = gamma0 * f_K_eff;
    x3 = sqrt2_gamma3 * f_K_eff;
    x4 = sqrt2_gamma4 * f_K_eff;
    
    H_b_K(1,4) = - x0;
    H_b_K(4,1) = - conj(x0);
    H_b_K(2,3) = - conj(x0);
    H_b_K(3,2) = - x0;
    
    H_b_K(1,2) = - conj(x3);
    H_b_K(2,1) = - x3;
    
    H_b_K(1,3) = x4;
    H_b_K(3,1) = conj(x4);
    H_b_K(2,4) = conj(x4);
    H_b_K(4,2) = x4;
    % ---------------------------------------------------------------------
    H_b_Kp = zeros(4,4);
    H_b_Kp(1,1) = on_sites(1);
    H_b_Kp(2,2) = on_sites(2);
    H_b_Kp(3,3) = on_sites(3);
    H_b_Kp(4,4) = on_sites(4);
    
    H_b_Kp(3,4) = sqrt2_gamma1;
    H_b_Kp(4,3) = sqrt2_gamma1;

    f_Kp_eff = (akx + Im * aky);
    x0 = gamma0 * f_Kp_eff;
    x3 = sqrt2_gamma3 * f_Kp_eff;
    x4 = sqrt2_gamma4 * f_Kp_eff;

    H_b_Kp(1,4) = - x0;
    H_b_Kp(4,1) = - conj(x0);
    H_b_Kp(2,3) = - conj(x0);
    H_b_Kp(3,2) = - x0;
    
    H_b_Kp(1,2) = - conj(x3);
    H_b_Kp(2,1) = - x3;
    
    H_b_Kp(1,3) = x4;
    H_b_Kp(3,1) = conj(x4);
    H_b_Kp(2,4) = conj(x4);
    H_b_Kp(4,2) = x4;
    
    H_b_Kp(1,1) = on_sites(1);
    H_b_Kp(2,2) = on_sites(2);
    H_b_Kp(1,1) = on_sites(1);
    H_b_Kp(2,2) = on_sites(2);
    
    H_b_Kp(3,4) = sqrt2_gamma1;
    H_b_Kp(4,3) = sqrt2_gamma1;
    
    hem1 = check_hermite(H_b_K, 1e-8);
    hem2 = check_hermite(H_b_Kp, 1e-8);
    if hem1 == 0 || hem2 == 0
       disp("monolayer Ham is not hermitian")
    end
end

