function [H_m_K, H_m_Kp] = construct_monolayer_effective_Ham(Ea, Eb, gamma0, akx, aky)
    Im = 1i;
    % K处的低能有效模型
    H_m_K = zeros(2,2);
    H_m_K(1,1) = Ea;
    H_m_K(2,2) = Eb;
    H_m_K(1,2) = - sqrt(3) / 2 * gamma0 * (-akx + Im * aky);
    % H_m_K(1,2) = - gamma0 * (-akx + Im * aky);
    H_m_K(2,1) = conj(H_m_K(1,2));
    
    % K^prime处的低能有效模型
    H_m_Kp = zeros(2,2);
    H_m_Kp(1,1) = Ea;
    H_m_Kp(2,2) = Eb;
    H_m_Kp(1,2) = - sqrt(3) / 2 * gamma0 * (akx + Im * aky);
    % H_m_Kp(1,2) = - gamma0 * (akx + Im * aky);
    H_m_Kp(2,1) = conj(H_m_Kp(1,2));
    
    hem1 = helper_check_hermite(H_m_K, 1e-8);
    hem2 = helper_check_hermite(H_m_Kp, 1e-8);
    if hem1 == 0 || hem2 == 0
       disp("monolayer Ham is not hermitian")
    end
end