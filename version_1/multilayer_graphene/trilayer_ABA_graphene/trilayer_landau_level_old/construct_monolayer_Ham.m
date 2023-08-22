function H_m_k = construct_monolayer_Ham(Ea, Eb, gamma0, f_k)
    H_m_k = zeros(2,2);
    H_m_k(1,1) = Ea;
    H_m_k(2,2) = Eb;
    H_m_k(1,2) = - gamma0 * f_k;
    H_m_k(2,1) = conj(H_m_k(1,2));
    
    hem = helper_check_hermite(H_m_k, 1e-8);
    if hem == 0
       disp("monolayer Ham is not hermitian")
    end
end


% example
