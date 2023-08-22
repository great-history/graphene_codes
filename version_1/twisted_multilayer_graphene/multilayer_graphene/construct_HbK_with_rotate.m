function H_b_K = construct_HbK_with_rotate(gamma0, gamma1, gamma3, gamma4, delta_dimer, akx, aky, angle, flag_chiral)
    % 这个代码并没有考虑面内的次紧邻跃迁
    % vf = 2.1354; % eV为单位
    % flag_chiral : 描述了手性:是AB还是BA
    vf = gamma0 * sqrt(3) / 2;
    v3 = gamma3 * sqrt(3) / 2;
    v4 = gamma4 * sqrt(3) / 2;
    
    phase = exp(-1j*angle); % angle已经转化为弧度制了
    
    %% 单层
    H_m_K_1 = zeros(2,2);
    H_m_K_1(1,2) = vf * (akx + 1j * aky) * phase;
    H_m_K_1(2,1) = conj(H_m_K_1(1,2));
    H_m_K_1(2,2) = delta_dimer;
    
    H_m_K_2 = zeros(2,2);
    H_m_K_2(1,2) = vf * (akx + 1j * aky) * phase;
    H_m_K_2(2,1) = conj(H_m_K_2(1,2));
    H_m_K_2(1,1) = delta_dimer;
    
    H_inter = zeros(2,2);
    H_inter(1,1) = - v4 * (akx - 1j * aky) * conj(phase);
    H_inter(1,2) = gamma1;
    H_inter(2,1) = - v3 * (akx + 1j * aky) * phase;
    H_inter(2,2) = - v4 * (akx - 1j * aky) * conj(phase);
    
    %% 得到哈密顿量
    H_b_K = zeros(4,4);
    if flag_chiral % AB-stacking
        H_b_K(1:2,1:2) = H_m_K_1;
        H_b_K(3:4,3:4) = H_m_K_2;
        
        H_b_K(3:4,1:2) = H_inter;
        H_b_K(1:2,3:4) = H_inter';
    else % BA-stacking
        H_b_K(1:2,1:2) = H_m_K_2;
        H_b_K(3:4,3:4) = H_m_K_1;
        
        H_b_K(3:4,1:2) = H_inter';
        H_b_K(1:2,3:4) = H_inter;
    end
end