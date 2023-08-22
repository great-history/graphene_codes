function [HK_ham, HKp_ham] = construct_bilayer_continuum_model(gamma0, gamma1, gamma3, gamma4, delta, ext_mat_K, ext_mat_Kp, akx, aky)
    % ext_mat : 描述外界引入的项，包括 电位移场 / Ising SOCs
    %% 构造哈密顿量（紧束缚模型）
    % 不包括除了delta之外的onsite potential
    % basis : A1, B1, A2, B2
    h0 = - sqrt(3) / 2 * gamma0 * (akx - 1j * aky);
    h3 = - sqrt(3) / 2 * gamma3 * (akx - 1j * aky);
    h4 = - sqrt(3) / 2 * gamma4 * (akx - 1j * aky);
    
    %% vally K = (4 * pi / (3 * sqrt(3)), 0)
    HK_ham = zeros(4,4);
    HK_ham(2,2) = delta;
    HK_ham(3,3) = delta;
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
    
    HK_ham = HK_ham + ext_mat_K; % 加入external matrix accounts for out-of-plane electric field || Ising SOC induced by proximity
    
    %% 检查厄米性
    hem = helper_check_hermite(HK_ham, 1e-8);
    if hem == 0
      disp("monolayer Ham is not hermitian")
    end
    
    %% vally Kp = (-4 * pi / (3 * sqrt(3)), 0)
    HKp_ham = zeros(4,4);
    HKp_ham(2,2) = delta;
    HKp_ham(3,3) = delta;
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
    
    HKp_ham = HKp_ham + ext_mat_Kp;
end


% function [HK_ham, HKp_ham] = construct_bilayer_continuum_model(gamma0, gamma1, gamma3, gamma4, delta_dimer, akx, aky)
%     %% 构造哈密顿量（紧束缚模型）
%     % 不包括除了delta之外的onsite potential
%     % basis : A1, B1, A2, B2
%     h0 = - sqrt(3) / 2 * gamma0 * (akx - 1j * aky);
%     h3 = - sqrt(3) / 2 * gamma3 * (akx - 1j * aky);
%     h4 = - sqrt(3) / 2 * gamma4 * (akx - 1j * aky);
%     
%     %% vally K = (4 * pi / (3 * sqrt(3)), 0)
%     HK_ham = zeros(4,4);
%     HK_ham(2,2) = delta_dimer;
%     HK_ham(3,3) = delta_dimer;
%     HK_ham(2,3) = gamma1;
%     HK_ham(3,2) = gamma1;
%     
%     HK_ham(1,2) = conj(h0);
%     HK_ham(1,3) = conj(h4);
%     HK_ham(1,4) = h3;
%     
%     HK_ham(2,1) = h0;
%     HK_ham(2,4) = conj(h4);
%     
%     HK_ham(3,1) = h4;
%     HK_ham(3,4) = conj(h0);
%     
%     HK_ham(4,1) = conj(h3);
%     HK_ham(4,2) = h4;
%     HK_ham(4,3) = h0;
%     
%     
%     %% 检查厄米性
%     % hem = helper_check_hermite(HK_ham, 1e-8);
%     % if hem == 0
%     %   disp("monolayer Ham is not hermitian")
%     % end
%     
%     %% vally Kp = (-4 * pi / (3 * sqrt(3)), 0)
%     HKp_ham = zeros(4,4);
%     HKp_ham(2,2) = delta_dimer;
%     HKp_ham(3,3) = delta_dimer;
%     HKp_ham(2,3) = gamma1;
%     HKp_ham(3,2) = gamma1;
%     
%     HKp_ham(1,2) = -h0;
%     HKp_ham(1,3) = -h4;
%     HKp_ham(1,4) = -conj(h3);
%     
%     HKp_ham(2,1) = -conj(h0);
%     HKp_ham(2,4) = -h4;
%     
%     HKp_ham(3,1) = -conj(h4);
%     HKp_ham(3,4) = -h0;
%     
%     HKp_ham(4,1) = -h3;
%     HKp_ham(4,2) = -conj(h4);
%     HKp_ham(4,3) = -conj(h0);
% end