% band structure parameters
gamma0 = 3100;
gamma1 = 390;
gamma2 = -28;
gamma3 = 315;  
gamma4 = 41;
gamma5 = 50;
delta = 46;

Delta1 = 0; % 0 / 25 / 50 / 150 / 180 / 250
Delta2 = 0;

%　kxy_mesh
Nx = 201;
Ny = 201;
ak_start = [-0.125,-0.125];
ak_end = [0.125,0.125];

akxs = linspace(ak_start(1), ak_end(1), Nx);
akys = linspace(ak_start(2), ak_end(2), Ny);

for i = 1:1
    akx = akxs(i);
    aky = akys(1);
    [H_t_K, H_t_Kp] = construct_trilayer_effective_Ham_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, akx, aky);
    
    [HK_m_ham, HKp_m_ham, HK_b_ham, HKp_b_ham] = construct_trilayer_ABA_six_bands_continuum_model(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, akx, aky);
end