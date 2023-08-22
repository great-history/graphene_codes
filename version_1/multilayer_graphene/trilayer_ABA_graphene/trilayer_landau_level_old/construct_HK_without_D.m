% 构造没有D field下的trilayer graphene H(B) under perpendicular magnetic field
function [HK_b, HK_m, HKp_b, HKp_m] = construct_HK_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, B_field, N_LL)
    % bilayer
    U1 = Delta2 + gamma2 / 2;
    U2 = Delta2 + delta + gamma5 / 2;
    U3 = -2*Delta2 + delta;
    U4 = -2*Delta2;
    % monolayer
    U5 = Delta2 - gamma2 / 2;
    U6 = Delta2 + delta - gamma5 / 2;
   
    E0 = (sqrt(2 * B_field) * gamma0 * 0.142 * 3) / (2 * 25.66);
%     E0 = (sqrt(B_field) * gamma0 * 0.142 * 3) / (2 * 25.66);
    E3 = E0 * sqrt(2) * gamma3 / gamma0;
    E4 = E0 * sqrt(2) * gamma4 / gamma0;
    
    % 因为D = 0,所以可以把HK分成两个之间完全无耦合的矩阵H_k_b和H_k_m,对HKp也是同理
    
    % construct the matrix for valley K
    % basis: bilayer basis
    %        |1>=(0,|0>,0,0)
    %        |2>=(0,0,0,|0>);
    %        |n>=1/sqrt(2)(|n-3>,|n-2>,0,0);   3<=n<=N_LL+3
    %        |N_LL+1+n>=1/sqrt(2)(|n-3>,-|n-2>,0,0);   3<=n<=N_LL+3
    %        |2N+2+n>=1/sqrt(2)(0,0,|n-3>,|n-2>);   3<=n<=N_LL+3
    %        |3N+3+n>=1/sqrt(2)(0,0,|n-3>,-|n-2>);   3<=n<=N_LL+3
    %        monolayer basis
    %        |4N+7>=(0,|0>)
    %        |4N+5+n>=1/sqrt(2)(|n-3>,|n-2>);   3<=n<=N_LL+3
    %*********************************
    % construct the matrix for valley K
    [HK_m, HKp_m] = construct_HK_monolayer(E0, U5, U6, N_LL);
    [HK_b, HKp_b] = construct_HK_bilayer(sqrt(2) * gamma1, E0, E3, E4, U1, U2, U3, U4, N_LL);
    
    % hem = check_hermite(HK_b, 1e-8);
    % if hem == 1
    %    disp("HK_b is hermitian")
    % end
    
    % don't need to construct the matrix between monolayer and bilayer
end

% example:
% gamma0 = 3100;
% gamma1 = 355;
% gamma2 = -20.7904;
% gamma3 = 346.7153;  
% gamma4 = 91.2578;
% gamma5 = 35.7223;
% D_field = 0;
% delta = 1.5+(gamma5-gamma2)/2;
% Delta1 = D_field * 0.1 * 1000;
% Delta2 = 1.8;
% N_LL = 30;
% B_field = 0.05;
% call the function
% [HK_b, HK_m, HKp_b, HKp_m] = construct_HK_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, B_field, N_LL)