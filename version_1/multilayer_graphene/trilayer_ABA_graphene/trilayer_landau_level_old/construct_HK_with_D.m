function [HK,HKp] = construct_HK_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, B_field, N_LL)
    U1 = Delta2 + gamma2 / 2;
    U2 = Delta2 + delta + gamma5 / 2;
    U3 = -2*Delta2 + delta;
    U4 = -2*Delta2;
    U5 = Delta2 - gamma2 / 2;
    U6 = Delta2 + delta - gamma5 / 2;
   
    E0 = (sqrt(2 * B_field) * gamma0 * 0.142 * 3) / (2 * 25.66);
%     E0 = (sqrt(B_field) * gamma0 * 0.142 * 3) / (2 * 25.66);
    E3 = E0 * sqrt(2) * gamma3 / gamma0;
    E4 = E0 * sqrt(2) * gamma4 / gamma0;
    
    % 因为D = 0,所以可以把HK分成两个之间完全无耦合的矩阵H_k_b和H_k_m,对HKp也是同理
    
    % construct the matrix for valley K
    % basis: bilayer basis
    % |1>=(0,|0>,0,0);
    % |n>=1/sqrt(2)(|n-2>,|n-1>,0,0);   2<=n<=N_LL+2
    % |m+n>=1/sqrt(2)(|n-2>,-|n-1>,0,0);   2<=n<=N_LL+2
    
    % |p+1> = (0,0,0,|0>)
    % |p+n>=1/sqrt(2)(0,0,|n-2>,|n-1>);   2<=n<=N_LL+2
    
    % |q+n>=1/sqrt(2)(0,0,|n-2>,-|n-1>);   2<=n<=N_LL+2
    
    % monolayer basis
    % |1>=(0,|0>)
    % |n>=1/sqrt(2)(|n-2>,|n-1>);   2<=n<=N_LL+2
    % |N_LL+1+n>=1/sqrt(2)(|n-2>,-|n-1>);   2<=n<=N_LL+2
    %*********************************
    % construct the matrix for valley K
    [HK_m, HKp_m] = construct_HK_monolayer(E0, U5, U6, N_LL);
    [HK_b, HKp_b] = construct_HK_bilayer(sqrt(2) * gamma1, E0, E3, E4, U1, U2, U3, U4, N_LL);
    
    DK = DK_mb(Delta1, N_LL);
    DKp = DKp_mb(Delta1, N_LL);
    
    HK = [HK_m, DK; DK.',HK_b];
    HKp = [HKp_m, DKp; DKp.',HKp_b];
    
%     hem1 = check_hermite(HK,1e-8);
%     hem2 = check_hermite(HKp,1e-8);
%     if hem1 == 0 || hem2 == 0
%         disp("the Hamiltonian is not hermitian")
%     end
end

function DK = DK_mb(Delta1, N_LL)
    % basis: bilayer basis
    % |1>=(0,|0>,0,0);
    % |n>=1/sqrt(2)(|n-2>,|n-1>,0,0);   2<=n<=N_LL+2
    % |m+n>=1/sqrt(2)(|n-2>,-|n-1>,0,0);   2<=n<=N_LL+2
    % monolayer basis
    % |1>=(0,|0>)
    % |n>=1/sqrt(2)(|n-2>,|n-1>);   2<=n<=N_LL+2
    % |N_LL+1+n>=1/sqrt(2)(|n-2>,-|n-1>);   2<=n<=N_LL+2
    
    DK = zeros(2*N_LL+3,4*N_LL+6);
    DK(1,1) = Delta1;
    m = N_LL + 1;
    for n = 2:N_LL+2
        DK(n,n) = Delta1;
        DK(m+n,m+n) = Delta1;
    end
end

function DKp = DKp_mb(Delta1, N_LL)
    % bilayer basis
    % |p+1> = (0,0,0,|0>)
    % |p+n>=1/sqrt(2)(0,0,|n-2>,|n-1>);   2<=n<=N_LL+2
    % |q+n>=1/sqrt(2)(0,0,|n-2>,-|n-1>);   2<=n<=N_LL+2
    
    % monolayer basis
    % |1>=(0,|0>)
    % |n>=1/sqrt(2)(|n-2>,|n-1>);   2<=n<=N_LL+2
    % |N_LL+1+n>=1/sqrt(2)(|n-2>,-|n-1>);   2<=n<=N_LL+2
    
    DKp = zeros(2*N_LL+3,4*N_LL+6);
    m = N_LL + 1;
    p = 2 * N_LL + 3;
    q = 3 * N_LL + 4;
    
    DKp(1,p+1) = Delta1;
    for n = 2:N_LL+2
        DKp(n,p+n) = Delta1;
        DKp(m+n,q+n) = Delta1;
    end
end


% example:
% gamma0 = 3100;
% gamma1 = 355;
% gamma2 = -20.7904;
% gamma3 = 346.7153;  
% gamma4 = 91.2578;
% gamma5 = 35.7223;
% D_field = 1;
% delta = 1.5+(gamma5-gamma2)/2;
% Delta1 = D_field * 0.1 * 1000;
% Delta2 = 1.8;
% N_LL = 30;
% B_field = 0.05;
% % call the function
% [HK, HKp] = construct_HK_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, B_field, N_LL)