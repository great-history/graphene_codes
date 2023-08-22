function [HK_m, HKp_m] = construct_HK_monolayer(E0, U5, U6, N_LL)
    HK_m = construct_HK_m(E0, U5, U6, N_LL);
    
    % particle-hole transformation:-PHP : construct the matrix for valley K_prime
    HKp_m = construct_HKp_m(E0, U5, U6, N_LL);
end

function HK_m = construct_HK_m(E0, U5, U6, N_LL)
    HK_m = zeros(2*N_LL+3, 2*N_LL+3);
    % monolayer basis
    % |1>=(0,|0>)
    % |n>=1/sqrt(2)(|n-2>,|n-1>);   2<=n<=N_LL+2
    % |N_LL+1+n>=1/sqrt(2)(|n-2>,-|n-1>);   2<=n<=N_LL+2
    HK_m(1,1) = U6;
    
    x = 1/2*(U5+U6);
    y = 1/2*(U5-U6);
    for i = 2:N_LL + 2
        % construct monolayer ham
        % |n>=1/sqrt(2)(|n-2>,|n-1>);   2<=n<=N_LL+2
        HK_m(i,i) = x + E0*sqrt(i-1);
        HK_m(i,N_LL+1+i) = y;
        HK_m(N_LL+1+i,i) = y;

        % |N_LL+1+n>=1/sqrt(2)(|n-2>,-|n-1>);   2<=n<=N_LL+2
        HK_m(N_LL+1+i,N_LL+1+i) = x - E0*sqrt(i-1);
    end
    
end

function HKp_m = construct_HKp_m(E0, U5, U6, N_LL)
    HKp_m = construct_HK_m(E0, -U6, -U5, N_LL);
    HKp_m = - HKp_m;
end