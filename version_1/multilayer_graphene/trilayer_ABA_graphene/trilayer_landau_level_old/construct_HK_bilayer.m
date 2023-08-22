function [HK_b, HKp_b] = construct_HK_bilayer(gamma1, E0, E3, E4, U1, U2, U3, U4, N_LL)
    HK_b = construct_HK_b(gamma1, E0, E3, E4, U1, U2, U3, U4, N_LL);
    
    % particle-hole transformation:-PHP : construct the matrix for valley K_prime
    HKp_b = construct_HKp_b(gamma1, E0, E3, E4, U1, U2, U3, U4, N_LL);
end

function HK_b = construct_HK_b(gamma1, E0, E3, E4, U1, U2, U3, U4, N_LL)
    HK_b = zeros(4*N_LL+6, 4*N_LL+6);
    m = N_LL + 1;
    p = 2 * N_LL + 3;
    q = 3 * N_LL + 4;
    
    x = (U1+U2) / 2;
    z = (U1-U2) / 2;
    s = (U3 + U4) / 2;
    t = (U3 - U4) / 2;
    g = gamma1 / 2;
    e4 = E4 / 2;
    e3 = E3 / 2;
    
    % basis: bilayer basis
    % |1>=(0,|0>,0,0);
    % |n>=1/sqrt(2)(|n-2>,|n-1>,0,0);   2<=n<=N_LL+2
    % |m+n>=1/sqrt(2)(|n-2>,-|n-1>,0,0);   2<=n<=N_LL+2
    
    % |p+1> = (0,0,0,|0>)
    % |p+n>=1/sqrt(2)(0,0,|n-2>,|n-1>);   2<=n<=N_LL+2
    
    % |q+n>=1/sqrt(2)(0,0,|n-2>,-|n-1>);   2<=n<=N_LL+2
    
    % 先把头和尾处理掉
    % |1>=(0,|0>,0,0)
    HK_b(1,1) = U2;
    HK_b(1,p+2) = (gamma1 - E4) / sqrt(2);
    HK_b(p+2,1) = (gamma1 - E4) / sqrt(2);
    HK_b(1,q+2) = (gamma1 + E4) / sqrt(2);
    HK_b(q+2,1) = (gamma1 + E4) / sqrt(2);
    
    % |1>=(0,|0>,0,0)
    HK_b(p+1,p+1) = U4;
    HK_b(p+1,3) = E3 / sqrt(2);
    HK_b(3,p+1) = E3 / sqrt(2);
    HK_b(p+1,m+3) = E3 / sqrt(2);
    HK_b(m+3,p+1) = E3 / sqrt(2);
    
    % 中间部分的循环
    for n = 2:N_LL
        % |n>=1/sqrt(2)(|n-2>,|n-1>,0,0);   2<=n<=N_LL
        HK_b(n,n) = x + E0 * sqrt(n-1);
        
        HK_b(n,m+n) = z;
        HK_b(m+n,n) = z;
        
        HK_b(n,p+n+1) = g - e4 * (sqrt(n-1) + sqrt(n));
        HK_b(p+n+1,n) = g - e4 * (sqrt(n-1) + sqrt(n));
        
        HK_b(n,q+n+1) = g + e4 * (sqrt(n) - sqrt(n-1));
        HK_b(q+n+1,n) = g + e4 * (sqrt(n) - sqrt(n-1));
        
        % |m+n>=1/sqrt(2)(|n-2>,-|n-1>,0,0);   2<=n<=N_LL
        HK_b(m+n,m+n) = x - E0 * sqrt(n-1);
        
        HK_b(m+n,p+n+1) = - g + e4 * (sqrt(n) - sqrt(n - 1));
        HK_b(p+n+1,m+n) = - g + e4 * (sqrt(n) - sqrt(n - 1));
        
        HK_b(m+n,q+n+1) = - g - e4 * (sqrt(n) + sqrt(n - 1));
        HK_b(q+n+1,m+n) = - g - e4 * (sqrt(n) + sqrt(n - 1));
        
        % |p+n>=1/sqrt(2)(0,0,|n-2>,|n-1>);   2<=n<=N_LL
        HK_b(p+n,p+n) = s + E0 * sqrt(n-1);
        
        HK_b(p+n,q+n) = t;
        HK_b(q+n,p+n) = t;
        
        HK_b(p+n,n+2) = e3 * sqrt(n);
        HK_b(n+2,p+n) = e3 * sqrt(n);
        
        HK_b(p+n,m+n+2) = e3 * sqrt(n);
        HK_b(m+n+2,p+n) = e3 * sqrt(n);
        
        % |q+n>=1/sqrt(2)(0,0,|n-2>,-|n-1>);   2<=n<=N_LL
        HK_b(q+n,q+n) = s - E0 * sqrt(n-1);
        
        HK_b(q+n,n+2) = - e3 * sqrt(n);
        HK_b(n+2,q+n) = - e3 * sqrt(n);
        
        HK_b(q+n,m+n+2) = - e3 * sqrt(n);
        HK_b(m+n+2,q+n) = - e3 * sqrt(n);
        
    end
    
    % 再把尾巴处理掉
    % |N_LL+1> = 1/sqrt(2)(|N_LL-1>,|N_LL>,0,0)
    HK_b(N_LL+1,N_LL+1) = x + E0 * sqrt(N_LL);
    
    HK_b(N_LL+1,m+N_LL+1) = z;
    HK_b(m+N_LL+1,N_LL+1) = z;
    
    HK_b(N_LL+1,p+N_LL+2) = g - e4 * (sqrt(N_LL+1) + sqrt(N_LL));
    HK_b(p+N_LL+2,N_LL+1) = g - e4 * (sqrt(N_LL+1) + sqrt(N_LL));
    
    HK_b(N_LL+1,q+N_LL+2) = g + e4 * (sqrt(N_LL+1)-sqrt(N_LL));
    HK_b(q+N_LL+2,N_LL+1) = g + e4 * (sqrt(N_LL+1)-sqrt(N_LL));
    
    % |m+N_LL+1> = 1/sqrt(2)(|N_LL-1>,-|N_LL>,0,0)
    HK_b(m+N_LL+1,m+N_LL+1) = x - E0 * sqrt(N_LL);
    
    HK_b(m+N_LL+1,p+N_LL+2) = - g + e4 * (sqrt(N_LL+1)-sqrt(N_LL));
    HK_b(p+N_LL+2,m+N_LL+1) = - g + e4 * (sqrt(N_LL+1)-sqrt(N_LL));
    
    HK_b(m+N_LL+1,q+N_LL+2) = - g - e4 * (sqrt(N_LL)+sqrt(N_LL+1));
    HK_b(q+N_LL+2,m+N_LL+1) = - g - e4 * (sqrt(N_LL)+sqrt(N_LL+1));
    
    % |p+N_LL+1> = (0,0,|N_LL-1>,|N_LL>)
    HK_b(p+N_LL+1,p+N_LL+1) = s + E0 * sqrt(N_LL);
    
    HK_b(p+N_LL+1,q+N_LL+1) = t;
    HK_b(q+N_LL+1,p+N_LL+1) = t;
    
    % |q+N_LL+1> = (0,0,|N_LL-1>,-|N_LL>)
    HK_b(q+N_LL+1,q+N_LL+1) = s - E0 * sqrt(N_LL);
        
    % |N_LL+2> = 1/sqrt(2)(|N_LL>,|N_LL+1>,0,0)
    HK_b(N_LL+2,N_LL+2) = x + E0 * sqrt(N_LL+1);
    
    HK_b(N_LL+2,m+N_LL+2) = z;
    HK_b(m+N_LL+2,N_LL+2) = z;
    
    % |m+N_LL+2> = (|N_LL>,-|N_LL+1>,0,0)
    HK_b(m+N_LL+2,m+N_LL+2) = x - E0 * sqrt(N_LL+1);
    
    % |p+N_LL+2> = (0,0,|N_LL>,|N_LL+1>)
    HK_b(p+N_LL+2,p+N_LL+2) = s + E0 * sqrt(N_LL+1);
    
    HK_b(p+N_LL+2,q+N_LL+2) = t;
    HK_b(q+N_LL+2,p+N_LL+2) = t;
    
    % |q+N_LL+2> = (0,0,|N_LL>,-|N_LL+1>)
    HK_b(q+N_LL+2,q+N_LL+2) = s - E0 * sqrt(N_LL+1);
end

function HKp_b = construct_HKp_b(gamma1, E0, E3, E4, U1, U2, U3, U4, N_LL)
    HKp_b = construct_HK_b(-gamma1, E0, E3, E4, -U4, -U3, -U2, -U1, N_LL);
    HKp_b = - HKp_b;
end