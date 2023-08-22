gamma0 = 3165;
gamma1 = 391;
gamma2 = -20.2;
gamma3 = 315.15;
gamma4 = 44.24;

gamma5 = 38.5;
delta = 37.5;
delta1 = 0;
delta2 = 0;

U1 = delta2 + gamma2 / 2;
U2 = delta2 + delta + gamma5 / 2;
U3 = -2*delta2 + delta;
U4 = -2*delta2;

U5 = delta2 - gamma2 / 2;
U6 = delta2 + delta - gamma5 / 2;
vel = gamma0 * 0.142 * 3 / 2;

% band structure calculation
% construct the Hamiltonian

B_start = 1;
B_final = 1;
B_steps = 1;
B_fields = linspace(B_start, B_final, B_steps);

N_LL = 21;
N_dim = N_LL * 6 + 9;
LL_level_K = zeros(B_steps, N_dim);
LL_level_Kp = zeros(B_steps, N_dim);

for B_index = 1:B_steps
    B = B_fields(B_index);
    E0 = sqrt(2) * vel * sqrt(B) / 25.66;
    E3 = E0 * gamma3 / gamma0;
    E4 = E0 * gamma4 / gamma0;
    
    % construct the matrix for valley K
    % basis: bilayer basis
    %        |1>=(0,|0>,0,0,0,0)
    %        |2>=(0,0,0,|0>,0,0);
    %        |n>=1/sqrt(2)(|n-3>,|n-2>,0,0,0,0);   3<=n<=N_LL+3
    %        |N_LL+1+n>=1/sqrt(2)(|n-3>,-|n-2>,0,0,0,0);   3<=n<=N_LL+3
    %        |2N+2+n>=1/sqrt(2)(0,0,|n-3>,|n-2>,0,0);   3<=n<=N_LL+3
    %        |3N+3+n>=1/sqrt(2)(0,0,|n-3>,-|n-2>,0,0);   3<=n<=N_LL+3
    %        monolayer basis
    %        |4N+7>=(0,0,0,0,0,|0>)
    %        |4N+5+n>=1/sqrt(2)(0,0,0,0,|n-3>,|n-2>);   3<=n<=N_LL+3
    %*********************************
    HK = zeros(6*N_LL+9);
    %        |4N+7>=(0,0,0,0,0,|0>)
    HK(4*N_LL+7,4*N_LL+7) = U6;
    
    %        |1>=(0,|0>,0,0,0,0)
    %        |2>=(0,0,0,|0>,0,0);
    HK(2,2) = U2;
    HK(2,2*N_LL+2+3) = (gamma1 - E4);
    HK(2*N_LL+2+3,2) = (gamma1 - E4);
    HK(2,3*N_LL+3+3) = (gamma1 + E4);
    HK(3*N_LL+3+3,2) = (gamma1 + E4);
    HK(1,1) = U4;
    HK(1,4) = E3;
    HK(4,1) = E3;
    HK(1,N_LL+1+4) = E3;
    HK(N_LL+1+4,1) = E3;

    HK(2,4*N_LL+7) = delta1;
    
    for i = 3:N_LL+3
        %        |n>=1/sqrt(2)(|n-3>,|n-2>,0,0,0,0);   3<=n<=N_LL+3
        HK(i,i) = 1/2*(U1+U2) + E0 * sqrt(i-2);
        HK(i,N_LL+1+i) = 1/2*(U1-U2);
        HK(N_LL+1+i,i) = 1/2*(U1-U2);
        if i < N_LL+3
            HK(i,2*N_LL+3+i) = 1/sqrt(2) * (gamma1 - E4 *(sqrt(i-2)+sqrt(i-1)));
            HK(i,3*N_LL+4+i) = 1/sqrt(2) * (gamma1 - E4 *(sqrt(i-2)-sqrt(i-1)));
            HK(2*N_LL+3+i,i) = 1/sqrt(2) * (gamma1 - E4 *(sqrt(i-2)+sqrt(i-1)));
            HK(3*N_LL+4+i,i) = 1/sqrt(2) * (gamma1 - E4 *(sqrt(i-2)-sqrt(i-1)));
        end
        if i >= 5
            HK(i,2*N_LL+i) = 1/sqrt(2) * E3 * sqrt(i-3);
            HK(i,3*N_LL+1+i) = -1/sqrt(2) * E3 * sqrt(i-3);
            HK(2*N_LL+i,i) = 1/sqrt(2) * E3 * sqrt(i-3);
            HK(3*N_LL+1+i,i) = -1/sqrt(2) * E3 * sqrt(i-3);
        end
        
        %        |N_LL+1+n>=1/sqrt(2)(|n-3>,-|n-2>,0,0,0,0);   3<=n<=N_LL+3
        HK(N_LL+1+i, N_LL+1+i) = 1/2*(U1+U2) - E0*sqrt(i-2);
        if i < N_LL+3
            HK(N_LL+1+i, 2*N_LL+3+i) = -1/sqrt(2) * (gamma1 + E4*(sqrt(i-2)-sqrt(i-1)));
            HK(N_LL+1+i, 3*N_LL+4+i) = -1/sqrt(2) * (gamma1 + E4*(sqrt(i-2)+sqrt(i-1)));
            HK(2*N_LL+3+i, N_LL+1+i) = -1/sqrt(2) * (gamma1 + E4*(sqrt(i-2)-sqrt(i-1)));
            HK(3*N_LL+4+i, N_LL+1+i) = -1/sqrt(2) * (gamma1 + E4*(sqrt(i-2)+sqrt(i-1)));
        end
        if i >= 5
            HK(N_LL+1+i, 2*N_LL+i) = 1/sqrt(2) * E3 * sqrt(i-3);
            HK(N_LL+1+i, 3*N_LL+1+i) = -1/sqrt(2) * E3 * sqrt(i-3);
            HK(2*N_LL+i, N_LL+1+i) = 1/sqrt(2) * E3 * sqrt(i-3);
            HK(3*N_LL+1+i, N_LL+1+i) = -1/sqrt(2) * E3 * sqrt(i-3);
        end
        
        %        |2N+2+n>=1/sqrt(2)(0,0,|n-3>,|n-2>,0,0);   3<=n<=N_LL+3
        HK(2*N_LL+2+i,2*N_LL+2+i) = 1/2*(U3+U4) + E0*sqrt(i-2);
        HK(2*N_LL+2+i,3*N_LL+3+i) = 1/2*(U3-U4);
        HK(3*N_LL+3+i,2*N_LL+2+i) = 1/2*(U3-U4);
        
        %        |3N+3+n>=1/sqrt(2)(0,0,|n-3>,-|n-2>,0,0);   3<=n<=N_LL+3
        HK(3*N_LL+3+i,3*N_LL+3+i) = 1/2*(U3+U4) - E0*sqrt(i-2);
        
        %        |4N+5+n>=1/sqrt(2)(0,0,0,0,|n-3>,|n-2>);   3<=n<=N_LL+3
        HK(4*N_LL+5+i,4*N_LL+5+i) = 1/2*(U5+U6) + E0*sqrt(i-2);
        HK(4*N_LL+5+i,5*N_LL+6+i) = 1/2*(U5-U6);
        HK(5*N_LL+6+i,4*N_LL+5+i) = 1/2*(U5-U6);
        %        |4N+5+n>=1/sqrt(2)(0,0,0,0,|n-3>,|n-2>);   3<=n<=N_LL+3
        HK(5*N_LL+6+i,5*N_LL+6+i) = 1/2*(U5+U6) - E0*sqrt(i-2);
        
        % bilayer-monolayer hyb
        HK(i,4*N_LL+5+i) = delta1;
        HK(4*N_LL+5+i,i) = delta1;
        HK(N_LL+1+i,5*N_LL+6+i) = delta1;
        HK(5*N_LL+6+i,N_LL+1+i) = delta1;
    end
    
     %Diagnolize
    if delta1 == 0
        Q_now = zeros(6*N_LL+9);
        D_now = zeros(6*N_LL+9,1);
        [Q, D] = eig(HK(1:4*N_LL+6,1:4*N_LL+6));
        % sort eigen states to follow the continuous evolution
        Q_now(1:4*N_LL+6,1:4*N_LL+6) = Q;
        D_now(1:4*N_LL+6,1) = diag(D);
        clear Q D

        if B_index >= 2
            for q = 1:4*N_LL+6
                for j = q:4*N_LL+6
                    in_prod = abs(dot(Q_last(:,q), Q_now(:,j)));
                    if in_prod > 0.8
                        vec_ = Q_now(:,j);
                        Q_now(:,j) = Q_now(:,q);
                        Q_now(:,q) = vec_;
                        val_ = D_now(j);
                        D_now(j) = D_now(q);
                        D_now(q) = val_;
                        break
                    end
                end
            end
        end
        
        [Q, D] = eig(HK(4*N_LL+7:6*N_LL+9,4*N_LL+7:6*N_LL+9));
        % sort eigen states to follow the continuous evolution
        Q_now(4*N_LL+7:6*N_LL+9,4*N_LL+7:6*N_LL+9) = Q;
        D_now(4*N_LL+7:6*N_LL+9,1) = diag(D);
        clear Q D
        
        if B_index >= 2
            for q = 4*N_LL+7:6*N_LL+9
                for j = q:6*N_LL+9
                    in_prod = abs(dot(Q_last(:,q), Q_now(:,j)));
                    if in_prod > 0.8
                        vec_ = Q_now(:,j);
                        Q_now(:,j) = Q_now(:,q);
                        Q_now(:,q) = vec_;
                        val_ = D_now(j);
                        D_now(j) = D_now(q);
                        D_now(q) = val_;
                        break
                    end
                end
            end
        end
        
        Q_last = Q_now;
        LL_level_K(B_index,:) = D_now;
        clear Q_now D_now
    else
        [Q, D] = eig(HK);
        % sort eigen states to follow the continuous evolution
        Q_now = Q;
        D_now = diag(D);
        clear Q D

        if B_index >= 2
            for q = 1:6*N_LL+9
                for j = q:6*N_LL+9
                    in_prod = abs(dot(Q_last(:,q), Q_now(:,j)));
                    if in_prod > 0.8
                        vec_ = Q_now(:,j);
                        Q_now(:,j) = Q_now(:,q);
                        Q_now(:,q) = vec_;
                        val_ = D_now(j);
                        D_now(j) = D_now(q);
                        D_now(q) = val_;
                        break
                    end
                end
            end
        end

        Q_last = Q_now;
        LL_level_K(B_index,:) = D_now;
        clear Q_now D_now
    end
    
    
    % construct the matrix for valley K_prime
    % particle-hole transformation:-PHP
    HKp = zeros(N_dim);
    %        |4N+7>=(0,0,0,0,0,|0>)
    HKp(4*N_LL+7,4*N_LL+7) = U5;
    
    %        |1>=(0,|0>,0,0,0,0)
    %        |2>=(0,0,0,|0>,0,0);
    HKp(2,2) = U2;
    HKp(2,2*N_LL+2+3) = (gamma1 + E4);
    HKp(2*N_LL+2+3,2) = (gamma1 + E4);
    HKp(2,3*N_LL+3+3) = (gamma1 - E4);
    HKp(3*N_LL+3+3,2) = (gamma1 - E4);
    HKp(1,1) = U1;
    HKp(1,4) = -E3;
    HKp(4,1) = -E3;
    HKp(1,N_LL+1+4) = -E3;
    HKp(N_LL+1+4,1) = -E3;
    
    HKp(1,4*N_LL+7) = delta1;
    for i = 3:N_LL+3
        %        |n>=1/sqrt(2)(|n-3>,|n-2>,0,0,0,0);   3<=n<=N_LL+3
        HKp(i,i) = 1/2*(U3+U4) - E0 * sqrt(i-2);
        HKp(i,N_LL+1+i) = 1/2*(U4-U3);
        HKp(N_LL+1+i,i) = 1/2*(U4-U3);
        if i < N_LL+3
            HKp(i,2*N_LL+3+i) = 1/sqrt(2) * (gamma1 + E4 *(sqrt(i-2)+sqrt(i-1)));
            HKp(i,3*N_LL+4+i) = 1/sqrt(2) * (gamma1 + E4 *(sqrt(i-2)-sqrt(i-1)));
            HKp(2*N_LL+3+i,i) = 1/sqrt(2) * (gamma1 + E4 *(sqrt(i-2)+sqrt(i-1)));
            HKp(3*N_LL+4+i,i) = 1/sqrt(2) * (gamma1 + E4 *(sqrt(i-2)-sqrt(i-1)));
        end
        if i >= 5
            HKp(i,2*N_LL+i) = -1/sqrt(2) * E3 * sqrt(i-3);
            HKp(i,3*N_LL+1+i) = 1/sqrt(2) * E3 * sqrt(i-3);
            HKp(2*N_LL+i,i) = -1/sqrt(2) * E3 * sqrt(i-3);
            HKp(3*N_LL+1+i,i) = 1/sqrt(2) * E3 * sqrt(i-3);
        end
        
        %        |N_LL+1+n>=1/sqrt(2)(|n-3>,-|n-2>,0,0,0,0);   3<=n<=N_LL+3
        HKp(N_LL+1+i, N_LL+1+i) = 1/2*(U4+U3) + E0*sqrt(i-2);
        if i < N_LL+3
            HKp(N_LL+1+i, 2*N_LL+3+i) = -1/sqrt(2) * (gamma1 - E4*(sqrt(i-2)-sqrt(i-1)));
            HKp(N_LL+1+i, 3*N_LL+4+i) = -1/sqrt(2) * (gamma1 - E4*(sqrt(i-2)+sqrt(i-1)));
            HKp(2*N_LL+3+i, N_LL+1+i) = -1/sqrt(2) * (gamma1 - E4*(sqrt(i-2)-sqrt(i-1)));
            HKp(3*N_LL+4+i, N_LL+1+i) = -1/sqrt(2) * (gamma1 - E4*(sqrt(i-2)+sqrt(i-1)));
        end
        if i >= 5
            HKp(N_LL+1+i, 2*N_LL+i) = -1/sqrt(2) * E3 * sqrt(i-3);
            HKp(N_LL+1+i, 3*N_LL+1+i) = 1/sqrt(2) * E3 * sqrt(i-3);
            HKp(2*N_LL+i, N_LL+1+i) = -1/sqrt(2) * E3 * sqrt(i-3);
            HKp(3*N_LL+1+i, N_LL+1+i) = 1/sqrt(2) * E3 * sqrt(i-3);
        end
        
        %        |2N+2+n>=1/sqrt(2)(0,0,|n-3>,|n-2>,0,0);   3<=n<=N_LL+3
        HKp(2*N_LL+2+i,2*N_LL+2+i) = 1/2*(U2+U1) - E0*sqrt(i-2);
        HKp(2*N_LL+2+i,3*N_LL+3+i) = 1/2*(U2-U1);
        HKp(3*N_LL+3+i,2*N_LL+2+i) = 1/2*(U2-U1);
        
        %        |3N+3+n>=1/sqrt(2)(0,0,|n-3>,-|n-2>,0,0);   3<=n<=N_LL+3
        HKp(3*N_LL+3+i,3*N_LL+3+i) = 1/2*(U2+U1) + E0*sqrt(i-2);
        
        %        |4N+5+n>=1/sqrt(2)(0,0,0,0,|n-3>,|n-2>);   3<=n<=N_LL+3
        HKp(4*N_LL+5+i,4*N_LL+5+i) = 1/2*(U6+U5) - E0*sqrt(i-2);
        HKp(4*N_LL+5+i,5*N_LL+6+i) = 1/2*(U6-U5);
        HKp(5*N_LL+6+i,4*N_LL+5+i) = 1/2*(U6-U5);
        %        |4N+5+n>=1/sqrt(2)(0,0,0,0,|n-3>,|n-2>);   3<=n<=N_LL+3
        HKp(5*N_LL+6+i,5*N_LL+6+i) = 1/2*(U6+U5) + E0*sqrt(i-2);
        
        % bilayer-monolayer hyb
        HKp(2*N_LL+2+i,4*N_LL+5+i) = delta1;
        HKp(4*N_LL+5+i,2*N_LL+2+i) = delta1;
        HKp(3*N_LL+3+i,5*N_LL+6+i) = delta1;
        HKp(5*N_LL+6+i,3*N_LL+3+i) = delta1;
    end
    
    %Diagnolize
    if delta1 == 0
        Qp_now = zeros(6*N_LL+9);
        Dp_now = zeros(6*N_LL+9,1);
        [Qp, Dp] = eig(HKp(1:4*N_LL+6,1:4*N_LL+6));
        % sort eigen states to follow the continuous evolution
        Qp_now(1:4*N_LL+6,1:4*N_LL+6) = Qp;
        Dp_now(1:4*N_LL+6,1) = diag(Dp);
        clear Qp Dp

        if B_index >= 2
            for q = 1:4*N_LL+6
                for j = q:4*N_LL+6
                    in_prod = abs(dot(Qp_last(:,q), Qp_now(:,j)));
                    if in_prod > 0.8
                        vec_ = Qp_now(:,j);
                        Qp_now(:,j) = Qp_now(:,q);
                        Qp_now(:,q) = vec_;
                        val_ = Dp_now(j);
                        Dp_now(j) = Dp_now(q);
                        Dp_now(q) = val_;
                        break
                    end
                end
            end
        end
        
        [Qp, Dp] = eig(HKp(4*N_LL+7:6*N_LL+9,4*N_LL+7:6*N_LL+9));
        % sort eigen states to follow the continuous evolution
        Qp_now(4*N_LL+7:6*N_LL+9,4*N_LL+7:6*N_LL+9) = Qp;
        Dp_now(4*N_LL+7:6*N_LL+9,1) = diag(Dp);
        clear Qp Dp
        
        if B_index >= 2
            for q = 4*N_LL+7:6*N_LL+9
                for j = q:6*N_LL+9
                    in_prod = abs(dot(Qp_last(:,q), Qp_now(:,j)));
                    if in_prod > 0.8
                        vec_ = Qp_now(:,j);
                        Qp_now(:,j) = Qp_now(:,q);
                        Qp_now(:,q) = vec_;
                        val_ = Dp_now(j);
                        Dp_now(j) = Dp_now(q);
                        Dp_now(q) = val_;
                        break
                    end
                end
            end
        end
        
        Qp_last = Qp_now;
        LL_level_Kp(B_index,:) = Dp_now;
        clear Q_now D_now
    else
        [Qp_now, D] = eig(HKp);
        % sort eigen states to follow the continuous evolution
        Dp_now = diag(D);
        clear D

        if B_index >= 2
            for q = 1:6*N_LL+9
                for j = q:6*N_LL+9
                    in_prod = abs(dot(Qp_last(:,q), Qp_now(:,j)));
                    if in_prod > 0.9
                        vec_ = Qp_now(:,j);
                        Qp_now(:,j) = Qp_now(:,q);
                        Qp_now(:,q) = vec_;
                        val_ = Dp_now(j);
                        Dp_now(j) = Dp_now(q);
                        Dp_now(q) = val_;
                        break
                    end
                end
            end
        end
        Qp_last = Qp_now;
        LL_level_Kp(B_index,:) = Dp_now;
        clear Qpnow Dp_now
    end

end


dims_m = 2 * N_LL + 3;
dims_b = 4 * N_LL + 6;
LL_K_m = zeros(B_steps, dims_m);
LL_Kp_m = zeros(B_steps, dims_m);
LL_K_b = zeros(B_steps, dims_b);
LL_Kp_b = zeros(B_steps, dims_b);

% 先将第一个磁场先算出来，作为下一步的一个参考标准
B_field = B_fields(1);
[HK_b, HK_m, HKp_b, HKp_m] = construct_HK_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, delta2, B_field, N_LL);
[eigvec_HK_m_last, eigval_HK_m] = eig(HK_m);
[eigvec_HKp_m_last, eigval_HKp_m] = eig(HKp_m);
[eigvec_HK_b_last, eigval_HK_b] = eig(HK_b);
[eigvec_HKp_b_last, eigval_HKp_b] = eig(HKp_b);

eigval_HK_m_diag_last = diag(eigval_HK_m);
eigval_HKp_m_diag_last = diag(eigval_HKp_m);
eigval_HK_b_diag_last = diag(eigval_HK_b);
eigval_HKp_b_diag_last = diag(eigval_HKp_b);

LL_K_m(1,:) = eigval_HK_m_diag_last;
LL_Kp_m(1,:) = eigval_HKp_m_diag_last;
LL_K_b(1,:) = eigval_HK_b_diag_last;
LL_Kp_b(1,:) = eigval_HKp_b_diag_last;


A = HK(4*N_LL+7:6*N_LL+9,4*N_LL+7:6*N_LL+9) - HK_m;

B = HK(:,1);
C = HK(:,2);
HK(:,1) = C;
HK(:,2) = B;
D = HK(1,:);
E = HK(2,:);
HK(1,:) = E;
HK(2,:) = D;

B = HK(:, 3:2*N_LL+4);
C = HK(:,2);
HK(:,2:2*N_LL+3) = B;
HK(:,2*N_LL+4) = C;

D = HK(3:2*N_LL+4, :);
E = HK(2,:);
HK(2:2*N_LL+3,:) = D;
HK(2*N_LL+4,:) = E;

AA = HK(1:4*N_LL+6,1:4*N_LL+6) - HK_b;