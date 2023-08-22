%% armchair monolayer graphene
% 构造Ha, Hb, Hc这三个矩阵
a = 0.142; % nm
block_length = sqrt(3) * a;
gamma0 = 3.1;
Ha = zeros(4,4);
for i = 2:4
    Ha(i-1,i) = gamma0;
    Ha(i,i-1) = gamma0;
end

Hb = zeros(4,4);
Hb(2,1) = gamma0;
Hb(1,2) = gamma0;
Hb_adj = Hb';

Hc = zeros(4,4);
Hc(4,1) = gamma0;

% 构造block H
N_rows = 100;  % 有11行共有5.5个原胞的厚度
N_u = floor(N_rows / 2);
half_u = mod(N_rows, 2);

dims = 4 * N_u + half_u * 2;
H0 = zeros(dims, dims);
H0(1:4,1:4) = Ha;
for i = 2:N_u
    x = (i-2)*4;
    y = (i-1)*4;
    H0((x+1):(x+4), (y+1):(y+4)) = Hb_adj;
    H0((y+1):(y+4), (y+1):(y+4)) = Ha;
    H0((y+1):(y+4), (x+1):(x+4)) = Hb;
end

if half_u == 1
    y = (N_u-1)*4;
    z = N_u * 4;
    H0((z+1):(z+2), (y+1):(y+4)) = Hb(2:3,1:4);
    H0((y+1):(y+4), (z+1):(z+2)) = Hb_adj(1:4,2:3);
    H0((z+1):(z+2), (z+1):(z+2)) = Ha(2:3,2:3);
end

H1 = zeros(dims, dims);
for i = 1:N_u
    y = (i-1)*4;
    H1((y+1):(y+4), (y+1):(y+4)) = Hc;
end

hem = helper_check_hermite(H0, 1e-8);

% 生成k格点(一维)
k_points = 1000;
aks = linspace(-pi, pi, k_points);
eig_enes = zeros(dims, k_points);

Im = 1j;

idx = 1;
Hk = H0 + H1 * exp(Im * aks(idx)) + H1' * exp(-Im * aks(idx));
[eigvecs_last, eigvals] = eig(Hk);
eigvals_last = diag(eigvals);
eig_enes(:,idx) = eigvals_last;

eps = 0.01 * gamma0;
for idx = 1:k_points
    Hk = H0 + H1 * exp(Im * aks(idx)) + H1' * exp(-Im * aks(idx));
    [eigvecs_now, eigvals] = eig(Hk);
    eigvals_now = diag(eigvals);
    [eigvecs_now, eigvals_now, eig_num] = helper_re_order_states(eigvecs_last, eigvals_last, eigvecs_now, eigvals_now, dims, eps);
    eig_enes(:,idx) = eigvals_now;
    
    if (eig_num == dims)
       eigvecs_last = eigvecs_now;
       eigvals_last = eigvals_now;
    end
end

figure
for i = 1:dims
    plot(aks, eig_enes(i,:))
    hold on
end

%% zigzag monolayer graphene
a = 0.142; % nm
block_length = a;
gamma0 = 3.1;

% 构造Ha, Hb, Hc这三个矩阵
Ha = zeros(4,4);
for i = 2:4
    Ha(i-1,i) = gamma0;
    Ha(i,i-1) = gamma0;
end

Hb = zeros(4,4);
Hb(1,4) = gamma0;
Hb_adj = Hb';

Hc = zeros(4,4);
Hc(1,2) = gamma0;
Hc(4,3) = gamma0;

% 构造block H
N_rows = 40;  % 有11行共有1.75个原胞的厚度
N_u = floor(N_rows / 4);
half_u = mod(N_rows, 4);

dims = 4 * N_u + half_u;
H0 = zeros(dims, dims);
H0(1:4,1:4) = Ha;

H0(1:4,1:4) = Ha;
for i = 2:N_u
    x = (i-2)*4;
    y = (i-1)*4;
    H0((y+1):(y+4), (y+1):(y+4)) = Ha;
    H0((x+1):(x+4), (y+1):(y+4)) = Hb_adj;
    H0((y+1):(y+4), (x+1):(x+4)) = Hb;
end

if ~(half_u == 0)
    y = (N_u-1)*4;
    z = N_u * 4;
    H0((z+1):(z+half_u), (z+1):(z+half_u)) = Ha(1:half_u, 1:half_u);
    
    H0((z+1):(z+half_u), (y+1):(y+4)) = Hb(1:half_u, 1:4);
    H0((y+1):(y+4), (z+1):(z+half_u)) = Hb_adj(1:4, 1:half_u);
end


H1 = zeros(dims, dims);
for i = 1:N_u
    y = (i-1)*4;
    H1((y+1):(y+4), (y+1):(y+4)) = Hc;
end

if ~(half_u == 0)
    z = N_u * 4;
    H1((z+1):(z+half_u), (z+1):(z+half_u)) = Hc(1:half_u, 1:half_u);
end
hem = helper_check_hermite(H0, 1e-8);

% 生成k格点(一维)
k_points = 1000;
aks = linspace(-pi, pi, k_points);
eig_enes = zeros(dims, k_points);

Im = 1j;
for idx = 1:k_points
    Hk = H0 + H1 * exp(Im * aks(idx)) + H1' * exp(-Im * aks(idx));
    [eigvec_Hk_now, eigval_Hk_now] = eig(Hk);
    eig_enes(:, idx) = diag(eigval_Hk_now);
end

figure
for i = 1:dims
    plot(aks, eig_enes(i,:))
    hold on
end