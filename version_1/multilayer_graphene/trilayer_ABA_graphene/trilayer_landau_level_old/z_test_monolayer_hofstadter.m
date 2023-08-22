%% add magnetic field on zigzag
a = 0.246; % nm
block_length = sqrt(3) * a;
gamma0 = 3.1; % eV

t = 1;
B_field = 0;

% 加入磁场之后，我们可以找出一个面积S, 其对应的磁通\Phi = BS与量子磁通\Phi_0 = BS的比值\Phi / \Phi_0
length_ratio = (a / 25.66)^2;
% BS_BS0 = 2 * pi * length_ratio * B_field * sqrt(3) / (8 * pi);  % B_field = 10T时，如果宽度只有100个原胞的话，朗道能级是看不出来的
% 如果要看到朗道能级，要么就是宽度要足够大，要么就是磁场足够高

N_rows = 200;  % 有11行共有5.5个原胞的厚度
BS_BS0 = 0;
N_block = floor(N_rows / 4);
half_block = mod(N_rows, 4);

dims = N_rows;
H0 = zeros(dims, dims);
H1 = zeros(dims, dims);

% construct the tight-binding model under a magnetic field
mass = 0;

Im = 1j;
for n = 1:N_block
    phi1 = BS_BS0 * (2 * n - 1);  % 12
    phi2 = BS_BS0 * (2 * n - 2);  % 34 或 43
    
    H0(4*(n-1)+4, 4*(n-1)+4) = mass;
    H0(4*(n-1)+3, 4*(n-1)+3) = -mass;
    H0(4*(n-1)+2, 4*(n-1)+2) = mass;
    H0(4*(n-1)+1, 4*(n-1)+1) = -mass;
    
    H0(4*(n-1)+3, 4*(n-1)+4) = t * exp(-Im * phi1);
    H0(4*(n-1)+4, 4*(n-1)+3) = conj(H0(4*(n-1)+3, 4*(n-1)+4));
    
    H0(4*(n-1)+2, 4*(n-1)+3) = t;
    H0(4*(n-1)+3, 4*(n-1)+2) = conj(H0(4*(n-1)+2, 4*(n-1)+3));
    
    H0(4*(n-1)+1, 4*(n-1)+2) = t * exp(Im * phi2);
    H0(4*(n-1)+2, 4*(n-1)+1) = conj(H0(4*(n-1)+1, 4*(n-1)+2));
    
    if ~(n == N_block)
        H0(4*n + 1, 4*(n-1) + 4) = t;
        H0(4*(n-1) + 4, 4*n + 1) = conj(H0(4*n + 1, 4*(n-1)+4));
%     else
%         H0(4*(n-1) + 4, 1) = t;
%         H0(1, 4*(n-1) + 4) = t;
    end
    
    % H1不是厄密的矩阵
    H1(4*(n-1)+4, 4*(n-1)+3) = t * exp(-Im * phi1); 
    H1(4*(n-1)+1, 4*(n-1)+2) = t * exp(-Im * phi2);
    
    % H_{-1}是H1的厄密，故不用创建
end

if ~(half_block == 0)
    phi1 = BS_BS0 * (2 * N_block + 1);  % 12
    phi2 = BS_BS0 * (2 * N_block);  % 34 或 43
end


H_1 = H1'; % H_{-1}是H1的厄密


% hem = helper_check_hermite(H0, 1e-8);


% 生成k格点(一维)
k_points = 5000;
aks = linspace(-pi, pi, k_points);
eig_enes = zeros(dims, k_points);

idx = 1;
Hk = H0 + H1 * exp(Im * aks(idx)) + H_1 * exp(-Im * aks(idx));
[eigvecs_last, eigvals] = eig(Hk);
eigvals_last = diag(eigvals);
eig_enes(:,idx) = eigvals_last;

eps = 0.1;
for idx = 2:k_points
    Hk = H0 + H1 * exp(Im * aks(idx)) + H_1 * exp(-Im * aks(idx));
    [eigvecs_now, eigvals] = eig(Hk);
    eigvals_now = diag(eigvals);

    eig_enes(:,idx) = eigvals_now;
end

figure
axis([aks(1) aks(end) -1 1])
hold on
for i = 1:dims
    plot(aks, eig_enes(i,:))
    hold on
end