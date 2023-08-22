%% add magnetic field on zigzag
a = 0.142; % nm
block_length = sqrt(3) * a;
gamma0 = 3.1;
B_field = 1000;

Im = 1j;
length_ratio = (a / 25.66)^2;
phase_coeff = length_ratio * (sqrt(3) / 2);
phase_factor1 = phase_coeff / 4 * B_field;
phase_factor2 = phase_coeff * 3 / 4 * B_field;

Ha = zeros(4,4);

Hb = zeros(4,4);
Hb(1,4) = gamma0;
Hb_adj = Hb';

% 构造block H
N_rows = 200;  % 有11行共有5.5个原胞的厚度
N_u = floor(N_rows / 2);
half_u = mod(N_rows, 2);

dims = 4 * N_u + half_u * 2;
H0 = zeros(dims, dims);
H0(1:4,1:4) = Ha;

for i = 1:(N_u-1)
    x = (i-1)*4;
    y = i*4;
    phase1 = phase_factor1 + (i-1) * phase_coeff;
    phase2 = phase_factor2 + (i-1) * phase_coeff;
    
    Ha(2,1) = gamma0 * exp(-Im * phase1);
	Ha(1,2) = conj(Ha(2,1));
    Ha(3,2) = gamma0;
    Ha(2,3) = conj(Ha(3,2));
    Ha(4,3) = gamma0 * exp(-Im * phase2);
    Ha(3,4) = conj(Ha(4,3));
    
    H0((x+1):(x+4), (x+1):(x+4)) = Ha;
    H0((y+1):(y+4), (x+1):(x+4)) = Hb;
    H0((x+1):(x+4), (y+1):(y+4)) = Hb_adj; 
end

i = N_u;
x = (i-1)*4;
y = i*4;
phase1 = phase_factor1 + (i-1) * phase_coeff;
phase2 = phase_factor2 + (i-1) * phase_coeff;

Ha(2,1) = gamma0 * exp(-Im * phase1);
Ha(1,2) = conj(Ha(2,1));
Ha(3,2) = gamma0;
Ha(2,3) = conj(Ha(3,2));
Ha(4,3) = gamma0 * exp(-Im * phase2);
Ha(3,4) = conj(Ha(4,3));

H0((x+1):(x+4), (x+1):(x+4)) = Ha;

% if half_u == 1
%     y = (N_u-1)*4;
%     z = N_u * 4;
%     H0((z+1):(z+2), (y+1):(y+4)) = Hb(2:3,1:4);
%     H0((y+1):(y+4), (z+1):(z+2)) = Hb_adj(1:4,2:3);
%     H0((z+1):(z+2), (z+1):(z+2)) = Ha(2:3,2:3);
% end

Hc = zeros(4,4);

H1 = zeros(dims, dims);
for i = 1:N_u
    phase1 = phase_factor1 + (i-1) * phase_coeff;
    phase2 = phase_factor2 + (i-1) * phase_coeff;
    
    y = (i-1)*4;
    Hc(1,2) = gamma0 * exp(-Im * phase_factor1);
    Hc(4,3) = gamma0 * exp(-Im * phase_factor2);
    
    H1((y+1):(y+4), (y+1):(y+4)) = Hc;
end

hem = helper_check_hermite(H0, 1e-8);

% 生成k格点(一维)
k_points = 1000;
aks = linspace(0, 2*pi, k_points);
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
%     [eigvecs_now, eigvals_now, eig_num] = helper_re_order_states(eigvecs_last, eigvals_last, eigvecs_now, eigvals_now, dims, eps);
    eig_enes(:,idx) = eigvals_now;
    
%     if (eig_num == dims)
%        eigvecs_last = eigvecs_now;
%        eigvals_last = eigvals_now;
%     end
end

figure
for i = 1:dims
    plot(aks, eig_enes(i,:))
    hold on
end

%% code from Guan JiHuan
N_u = 30;
