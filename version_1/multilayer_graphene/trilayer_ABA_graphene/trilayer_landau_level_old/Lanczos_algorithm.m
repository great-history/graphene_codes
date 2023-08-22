a = 0.246; % nm
block_length = sqrt(3) * a;
gamma0 = 3.1; % eV
mass = 0;

t = 1;
B_field = 150;

% 加入磁场之后，我们可以找出一个面积S, 其对应的磁通\Phi = BS与量子磁通\Phi_0 = BS的比值\Phi / \Phi_0
length_ratio = (a / 25.66)^2;
N_rows = 4; % 有11行共有5.5个原胞的厚度
dims = N_rows;
Im = 1j;
y_shift = 0;
BS_BS0 = 2 * pi * length_ratio * B_field * sqrt(3) / (8 * pi);

[H0, H1, H_1, y_coords] = construct_monolayer_GNR_zigzag(t, [0.0,0.0], N_rows, y_shift, BS_BS0);

% 直接用eig对H0进行求解
[eigvecs_last, eigvals_last] = eig(H0);

% 使用Lanczos test:对H0进行求解
v0 = randn(N_rows, 1);
norm2_v =(v0'*v0);
% v0 = v0 / norm_v;

n_iter = 4;
diag_elements = zeros(n_iter, 1);
off_diag_elements = zeros(n_iter - 1, 1);

a0 = v0' * H0 * v0 / norm2_v;
b0 = 0;
diag_elements(1) = a0;

vn_last = v0;
vn_last_last = 0;
an_last = a0;
bn_last = b0;
norm2_vn_last = norm2_v;

for i = 2:n_iter
    vn = H0 * vn_last - an_last * vn_last - bn_last * vn_last_last;
    norm2_vn =(vn'*vn);
    
    an = vn' * H0 * vn / norm2_vn;
    bn_2 = norm2_vn / norm2_vn_last;
    bn = sqrt(bn_2);
    
    diag_elements(i) = an;
    off_diag_elements(i-1) = bn;
    
    vn_last = vn;
    vn_last_last = vn_last;
    an_last = an;
    bn_last = bn;
    norm2_vn_last = norm2_vn;
end

% construct tri-diagonal matrix T
Tri_diag = diag(diag_elements) + diag(off_diag_elements, -1) + diag(off_diag_elements, 1);
Tri_diag = real(Tri_diag);
[eigvecs_now, eigvals_now] = eig(Tri_diag);