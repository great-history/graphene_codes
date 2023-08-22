% syms x
% fplot(laguerreL(1:4,x))
% axis([-2 10 -10 10])
% grid on
% 
% ylabel('L_n(x)')
% title('Laguerre polynomials of orders 1 through 4')
% legend('1','2','3','4','Location','best')


d = 0.142;
a = sqrt(3) * d;
theta = 1.05;

k_d = 4 * pi / (3 * a);
k_theta = 2 * k_d * sin(theta / 2);

qb = k_theta * [0, -1];
qr = k_theta * [sqrt(3) / 2, 1 / 2];
ql = k_theta * [-sqrt(3) / 2, 1 / 2];

%% commensurate case
p = 1;
q = 10;
vel = 3 * 10^(6);
b = sqrt(3) * (k_theta)^2 * p / q / (8 * pi);
mag_len = 1 / sqrt(b);

z_sqaure = (k_theta * mag_len / sqrt(2))^2;
zb = qb * mag_len / sqrt(2);
zl = ql * mag_len / sqrt(2);
zr = qr * mag_len / sqrt(2);

%% LL's form factor
n = 3;
m = 1;

syms x
fplot(laguerreL(n,m,x))
axis([-2 10 -10 10])
grid on

syms x
fplot(laguerreL(m,n-m,x))
axis([-2 10 -10 10])
grid on

format long
fnm_b = get_LL_form_factor(n, m, zb(1), zb(2));
fnm_r = get_LL_form_factor(n, m, zr(1), zr(2));
fnm_l = get_LL_form_factor(n, m, zl(1), zl(2));
    
function fnm = get_LL_form_factor(n, m, zx, zy)
    z_square = (zx)^2 + (zy)^2;
    fnm = sqrt(factorial(m) / factorial(n)) * (-zx + 1j*zy)^(n - m) * ...   
            exp(- z_square / 2) * laguerreL(m, n - m, z_square);
end

