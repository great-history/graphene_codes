%% 参数配置
a = 0.142; % C-C bond length (nm)

gamma0 = 2.8;
gamma1 = 0.4; % 是否应该有负号??，根据Slater Koster公式
gamma3 = 0.3;
gamma4 = 0.04;
u = gamma1 * 0.2;  % 0

v0 = 3 * a * gamma0 / 2; % hv = 3 / 2 * a * t
u_bar = u / gamma1;
lambda_bar = gamma3 / gamma0;

num_B_point = 600;
B_field_list = linspace(0, 6, num_B_point);
% lb = 26.5 / sqrt(B_field_list); % magnetic length (nm)

LL_index_trunc = 30; % 
dims = (LL_index_trunc + 1) * 2 - 2;
LL_K_mat = zeros(dims, dims);
LL_Kp_mat = zeros(dims, dims);
LL_K_array = zeros(num_B_point, dims);
LL_Kp_array = zeros(num_B_point, dims);

%% 在每个磁场下构造哈密顿量和对其对角化
for B_index = 1:num_B_point
    B_field = B_field_list(B_index);
    
    e0 = sqrt(2) * v0 * sqrt(B_field) / (gamma1 * 26.5) ;
    
    % construct diagonal block
    n = 1;
    LL_Kp_mat((2*n-1):2*n, (2*n-1):2*n) = [u_bar / 2, 0; 0, u_bar / 2 * (1 - 2 * e0^2)];
    
    for n = 2:LL_index_trunc
        LL_Kp_mat((2*n-1):2*n, (2*n-1):2*n) = [u_bar / 2 * (1 - 2 * n * e0^2), - sqrt(n * (n - 1)) * e0^2; ...
                                            - sqrt(n * (n - 1)) * e0^2,     u_bar / 2 * (1 - 2 * (n - 1) * e0^2)];
    end
    
    % construct off diagonal block
    n = 3;
    LL_Kp_mat((2*(n - 2)-1):2*(n - 2), (2*n-1):2*n) = [0, lambda_bar * e0; 0, 0];
    LL_Kp_mat((2*n-1):2*n, (2*(n - 2)-1):2*(n - 2)) = [0, 0; lambda_bar * e0, 0];
    
    n = 4;
    LL_Kp_mat((2*(n - 3)-1):2*(n - 3), (2*n-1):2*n) = [0, 0; 0, sqrt(2) * lambda_bar * e0];
    LL_Kp_mat((2*n-1):2*n, (2*(n - 3)-1):2*(n - 3)) = [0, 0; 0, sqrt(2) * lambda_bar * e0];
    
    for n = 5:LL_index_trunc
        LL_Kp_mat((2*(n - 3)-1):2*(n - 3), (2*n-1):2*n) = [0, sqrt(n + 1) * lambda_bar * e0; 0, 0];
        LL_Kp_mat((2*n-1):2*n, (2*(n - 3)-1):2*(n - 3)) = [0, 0; sqrt(n + 1) * lambda_bar * e0, 0];
    end
    
    [eigvecs_Kp, eigvals_Kp] = eig(LL_Kp_mat);
    LL_Kp_array(B_index,:) = diag(eigvals_Kp);    
end

figure
axis([B_field_list(1) B_field_list(end) min(min(LL_Kp_array))-0.01 max(max(LL_Kp_array))+0.01])
hold on

plot(B_field_list, LL_Kp_array(:, :), 'm*', 'LineWidth', 0.025, 'MarkerSize',2) 