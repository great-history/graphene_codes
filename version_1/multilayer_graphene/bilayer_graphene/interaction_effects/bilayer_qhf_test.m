% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 创建优化问题变量
% r0 = optimvar('r0');
% r1 = optimvar('r1');
% r2 = optimvar('r2');
% phi1 = optimvar('phi1');
% phi2 = optimvar('phi2');
% 
% % 创建使用优化变量的表达式作为目标函数
% obj = objfunx_exchange_energy(r0, r1, r2, phi1, phi2);
% 
% % 创建一个以 obj 为目标函数的优化问题
% prob = optimproblem('Objective', obj, 'ObjectiveSense','min');
% normlization = r0.^2 + r1.^2 + r2.^2 == 1;
% prob.Constraints.cons1 = normlization;
% phi1_cons1 = (phi1 >= - pi);
% phi1_cons2 = (phi1 <= pi);
% prob.Constraints.cons2 = phi1_cons1;
% prob.Constraints.cons3 = phi1_cons2;
% phi2_cons1 = (phi2 >= - pi);
% phi2_cons2 = (phi2 <= pi);
% prob.Constraints.cons4 = phi2_cons1;
% prob.Constraints.cons5 = phi2_cons2;
% 
% prob.Constraints.cons6 = (r0 >= 0);
% prob.Constraints.cons7 = (r1 >= 0);
% prob.Constraints.cons8 = (r2 >= 0);
% 
% % phi2_cons = (phi2 >= - pi) && (phi2 <= pi);
% % prob.Constraints.constr = phi2_cons;
% 
% % 创建一个结构体，将初始点表示为
% parameters.r0 = 0.3;
% parameters.r1 = 0.4;
% parameters.r2 = sqrt(1 - 0.3^2 - 0.4^2);
% parameters.phi1 = pi;
% parameters.phi2 = 0;
% 
% show(prob)
% [sol,fval] = solve(prob, parameters);
% 
% % 绘制exchange_energy关于角度(phi1, phi2)的一个分布
% f = @objfunx_exchange_energy_ang;
% rnge = [-pi pi 0 2 * pi];
% fcontour(f, rnge)
% colorbar
% hold off
% % disp(objfunx_exchange_energy(sol.r0, sol.r1, sol.r2, sol.phi1, sol.phi2))

%% 添加路径
addpath('.\interaction_effects\')
format long
%% test
S0000 = 0.64;
S0011 = 0.47;
S0212 = 0.29;
S0102 = -0.29;

S_tensor = zeros(3,3,3,3);
for alpha = 1:3
    for beta = 1:3
        for lambda = 1:3
            for sigma = 1:3
                % 非零必须至少有两个指标一样
                if (alpha == beta) && (alpha == lambda) && (alpha == sigma)
                    S_tensor(alpha, beta, lambda, sigma) = S0000;
                    continue % 一定要加continue!!!
                end
                
                if (alpha == beta) && (lambda == sigma)
                    S_tensor(alpha, beta, lambda, sigma) = S0011;
                    continue % 一定要加continue!!!
                end
                
                if (alpha == sigma) && (lambda == beta)
                    S_tensor(alpha, beta, lambda, sigma) = S0011;
                    continue % 一定要加continue!!!
                end
                
                if (beta == sigma)
                    if ~(alpha == lambda) && ~(alpha == beta) && ~(lambda == beta)
                        if beta == 3
                            S_tensor(alpha, beta, lambda, sigma) = S0212;
                        else
                            S_tensor(alpha, beta, lambda, sigma) = S0102;
                        end
                    end
                    
                    continue
                end
                
                if (alpha == lambda)
                    if ~(beta == sigma) && ~(beta == alpha) && ~(sigma == alpha)
                        if alpha == 3
                            S_tensor(alpha, beta, lambda, sigma) = S0212;
                        else
                            S_tensor(alpha, beta, lambda, sigma) = S0102;
                        end
                    end
                    
                    continue
                end
            end
        end
    end
end

%% energy scale
B_field = 1;
E_zeeman = 0.289 * B_field;
E_ee = 140 * sqrt(B_field) / 6;  % 6为相对介电常数
filling_factor = 4;

%% 构造一个初始的随机矩阵
% % 构造方法零
% if sol.r0 >=0 
%     eigvec_temp = [sol.r0; sol.r1 * exp(1j * sol.phi1); sol.r2 * exp(1j * sol.phi2)];
% else
%     eigvec_temp = (-1) * [sol.r0; sol.r1 * exp(1j * sol.phi1); sol.r2 * exp(1j * sol.phi2)];
% end
% density_matrix_temp = eigvec_temp * eigvec_temp';

% % 构造方法一
% density_matrix_down_temp = construct_random_density_matrix(1, 3);
% density_matrix_up_temp = construct_random_density_matrix(1, 3);

% % 构造方法二
% density_matrix_down_temp = construct_random_density_matrix(1, 3);
% density_matrix_up_temp = zeros(3);
% 
% % 构造方法三 ： 在构造方法一的基础上添加一个权重(这样的话容易出现局域极小值的解)，权重手动调节
% density_matrix_down_temp = construct_random_density_matrix(1, 3);
% density_matrix_up_temp = construct_random_density_matrix(1, 3);
% 
% weight1 = unifrnd(0,1);
% weight_list = [weight1; 1 - weight1];
% density_matrix_down_temp = 1.25 * density_matrix_down_temp;
% density_matrix_up_temp = 0.25 * density_matrix_up_temp;

% 构造方法四（最科学） ： 在构造方法一的基础上添加一个权重(这样的话容易出现局域极小值的解)，权重随机生成
filling_list = randi([0 3], 2, 1);
while filling_list == [0;0]
    filling_list = randi([0 3], 2, 1);
end
density_matrix_down_temp = construct_random_density_matrix(filling_list(1), 3, "complex");
density_matrix_up_temp = construct_random_density_matrix(filling_list(2), 3, "complex");

weight1 = unifrnd(0,1);
weight_list = [weight1; 1 - weight1];
density_matrix_down_temp = weight_list(1) * density_matrix_down_temp;
density_matrix_up_temp = weight_list(2) * density_matrix_up_temp;

%% Hartree-Fock 初始化（第一步）
[H_hf_up, H_hf_down] = construct_H_hf(density_matrix_up_temp, density_matrix_down_temp, S_tensor, E_ee, E_zeeman);
[error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, filling_factor, 3);
% disp(["误差为", error])
% disp("密度矩阵为：")
% density_matrix_temp

%% start self-consistent calculation
tic
steps = 0;
while error >= 1e-8
    [H_hf_up, H_hf_down] = construct_H_hf(density_matrix_up_temp, density_matrix_down_temp, S_tensor, E_ee, E_zeeman);
    H_hf_up = refine_H_hf(H_hf_up, 3);
    H_hf_down = refine_H_hf(H_hf_down, 3);
    [error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, filling_factor, 3);
    steps = steps + 1;
    % disp(["误差为", error])
end

if steps < 2
    steps = 10;
end

for ii = 1:steps*100
    [H_hf_up, H_hf_down] = construct_H_hf(density_matrix_up_temp, density_matrix_down_temp, S_tensor, E_ee, E_zeeman);
    H_hf_up = refine_H_hf(H_hf_up, 3);
    H_hf_down = refine_H_hf(H_hf_down, 3);
    [error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, filling_factor, 3);
    if ~mod(ii,100)
        disp(["误差为", error])
        % [eigvecs_temp, eigvals_temp] = eig(H_hf);
        % disp(["本征值为：", eigvals_temp(1,1), eigvals_temp(2,2), eigvals_temp(3,3)])
        disp("自旋向上密度矩阵为：")
        density_matrix_up_temp
        disp("自旋向下密度矩阵为：")
        density_matrix_down_temp
    end
end

[eigvecs_down_temp, eigvals_down_temp] = eig(H_hf_down);
[eigvecs_up_temp, eigvals_up_temp] = eig(H_hf_up);

disp("自旋向上本征态：")
real(diag(eigvals_up_temp))

disp("自旋向下本征态：")
real(diag(eigvals_down_temp))
toc


function ex_ene = objfunx_exchange_energy(r0, r1, r2, phi1, phi2)
    % S_tensor
    % fock energy functional @ B = 1T
    S0000 = 0.64;
    S0011 = 0.47;
    S0212 = 0.29;
    S0102 = -0.29;
    ex_ene = S0000 * (r0^4 + r1^4 + r2^4) + 4 * S0011 * (r0^2 * r1^2 + r0^2 * r2^2 + r1^2 * r2^2);
    ex_ene = ex_ene + 4 * S0212 * r0 * r1 * r2^2 * cos(2 * phi2 - phi1);
    ex_ene = ex_ene + 4 * S0102 * r0^2 * r1 * r2 * cos(phi1 + phi2);
    ex_ene = ex_ene + 4 * S0102 * r0 * r1^2 * r2 * cos(2 * phi1 - phi2);
    ex_ene = - ex_ene;
    % f = exp(x).*(4*x.^2 + 2*y.^2 + 4*x.*y + 2*y - 1);
end

function ex_ene = objfunx_exchange_energy_ang(phi1, phi2)
    r0 = 1 / sqrt(3);
    r1 = 1 / sqrt(3); 
    r2 = 1 / sqrt(3);
    ex_ene = objfunx_exchange_energy(r0, r1, r2, phi1, phi2);
end