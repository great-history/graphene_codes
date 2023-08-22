%% 该脚本是对ABC-stacked trilayer graphene的两能带模型的LLL进行Hartree-Fock simulation
% ABC-stacked trilayer graphene的两能带模型与bilayer的两能带模型在数学结构上很相似，这得益于它们相同的手性特征
% 该脚本分为两部分：
% 第一部分：LL能带结构的计算，主要验证LLL是否算错即可
% 第二部分：考虑LLL之间的相互作用(包括Hartree & Fock)

%% 添加路径
addpath('.\interaction_effects\')
addpath('..\utils\')

%% 参数准备
% 基本参数
one_electron_coulomb = 1.602176634 * 10^(-19); %单位是C
epsilon_0 = 8.854187817 * 10^(-12) * 10^(-9); % 单位是F/nm
epsilon_bn = 1; % 无量纲，the in-plane dielectric constant of BN
d_interlayer = 0.335; % 单位是nm

B_field = 20.0;
mag_length = 25.6 / sqrt(B_field); % 以nm为单位

% energy scale
E_zeeman = 0.116 * B_field / 1000; % 单位是eV
E_F = one_electron_coulomb / (4 * pi * epsilon_0 * mag_length); % 单位是 eV fock interaction strength
E_exchange = E_F / epsilon_bn; % 由于hBN存在介电常数（这里是in-plane dielectric constant），这里取为6.6，它可以有效减小交换相互作用的强度
E_H = (2 * d_interlayer) / (2 * mag_length) * E_F; % 单位是 eV hartree interaction strength

%% 计算涉及到的所有exchange integral
% 对于ABC三层石墨烯而言，其LLL具有特殊的性质：每个谷都发生了层极化，并且LL_index n = 0,1,2 ———— 总共有 2 * 3 * 3 = 12个态，即存在十二重简并，这个子空间也被称为duodectet subspace
% 因此构造S矩阵只需考虑n = 0,1,2这三个LL之间的层内或层间的交换相互作用：在PRB 85，165139（2021）TABLE I中已经给出了在任意磁场下的层内交换相互作用和B=20T/30T下的层间交换相互作用
% 下面就是构造出这样的S矩阵
format long
% 构造exchange integral
tic
indice_list = zeros(3*3, 2); % 将(n1,n2)按照三进制进行编码
for ii = 0:2
    num1 = ii * 3;
    for jj = 0:2
        num2 = jj * 1;
        indice = num1 + num2 + 1;
        indice_list(indice, 1) = ii; % 三进制数的第一位(位权为3)
        indice_list(indice, 2) = jj; % 三进制数的第二位(位权为1)
    end
end

exchange_integrals_intralayer_cell = cell(3*3, 1); % 其实只需要三个指标即可，因为 X_{n1 n2, m1 m2}不为零需要满足条件：
exchange_integrals_interlayer_cell = cell(3*3, 1); % 其实只需要三个指标即可，因为 X_{n1 n2, m1 m2}不为零需要满足条件：
% 一种可以减少计算量的方法就是对指标进行编码，然后利用exchange integral的对称性关系(计算量可以减半)
dims = size(indice_list, 1);
count_nonzeros = 0;
for ii = dims:-1:1 %% 对应max([n1 n2]_3, [m1 m2]_3)，其中[n1 n2]_3和[m1 m2]_3是三进制数
    exchange_integrals_intralayer_list = zeros(ii, 1); %% 存放min([n1 n2]_3, [m1 m2]_3)
    exchange_integrals_interlayer_list = zeros(ii, 1); %% 存放min([n1 n2]_3, [m1 m2]_3)
    
    for jj = 1:ii
        n1 = indice_list(ii, 1);
        n2 = indice_list(ii, 2);
        m1 = indice_list(jj, 1);
        m2 = indice_list(jj, 2);
        
        % intralayer (Any B)
        X_nm_intralayer = get_exchange_integral(n1, n2, m1, m2);
        if X_nm_intralayer
           count_nonzeros = count_nonzeros + 1; 
        end
        
        % Interlayer(B = 20T)
        kd_interlayer = 2 * d_interlayer / mag_length; % 计算层间交换相互作用时需要用到
        X_nm_interlayer = get_exchange_integral(n1, n2, m1, m2, kd_interlayer);
        
        % % Interlayer(B = 30T)
        % kd_interlayer = 2 * d_interlayer / 25.6 * sqrt(30); % 乘以2是因为这里计算的是第一层和第三层之间的exchange
        % X_nm_interlayer = get_exchange_integral(n1, n2, m1, m2, d_interlayer);
        
        % 将交换积分存入list之中
        exchange_integrals_intralayer_list(jj) = X_nm_intralayer;
        exchange_integrals_interlayer_list(jj) = X_nm_interlayer;
    end
    
    % 将list存入cell之中
    exchange_integrals_intralayer_cell{ii} = exchange_integrals_intralayer_list; %% 存放max([n1 n2]_3, [m1 m2]_3)
    exchange_integrals_interlayer_cell{ii} = exchange_integrals_interlayer_list; %% 存放max([n1 n2]_3, [m1 m2]_3)
end
toc

%% 根据exchange integral计算出S矩阵
% 在这里不需要这一步 。。。。。。。。。。

%% 构造一个随机矩阵:Initialization
dim_subspace = 6; % ABC trilayer的LLL的LL_index可以是0,1,2，两个谷都是这样，因此共有2*3=6个
filling_factor = 1; % (整数)填充数
% 构造方法四（最科学） ： 在构造方法一的基础上添加一个权重(这样的话容易出现局域极小值的解)，权重随机生成
filling_list = randi([0 dim_subspace], 2, 1);
while filling_list == [0;0]
    filling_list = randi([0 dim_subspace], 2, 1);
end

% density_matrix_down_temp = construct_random_density_matrix(filling_list(1), dim_subspace, "complex");
density_matrix_down_temp = construct_random_density_matrix(filling_list(1), dim_subspace, "complex");
% density_matrix_up_temp = construct_random_density_matrix(filling_list(2), dim_subspace, "complex");
density_matrix_up_temp = zeros(6);

weight1 = unifrnd(0,1);
weight_list = [weight1; 1 - weight1];
density_matrix_down_temp = weight_list(1) * density_matrix_down_temp;
density_matrix_up_temp = weight_list(2) * density_matrix_up_temp;
disp("自旋向上密度矩阵为：")
density_matrix_up_temp
disp("自旋向下密度矩阵为：")
density_matrix_down_temp

%% Hartree-Fock 初始化
% % 只能是哈密顿量的构造存在问题
[H_hf_up, H_hf_down] = construct_ABC_trilayer_H_hf(density_matrix_up_temp, density_matrix_down_temp, ...
                                                   exchange_integrals_intralayer_cell, exchange_integrals_interlayer_cell, E_exchange, E_H, E_zeeman);
if ~helper_check_hermite(H_hf_down, 1e-3)
    disp("H_hf_down不是厄米矩阵")
end
H_hf_up = refine_H_hf(H_hf_up, dim_subspace);
H_hf_down = refine_H_hf(H_hf_down, dim_subspace);
[error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, filling_factor, dim_subspace);
[eigvecs_down_last, eigvals_down_last] = eig(H_hf_down);
[eigvecs_up_last, eigvals_up_last] = eig(H_hf_up);
eigvals_down_last = diag(eigvals_down_last);
eigvals_up_last = diag(eigvals_up_last);
disp(["误差为", error])

%% start self-consistent calculation
tic

steps = 0;
while error >= 1e-8
    [H_hf_up, H_hf_down] = construct_ABC_trilayer_H_hf(density_matrix_up_temp, density_matrix_down_temp, ...
                                                       exchange_integrals_intralayer_cell, exchange_integrals_interlayer_cell, E_exchange, E_H, E_zeeman);
    H_hf_up = refine_H_hf(H_hf_up, 6);
    H_hf_down = refine_H_hf(H_hf_down, 6);
    [error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, filling_factor, 6);
    steps = steps + 1;
    
    [eigvecs_down_temp, eigvals_down_temp] = eig(H_hf_down);
    [eigvecs_up_temp, eigvals_up_temp] = eig(H_hf_up);
    eigvals_down_temp = diag(eigvals_down_temp);
    eigvals_up_temp = diag(eigvals_up_temp);
    
    error_down = max(abs(eigvals_down_temp - eigvals_down_last));
    error_up = max(abs(eigvals_up_temp - eigvals_up_last));
    error = max(error_down, error_up);
    eigvals_down_last = eigvals_down_temp;
    eigvals_up_last = eigvals_up_temp;
    
    if ~mod(steps, 10000)
        fprintf("进行到第%d步，误差为%10.9f\n", steps, error)
        % disp(["进行到第", steps, "步，误差为", error]) % 当收敛比较慢时，可以通过这个来查看收敛的速度
    end
end

if steps == 0
    steps = 10;
end

for ii = 1:steps*100
    [H_hf_up, H_hf_down] = construct_ABC_trilayer_H_hf(density_matrix_up_temp, density_matrix_down_temp, ...
                                                       exchange_integrals_intralayer_cell, exchange_integrals_interlayer_cell, E_exchange, E_H, E_zeeman);
    if ~helper_check_hermite(H_hf_down, 1e-3)
        disp("H_hf_down不是厄米矩阵")
    end                                             
    H_hf_up = refine_H_hf(H_hf_up, dim_subspace);
    H_hf_down = refine_H_hf(H_hf_down, dim_subspace);
    [error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, filling_factor, 6);
    if ~mod(ii,100)
        fprintf("进行到第%d步，误差为%10.9f\n", steps, error)
        % disp(["误差为", error])
        % disp("自旋向上密度矩阵为：")
        % density_matrix_up_temp
        % disp("自旋向下密度矩阵为：")
        % density_matrix_down_temp
        
        % [eigvecs_temp, eigvals_temp] = eig(H_hf);
        % disp(["本征值为：", eigvals_temp(1,1), eigvals_temp(2,2), eigvals_temp(3,3)])
    end
end

density_matrix_down_list = zeros(6, 6, 10);
eigvals_down_list = zeros(6, 10);
for ii = 1:10
    [H_hf_up, H_hf_down] = construct_ABC_trilayer_H_hf(density_matrix_up_temp, density_matrix_down_temp, ...
                                                       exchange_integrals_intralayer_cell, exchange_integrals_interlayer_cell, E_exchange, E_H, E_zeeman);
    if ~helper_check_hermite(H_hf_down, 1e-3)
        disp("H_hf_down不是厄米矩阵")
    end                                               
    H_hf_up = refine_H_hf(H_hf_up, dim_subspace);
    H_hf_down = refine_H_hf(H_hf_down, dim_subspace);
    [error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_temp, density_matrix_down_temp, filling_factor, 6);

    disp(["误差为", error])
    % disp("自旋向上密度矩阵为：")
    % density_matrix_up_temp
    % disp("自旋向下密度矩阵为：")
    % density_matrix_down_temp
    
    [eigvecs_down_temp, eigvals_down_temp] = eig(H_hf_down);
    density_matrix_down_list(:, :, ii) = density_matrix_down_temp;
    eigvals_down_list(:, ii) = real(diag(eigvals_down_temp));
    % [eigvecs_temp, eigvals_temp] = eig(H_hf_down);
    % disp(["本征值为：", eigvals_temp(1,1), eigvals_temp(2,2), eigvals_temp(3,3)])
end

disp("自旋向上密度矩阵为：")
density_matrix_up_temp
disp("自旋向下密度矩阵为：")
density_matrix_down_temp

[eigvecs_down_temp, eigvals_down_temp] = eig(H_hf_down);
[eigvecs_up_temp, eigvals_up_temp] = eig(H_hf_up);
% 
% disp("自旋向上本征态：")
% real(diag(eigvals_up_temp))

disp("自旋向下本征态：")
eigvals_down_list(:,1) / (E_F * sqrt(pi / 2))
% sum(eigvals_down_list(:,1) / (E_F * sqrt(pi / 2))) % -8.430850167236729
% density_matrix_down_temp(1,1) + density_matrix_down_temp(2,2) + density_matrix_down_temp(3,3) + density_matrix_down_temp(4,4) + density_matrix_down_temp(5,5) + density_matrix_down_temp(6,6)
% density_matrix_up_temp(1,1) + density_matrix_up_temp(2,2) + density_matrix_up_temp(3,3) + density_matrix_up_temp(4,4) + density_matrix_up_temp(5,5) + density_matrix_up_temp(6,6)
toc