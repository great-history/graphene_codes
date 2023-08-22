%% 
function [H_hf_up, H_hf_down] = construct_ABA_trilayer_H_hf_LLL(density_matrix_up, density_matrix_down,...
                                                                exchange_integrals_intralayer_cell, exchange_integrals_interlayer_cell, ...
                                                                E_LL_list, E_exchange, E_H, E_zeeman)
    % 构造关于ABA trilayer graphene的full Hartree-Fock Hamiltonian
    % 尽管有相互作用，哈密顿量仍然是关于自旋对角的，因此我们可以分为自旋向上的块和自旋向下的块
    % 每个自旋块中基组为{|+,2,0>,         |+,2,1>,         |+,1,0>,         |-,1,0>,         |-,2,0>,   |-,2,1>}
    % 分别是            |A1,0> + |A3,0>, |A1,1> + |A3,1>, |A1,0> - |A3,0>, |B1,0> - |B3,0>, |B2,0>,    |B2,1>
    % 公式来自于  PRB 85，165139 ( 2012 )
    
    H_hf_up = zeros(6);
    H_hf_down = zeros(6);
    %% 第一项是近似简并的LL能量(这里LLL的能量都是零，因此不用写这一项)
    % E_LL_list中是{|+,2,0>, |+,2,1>, |+,1,0>, |-,1,0>, |-,2,0>, |-,2,1>}对应的能量
    % 在最简化的模型中E_LL_list为[gamma2 / 2, gamma2 / 2, -gamma2 / 2, delta' - gamma5 / 2, 0, 0]
    H_hf_up = H_hf_up + diag(E_LL_list);
    H_hf_down = H_hf_down + diag(E_LL_list);
    
    %% 第二项是Zeeman term
    H_hf_up = H_hf_up + E_zeeman * diag([1, 1, 1, 1, 1, 1]);
    H_hf_down = H_hf_down - E_zeeman * diag([1, 1, 1, 1, 1, 1]);
    
    %% 第三项是层内的交换相互作用 ( E_exchange )
    % 基组是 【|A1,0> + |A3,0>】, 【|A1,1> + |A3,1>】, 【|A1,0> - |A3,0>】, 【|B1,0> - |B3,0>】, |B2,0>, |B2,1>
    % 首先是【|A1,0> + |A3,0>】, 【|A1,1> + |A3,1>】
    
    
    % 然后是【|A1,0> - |A3,0>】, 【|B1,0> - |B3,0>】
    
    
    % 接着是|B2,0>, |B2,1>
    
    
    % 之后是【|A1,0> + |A3,0>】, 【|A1,1> + |A3,1>】与【|A1,0> - |A3,0>】, 【|B1,0> - |B3,0>】
    
    
    % 然后又是【|A1,0> + |A3,0>】, 【|A1,1> + |A3,1>】与|B2,0>, |B2,1>
    
    
    % 最后是【|A1,0> - |A3,0>】, 【|B1,0> - |B3,0>】与|B2,0>, |B2,1>
    
    
    %% 第五项是层间的capacitive electrostatic energy ( E_H )
    density_layer2 = density_matrix_up(5, 5) + density_matrix_up(6, 6) + density_matrix_down(5, 5) + density_matrix_down(6, 6);
    H_hf_up = H_hf_up - E_H * density_layer2 * diag([1, 1, 1, 1, -1, -1]); % 这里的E_H = d / mag_length * E_F
    H_hf_down = H_hf_down - E_H * density_layer2 * diag([1, 1, 1, 1, -1, -1]); % 这里的E_H = d / mag_length * E_F
    
end