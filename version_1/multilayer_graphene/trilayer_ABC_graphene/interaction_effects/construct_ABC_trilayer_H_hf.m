function [H_hf_up, H_hf_down] = construct_ABC_trilayer_H_hf(density_matrix_up, density_matrix_down, exchange_integrals_intralayer_cell, exchange_integrals_interlayer_cell, E_exchange, E_H, E_zeeman)
    % 构造关于ABC trilayer graphene的full Hartree-Fock Hamiltonian
    % 尽管有相互作用，哈密顿量仍然是关于自旋对角的，因此我们可以分为自旋向上的块和自旋向下的块
    % 每个自旋块中基组为{|A1,0>,|A1,1>,|A1,2>,|B3,0>,|B3,1>,|B3,2>}
    
    %% 第一项是近似简并的LL能量(这里LLL的能量都是零，因此不用写这一项)
    %% 第二项是Zeeman term
    H_hf_up = E_zeeman * diag([1, 1, 1, 1, 1, 1]);
    H_hf_down = - E_zeeman * diag([1, 1, 1, 1, 1, 1]);
    
    %% 第三项是层内的交换相互作用 ( E_exchange )
    for n1 = 0:2
        for m2 = 0:2
            val_A1_up = 0.0;
            val_A1_down = 0.0;
            val_B3_up = 0.0;
            val_B3_down = 0.0;
            
            for m1 = 0:2
                for n2 = 0:2
                    num1 = num_system_mechine([n1,n2], 3) + 1;
                    num2 = num_system_mechine([m1,m2], 3) + 1;
                    num_max = max(num1, num2);
                    num_min = min(num1, num2);
                    
                    exchange_integral_intralayer = exchange_integrals_intralayer_cell{num_max}(num_min);
                    
                    % layer A1
                    val_A1_up = val_A1_up + exchange_integral_intralayer * density_matrix_up(m1 + 1, n2 + 1);
                    val_A1_down = val_A1_down + exchange_integral_intralayer * density_matrix_down(m1 + 1, n2 + 1);
                    
                    % layer B3
                    val_B3_up = val_B3_up + exchange_integral_intralayer * density_matrix_up(m1 + 4, n2 + 4);
                    val_B3_down = val_B3_down + exchange_integral_intralayer * density_matrix_down(m1 + 4, n2 + 4);
                end
            end
            
            H_hf_up(n1 + 1, m2 + 1) = H_hf_up(n1 + 1, m2 + 1) - E_exchange * val_A1_up;
            H_hf_down(n1 + 1, m2 + 1) = H_hf_down(n1 + 1, m2 + 1) - E_exchange * val_A1_down;
            H_hf_up(n1 + 4, m2 + 4) = H_hf_up(n1 + 4, m2 + 4) - E_exchange * val_B3_up;
            H_hf_down(n1 + 4, m2 + 4) = H_hf_down(n1 + 4, m2 + 4) - E_exchange * val_B3_down;
        end
    end
    
    %% 第四项是层间的交换相互作用 ( E_exchange )
    for n1 = 0:2
        for m2 = 0:2
            val_A1B3_up = 0.0;
            val_A1B3_down = 0.0;
            val_B3A1_up = 0.0;
            val_B3A1_down = 0.0;
            
            for m1 = 0:2
                for n2 = 0:2
                    num1 = num_system_mechine([n1,n2], 3) + 1;
                    num2 = num_system_mechine([m1,m2], 3) + 1;
                    num_max = max(num1, num2);
                    num_min = min(num1, num2);
                    
                    exchange_integral_interlayer = exchange_integrals_interlayer_cell{num_max}(num_min);
                    
                    % A1产生 B3湮灭
                    val_A1B3_up = val_A1B3_up + exchange_integral_interlayer * density_matrix_up(m1 + 1, n2 + 4);
                    val_A1B3_down = val_A1B3_down + exchange_integral_interlayer * density_matrix_down(m1 + 1, n2 + 4);
                    
                    % B3产生 A1湮灭
                    val_B3A1_up = val_B3A1_up + exchange_integral_interlayer * density_matrix_up(m1 + 4, n2 + 1);
                    val_B3A1_down = val_B3A1_down + exchange_integral_interlayer * density_matrix_down(m1 + 4, n2 + 1);
                end
            end
            
            H_hf_up(n1 + 1, m2 + 4) = H_hf_up(n1 + 1, m2 + 4) - E_exchange * val_A1B3_up;
            H_hf_down(n1 + 1, m2 + 4) = H_hf_down(n1 + 1, m2 + 4) - E_exchange * val_A1B3_down;
            H_hf_up(n1 + 4, m2 + 1) = H_hf_up(n1 + 4, m2 + 1) - E_exchange * val_B3A1_up;
            H_hf_down(n1 + 4, m2 + 1) = H_hf_down(n1 + 4, m2 + 1) - E_exchange * val_B3A1_down;
        end
    end
    
    %% 第五项是层间的capacitive electrostatic energy ( E_H )
    density_layer1 = density_matrix_up(1,1) + density_matrix_up(2,2) + density_matrix_up(3,3);
    density_layer1 = density_layer1 + density_matrix_down(1,1) + density_matrix_down(2,2) + density_matrix_down(3,3);
    density_layer2 = density_matrix_up(4,4) + density_matrix_up(4,4) + density_matrix_up(4,4);
    density_layer2 = density_layer2 + density_matrix_down(4,4) + density_matrix_down(4,4) + density_matrix_down(4,4);
    density_diff = density_layer1 - density_layer2;
    H_hf_up(1:3,1:3) = H_hf_up(1:3,1:3) + E_H * density_diff * diag([1,1,1]);
    H_hf_down(1:3,1:3) = H_hf_down(1:3,1:3) + E_H * density_diff * diag([1,1,1]);
    
    H_hf_up(4:6,4:6) = H_hf_up(4:6,4:6) - E_H * density_diff * diag([1,1,1]);
    H_hf_down(4:6,4:6) = H_hf_down(4:6,4:6) - E_H * density_diff * diag([1,1,1]);
end


%% 构造ABA三层石墨烯
function [H_hf_up, H_hf_down] = construct_ABC_trilayer_H_hf_test(density_matrix_up, density_matrix_down, exchange_integrals_intralayer_cell, exchange_integrals_interlayer_cell, E_exchange, E_H, E_zeeman)
    % 构造关于ABC trilayer graphene的full Hartree-Fock Hamiltonian
    % 尽管有相互作用，哈密顿量仍然是关于自旋对角的，因此我们可以分为自旋向上的块和自旋向下的块
    % 每个自旋块中基组为{|A1,0>,|A1,1>,|A1,2>,|B3,0>,|B3,1>,|B3,2>}
    
    %% 第一项是近似简并的LL能量(这里LLL的能量都是零，因此不用写这一项)
    %% 第二项是Zeeman term
    H_hf_up = E_zeeman * diag([1, 1, 1, 1, 1, 1]);
    H_hf_down = - E_zeeman * diag([1, 1, 1, 1, 1, 1]);
    
    %% 第三项是层内的交换相互作用 ( E_exchange )
    for n1 = 0:2
        for m2 = 0:2
            val_A1_up = 0.0;
            val_A1_down = 0.0;
            val_B3_up = 0.0;
            val_B3_down = 0.0;
            
            for m1 = 0:2
                for n2 = 0:2
                    num1 = num_system_mechine([n1,n2], 3) + 1;
                    num2 = num_system_mechine([m1,m2], 3) + 1;
                    num_max = max(num1, num2);
                    num_min = min(num1, num2);
                    
                    exchange_integral_intralayer = exchange_integrals_intralayer_cell{num_max}(num_min);
                    
                    % layer A1
                    val_A1_up = val_A1_up + exchange_integral_intralayer * density_matrix_up(m1 + 1, n2 + 1);
                    val_A1_down = val_A1_down + exchange_integral_intralayer * density_matrix_down(m1 + 1, n2 + 1);
                    
                    % layer B3
                    val_B3_up = val_B3_up + exchange_integral_intralayer * density_matrix_up(m1 + 4, n2 + 4);
                    val_B3_down = val_B3_down + exchange_integral_intralayer * density_matrix_down(m1 + 4, n2 + 4);
                end
            end
            
            H_hf_up(n1 + 1, m2 + 1) = H_hf_up(n1 + 1, m2 + 1) - E_exchange * val_A1_up;
            H_hf_down(n1 + 1, m2 + 1) = H_hf_down(n1 + 1, m2 + 1) - E_exchange * val_A1_down;
            H_hf_up(n1 + 4, m2 + 4) = H_hf_up(n1 + 4, m2 + 4) - E_exchange * val_B3_up;
            H_hf_down(n1 + 4, m2 + 4) = H_hf_down(n1 + 4, m2 + 4) - E_exchange * val_B3_down;
        end
    end
    
    %% 第四项是层间的交换相互作用 ( E_exchange )
    for n1 = 0:2
        for m2 = 0:2
            val_A1B3_up = 0.0;
            val_A1B3_down = 0.0;
            val_B3A1_up = 0.0;
            val_B3A1_down = 0.0;
            
            for m1 = 0:2
                for n2 = 0:2
                    num1 = num_system_mechine([n1,n2], 3) + 1;
                    num2 = num_system_mechine([m1,m2], 3) + 1;
                    num_max = max(num1, num2);
                    num_min = min(num1, num2);
                    
                    exchange_integral_interlayer = exchange_integrals_interlayer_cell{num_max}(num_min);
                    
                    % A1产生 B3湮灭
                    val_A1B3_up = val_A1B3_up + exchange_integral_interlayer * density_matrix_up(m1 + 1, n2 + 4);
                    val_A1B3_down = val_A1B3_down + exchange_integral_interlayer * density_matrix_down(m1 + 1, n2 + 4);
                    
                    % B3产生 A1湮灭
                    val_B3A1_up = val_B3A1_up + exchange_integral_interlayer * density_matrix_up(m1 + 4, n2 + 1);
                    val_B3A1_down = val_B3A1_down + exchange_integral_interlayer * density_matrix_down(m1 + 4, n2 + 1);
                end
            end
            
            H_hf_up(n1 + 1, m2 + 4) = H_hf_up(n1 + 1, m2 + 4) - E_exchange * val_A1B3_up;
            H_hf_down(n1 + 1, m2 + 4) = H_hf_down(n1 + 1, m2 + 4) - E_exchange * val_A1B3_down;
            H_hf_up(n1 + 4, m2 + 1) = H_hf_up(n1 + 4, m2 + 1) - E_exchange * val_B3A1_up;
            H_hf_down(n1 + 4, m2 + 1) = H_hf_down(n1 + 4, m2 + 1) - E_exchange * val_B3A1_down;
        end
    end
    
    %% 第五项是层间的capacitive electrostatic energy ( E_H )
    density_layer1 = density_matrix_up(1,1) + density_matrix_up(2,2) + density_matrix_up(3,3);
    density_layer1 = density_layer1 + density_matrix_down(1,1) + density_matrix_down(2,2) + density_matrix_down(3,3);
    density_layer2 = density_matrix_up(4,4) + density_matrix_up(4,4) + density_matrix_up(4,4);
    density_layer2 = density_layer2 + density_matrix_down(4,4) + density_matrix_down(4,4) + density_matrix_down(4,4);
    density_diff = density_layer1 - density_layer2;
    H_hf_up(1:3,1:3) = H_hf_up(1:3,1:3) + E_H * density_diff * diag([1,1,1]);
    H_hf_down(1:3,1:3) = H_hf_down(1:3,1:3) + E_H * density_diff * diag([1,1,1]);
    
    H_hf_up(4:6,4:6) = H_hf_up(4:6,4:6) - E_H * density_diff * diag([1,1,1]);
    H_hf_down(4:6,4:6) = H_hf_down(4:6,4:6) - E_H * density_diff * diag([1,1,1]);
end