function [carrier_density_matrix] = helper_get_n_asfo_ED(DOS_ED_matrix, delta_ene, E_points, Delta1_steps)
    % 要得到载流子浓度作为 能量E 和 Delta1 的函数，首先我们要找到CNP，要确定CNP，我们要先确定CNP上下最近邻的两条LL，对于一般的石墨烯体系，这应该是确定的两条
    % 比如monolayer就是LLm0_K和LLm0_Kp这两条，而trilayer在Delta1 = 0时则是LLb0_K和LLb1_K这两条
    
    carrier_density_matrix = zeros(E_points, Delta1_steps);  % 载流子浓度作为能量E和磁场B的函数
    
    for Delta1_index = 1:Delta1_steps
        % % get carrier density as a function of Ene & B
        % 对E作积分算出carrier density(积分区间为CNP到E之间的区间)，如果E<CNP那么就是负的载流子，如果E>CNP，那么就是正的载流子
        for E_index = 1:E_points
            carrier_density_matrix(E_index, Delta1_index) = sum(DOS_ED_matrix(1:E_index, Delta1_index)) * delta_ene;    
        end
    end
end