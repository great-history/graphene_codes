function [carrier_density_matrix] = helper_get_n_asfo_EB(LL_u1, LL_d1, DOS_EB_matrix, ene_vecs, delta_ene, E_points, B_steps)
    % 要得到载流子浓度作为能量E和磁场B的函数，首先我们要找到CNP，要确定CNP，我们要先确定CNP上下最近邻的两条LL，对于一般的石墨烯体系，这应该是确定的两条
    % 比如monolayer就是LLm0_K和LLm0_Kp这两条，而trilayer在Delta1 = 0时则是LLb0_K和LLb1_K这两条
    
    % 算出charge neutrality point处的能量值
    E_cnp_vecs = (LL_u1 + LL_d1) / 2; % LLb0_K与LLb1_K的中间值在低磁场下不一定是CNP,这个其实只在高磁场下成立
    % E_cnp_indexs = zeros(B_steps,1);
    
    carrier_density_matrix = zeros(E_points, B_steps);  % 载流子浓度作为能量E和磁场B的函数
    
    for B_index = 1:B_steps
        ene_cnp = E_cnp_vecs(B_index);
        
        % % get CNP's index as a function of B 
        % 确定CNP在ene_vecs中的位置
        for i = 1:(E_points - 1)
            if ene_cnp >= ene_vecs(i) && ene_cnp <= ene_vecs(i+1)  % 如果ene_vecs是从大到小排列，则要用ene_cnp < ene_vecs(i) && ene_cnp > ene_vecs(i+1)
                ene_cnp_index = i;
                break
            end
        end

        % E_cnp_indexs(B_index) = ene_cnp_index;

        % % get carrier density as a function of Ene & B
        % 对E作积分算出carrier density(积分区间为CNP到E之间的区间)，如果E<CNP那么就是负的载流子，如果E>CNP，那么就是正的载流子
        for E_index = 1:E_points
            if E_index < ene_cnp_index % 负载流子浓度
                carrier_density_matrix(E_index, B_index) = -sum(DOS_EB_matrix(E_index:ene_cnp_index, B_index)) * delta_ene;  % 最好乘delta_ene,保证量纲，同时也保证不会受到meV还是eV单位的影响
            else % 正载流子浓度
                carrier_density_matrix(E_index, B_index) = sum(DOS_EB_matrix(ene_cnp_index:E_index, B_index)) * delta_ene;
            end    
        end
    end
end