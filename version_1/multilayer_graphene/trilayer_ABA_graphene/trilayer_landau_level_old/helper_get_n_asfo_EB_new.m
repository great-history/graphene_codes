function [carrier_density_matrix] = helper_get_n_asfo_EB_new(Ene_cnp_vecs, DOS_EB_matrix, ene_vecs, delta_ene, E_points, B_steps)
   
    carrier_density_matrix = zeros(E_points, B_steps);  % 载流子浓度作为能量E和磁场B的函数
    
    for B_index = 1:B_steps
        ene_cnp = Ene_cnp_vecs(B_index);

        % 确定CNP在ene_vecs中的位置
        for i = 1:(E_points - 1)
            if ene_cnp >= ene_vecs(i) && ene_cnp <= ene_vecs(i+1)  % 如果ene_vecs是从大到小排列，则要用ene_cnp < ene_vecs(i) && ene_cnp > ene_vecs(i+1)
                cnp_index = i;
                break
            end
        end

        % % get carrier density as a function of Ene & B
        % 确定未定常数(它对应的是从-infinity到E1_bound之间的态密度积分,需要通过CNP来确定)
        const_B = - sum(DOS_EB_matrix(1:cnp_index, B_index)) * delta_ene;
            
        % 对E作积分算出carrier density(积分区间为CNP到E之间的区间)，如果E<CNP那么就是负的载流子，如果E>CNP，那么就是正的载流子
        for E_index = 1:E_points
            carrier_density_matrix(E_index, B_index) = sum(DOS_EB_matrix(1:E_index, B_index)) * delta_ene; 
        end
        
        carrier_density_matrix(:, B_index) = carrier_density_matrix(:, B_index) + const_B;
    end
end