function [density_ED_mat, density_max, density_min] = get_n_asfo_ED(dos_ED_mat, energy_list, Ene_steps, Delta1_steps, delta_ene)
    % delta_ene是能量间隔
    % 由于dos_EB_mat是通过energy cnp shift之后得到的，而在cnp处载流子浓度是0,在这里则是E=0处载流子浓度为0
    density_ED_mat = zeros(Delta1_steps, Ene_steps);
    
%     if energy_list(1) > energy_list(end)
%         sign_plus = - 1; % 往能量大的方向是指标减1
%     else
%         sign_plus = + 1; % 往能量大的方向是指标加1
%     end
%     
%     [~, ind1] = min(abs(energy_list));
%     if energy_list(ind1) < 0
%         ind2 = ind1 + sign_plus;
%     else
%         ind2 = ind1 - sign_plus;
%     end
    
    cnp_index = find(energy_list == 0);
    % 对态密度作关于E的积分算出载流子浓度(积分区间为CNP到E之间的区间)，如果E<CNP那么就是负的载流子，如果E>CNP，那么就是正的载流子
    density_max = inf; % density_max是不同B下的最大值中的最小值
    density_min = -inf; % density_min是不同B下的最小值中的最大值
    for D_index = 1:Delta1_steps
        density_slice = cumsum(dos_ED_mat(D_index, :)) * delta_ene;
        density_slice = density_slice - density_slice(cnp_index);
        
        density_ED_mat(D_index, :) = density_slice; 
        
        if max(density_slice) < density_max
            density_max = max(density_slice);
        end

        if min(density_slice) > density_min
            density_min = min(density_slice);
        end
    end
    
end