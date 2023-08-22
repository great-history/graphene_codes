function [density_EB_mat, density_max, density_min] = get_n_asfo_EB(dos_EB_mat, energy_list, Ene_steps, B_steps, delta_ene)
    % delta_ene是能量间隔
    % 由于dos_EB_mat是通过energy cnp shift之后得到的，而在cnp处载流子浓度是0,在这里则是E=0处载流子浓度为0
    density_EB_mat = zeros(B_steps, Ene_steps);
    
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
    
    cnp_index = energy_list == 0;
    % 对态密度作关于E的积分算出载流子浓度(积分区间为CNP到E之间的区间)，如果E<CNP那么就是负的载流子，如果E>CNP，那么就是正的载流子
    density_max = inf; % density_max是不同B下的最大值中的最小值
    density_min = -inf; % density_min是不同B下的最小值中的最大值
    for B_index = 1:B_steps
        density_slice = cumsum(dos_EB_mat(B_index,:)) * delta_ene;
        density_slice = density_slice - density_slice(cnp_index);
        
        density_EB_mat(B_index, :) = density_slice; 
        
        if max(density_slice) < density_max
            density_max = max(density_slice);
        end

        if min(density_slice) > density_min
            density_min = min(density_slice);
        end
    end
    
end

% function [carrier_density_matrix] = helper_get_n_asfo_EB_new(Ene_cnp_vecs, DOS_EB_matrix, ene_vecs, delta_ene, E_points, B_steps)
%    
%     carrier_density_matrix = zeros(E_points, B_steps);  % 载流子浓度作为能量E和磁场B的函数
%     
%     for B_index = 1:B_steps
%         ene_cnp = Ene_cnp_vecs(B_index);
% 
%         % 确定CNP在ene_vecs中的位置
%         for i = 1:(E_points - 1)
%             if ene_cnp >= ene_vecs(i) && ene_cnp <= ene_vecs(i+1)  % 如果ene_vecs是从大到小排列，则要用ene_cnp < ene_vecs(i) && ene_cnp > ene_vecs(i+1)
%                 cnp_index = i;
%                 break
%             end
%         end
% 
%         % % get carrier density as a function of Ene & B
%         % 确定未定常数(它对应的是从-infinity到E1_bound之间的态密度积分,需要通过CNP来确定)
%         const_B = - sum(DOS_EB_matrix(1:cnp_index, B_index)) * delta_ene;
%             
%         % 对E作积分算出carrier density(积分区间为CNP到E之间的区间)，如果E<CNP那么就是负的载流子，如果E>CNP，那么就是正的载流子
%         for E_index = 1:E_points
%             carrier_density_matrix(E_index, B_index) = sum(DOS_EB_matrix(1:E_index, B_index)) * delta_ene; 
%         end
%         
%         carrier_density_matrix(:, B_index) = carrier_density_matrix(:, B_index) + const_B;
%     end
% end