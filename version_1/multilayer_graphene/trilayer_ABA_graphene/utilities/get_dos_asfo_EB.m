function dos_EB_mat = get_dos_asfo_EB(eig_info_H_select_cell, energy_list, Ene_steps, B_fields_list, B_steps, ene_width)
    % eigvals_LL_K / eigvals_LL_Kp 是经过energy shifting之后的朗道扇形图
    % energy_list是存放需要计算DOS的各个能量点
    % ene_width是能量展宽（LL broadening）
    % 累加每个朗道能级带来的能量展宽
    
    dos_EB_mat = zeros(B_steps, Ene_steps); % 态密度作为能量E和磁场B的函数，其实DOS与简并度(degeneracy)是一回事
    for B_index = 1:B_steps
        B_field = B_fields_list(B_index);
        dos_factor = 2 * B_field / (2 * pi * 25.6^2); % 2来自自旋简并
        
        dos_slice = zeros(1, Ene_steps);
        num_vec = eig_info_H_select_cell{B_index, 3};
        eigvals_temp = eig_info_H_select_cell{B_index, 5};
        for ii = 1:num_vec
            delta_ene_list = energy_list - eigvals_temp(ii); % 计算能量差
            delta_ene_list = delta_ene_list.^2;
            delta_ene_list = ene_width^2 + delta_ene_list;
            delta_ene_list = ene_width ./ delta_ene_list;
            dos_slice = dos_slice + dos_factor * delta_ene_list;
        end
        
        dos_EB_mat(B_index, :) = dos_slice;
    end
end

% function [ene_vecs, DOS_EB_matrix] = get_dos_asfo_EB_new(LL_ene_vecs, B_fields, E1_bounds, E2_bounds, E_points, B_steps, dims, Gamma)
%     % E1_bounds是需要计算的那些能量点, 但在E2_bounds范围内的LL的影响都应该被包括进去
%     % LL_ene_vecs包含m_K/m_Kp/b_K/b_Kp在内的所有LL
%     
%     ene_vecs = linspace(E1_bounds(1), E1_bounds(2), E_points);
%     DOS_EB_matrix = zeros(E_points, B_steps); % 态密度作为能量E和磁场B的函数，其实DOS与简并度(degeneracy)是一回事
%     
%     % get DOS as a function of Ene & B
%     for B_index = 1:B_steps
%         B_field = B_fields(B_index);
%         
%         for E_index = 1:E_points
%             % 得到DOS
%             ene = ene_vecs(E_index);
%             dos = 0.0;
%             
%             for ll = 1:dims
%                 if LL_ene_vecs(B_index, ll) < E2_bounds(1) && LL_ene_vecs(B_index, ll) > E2_bounds(2)
%                     dos = dos + sum((Gamma / 2) / ((Gamma / 2)^2 + (ene - LL_ene_vecs(B_index, ll))^2));
%                 end
%             end
%             
%             dos = dos * B_field / (pi * 25.66)^2;
%             DOS_EB_matrix(E_index, B_index) = dos;
%         end
%     end
% end