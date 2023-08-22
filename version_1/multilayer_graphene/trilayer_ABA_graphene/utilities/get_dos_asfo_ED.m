function dos_ED_mat = get_dos_asfo_ED(eig_info_H_select_cell, energy_list, Ene_steps, Delta1_list, Delta1_steps, B_field, ene_width)
    % eigvals_LL_K / eigvals_LL_Kp 是经过energy shifting之后的朗道扇形图
    % energy_list是存放需要计算DOS的各个能量点
    % ene_width是能量展宽（LL broadening）
    % 累加每个朗道能级带来的能量展宽
    
    dos_ED_mat = zeros(Delta1_steps, Ene_steps); % 态密度作为能量E和磁场B的函数，其实DOS与简并度(degeneracy)是一回事
    dos_factor = 2 * B_field / (2 * pi * 25.6^2); % 2来自自旋简并
    
    for D_index = 1:Delta1_steps
        dos_slice = zeros(1, Ene_steps);
        num_vec = eig_info_H_select_cell{D_index, 3};
        eigvals_temp = eig_info_H_select_cell{D_index, 5};
        for ii = 1:num_vec
            delta_ene_list = energy_list - eigvals_temp(ii); % 计算能量差
            delta_ene_list = delta_ene_list.^2;
            delta_ene_list = ene_width^2 + delta_ene_list;
            delta_ene_list = ene_width ./ delta_ene_list;
            dos_slice = dos_slice + dos_factor * delta_ene_list;
        end
        
        dos_ED_mat(D_index, :) = dos_slice;
    end
end