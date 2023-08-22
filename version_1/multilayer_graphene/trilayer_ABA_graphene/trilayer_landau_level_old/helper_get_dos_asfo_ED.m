function [ene_vecs, DOS_ED_matrix] = helper_get_dos_asfo_ED(LL_K, LL_Kp, E_upper_bound, E_lower_bound, E_points, Delta1_steps, B_field, Gamma)
    % 一次只能处理一个K和一个Kp，因此对于trilayer需要调用两次这个函数
    ene_vecs = linspace(E_lower_bound, E_upper_bound, E_points);
    
    DOS_ED_matrix = zeros(E_points, Delta1_steps); % 态密度作为能量E和磁场B的函数，其实DOS与简并度(degeneracy)是一回事
    % get DOS as a function of Ene & B
    for Delta1_index = 1:Delta1_steps
        LLs_K = LL_K(Delta1_index, :);
        LLs_Kp = LL_Kp(Delta1_index, :);

        for E_index = 1:E_points
            % 得到DOS
            ene = ene_vecs(E_index);
            dos = 0.0;
            
            dos_ll_K = sum((Gamma / 2) ./ ((Gamma / 2)^2 + (ene - LLs_K).^2));
            dos_ll_Kp = sum((Gamma / 2) ./ ((Gamma / 2)^2 + (ene - LLs_Kp).^2));
            dos = dos + dos_ll_K + dos_ll_Kp;
            
            dos = dos * B_field / (pi * 25.66)^2;
            DOS_ED_matrix(E_index, Delta1_index) = dos;
        end
    end
end