function [ene_vecs, DOS_EB_matrix] = helper_get_dos_asfo_EB_new(LL_ene_vecs, B_fields, E1_bounds, E2_bounds, E_points, B_steps, dims, Gamma)
    % E1_bounds是需要计算的那些能量点, 但在E2_bounds范围内的LL的影响都应该被包括进去
    % LL_ene_vecs包含m_K/m_Kp/b_K/b_Kp在内的所有LL
    
    ene_vecs = linspace(E1_bounds(1), E1_bounds(2), E_points);
    DOS_EB_matrix = zeros(E_points, B_steps); % 态密度作为能量E和磁场B的函数，其实DOS与简并度(degeneracy)是一回事
    
    % get DOS as a function of Ene & B
    for B_index = 1:B_steps
        B_field = B_fields(B_index);
        
        for E_index = 1:E_points
            % 得到DOS
            ene = ene_vecs(E_index);
            dos = 0.0;
            
            for ll = 1:dims
                if LL_ene_vecs(B_index, ll) < E2_bounds(1) && LL_ene_vecs(B_index, ll) > E2_bounds(2)
                    dos = dos + sum((Gamma / 2) / ((Gamma / 2)^2 + (ene - LL_ene_vecs(B_index, ll))^2));
                end
            end
            
            dos = dos * B_field / (pi * 25.66)^2;
            DOS_EB_matrix(E_index, B_index) = dos;
        end
    end
end