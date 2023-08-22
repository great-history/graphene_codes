function mf_grand_pot = mean_field_grand_potential_func(kinetic_flavor_list, chem_pot_flavor_list, Hubbard_U, chem_pot, w_bandwidth, type)
    % get the mean field grand potential per moire unit cell
    % inputs : 
    %   kinetic_flavor_list : kinetic energy per flavor
    %   chem_pot_flavor_list : chemical potential per flavor
    %   Hubbard_U : strength of on-site electron-electron interaction
    %   chem_pot : global chemical potential
    %   type : linear
    
    mf_grand_pot = 0.0;
    num_flavor = 4;
    filling_factor_list = zeros(1,4);
    if type == "linear"
        for jj = 1:num_flavor
            chem_pot_flavor = chem_pot_flavor_list(jj);
            filling_factor_list(jj) = filling_factor_func_linear(chem_pot_flavor, w_bandwidth);            
        end
        
        % matlab最好少用for循环
        % kinetic energy
        mf_grand_pot = mf_grand_pot + sum(kinetic_flavor_list);
        
        % hubbard interaction
        mf_grand_pot = mf_grand_pot + Hubbard_U / 2 * (sum(filling_factor_list))^2;
        mf_grand_pot = mf_grand_pot - Hubbard_U / 2 * dot(filling_factor_list, filling_factor_list); % substract terms of the same flavor
        
        % -\mu N (global chemical potential)
        mf_grand_pot = mf_grand_pot - chem_pot * sum(filling_factor_list);
    end
end