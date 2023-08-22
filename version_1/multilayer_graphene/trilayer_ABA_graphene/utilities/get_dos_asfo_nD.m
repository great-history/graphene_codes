function dos_nD_mat = get_dos_asfo_nD(dos_ED_mat, density_ED_mat, density_list, Delta1_steps, density_steps)
    dos_nD_mat = zeros(Delta1_steps, density_steps);
    for D_index = 1:Delta1_steps
        density_slice = density_ED_mat(D_index, :); % 作为x值
        dos_slice = dos_ED_mat(D_index, :); % 作为y值
        % 做插值
        dos_fit_list = interp1(density_slice, dos_slice, density_list);
        dos_nD_mat(D_index, :) = dos_fit_list;
    end
end
