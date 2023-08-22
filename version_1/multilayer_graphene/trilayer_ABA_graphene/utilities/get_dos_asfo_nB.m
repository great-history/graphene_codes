function dos_nB_mat = get_dos_asfo_nB(dos_EB_mat, density_EB_mat, density_list, B_steps, density_steps)
    dos_nB_mat = zeros(B_steps, density_steps);
    for B_index = 1:B_steps
        density_slice = density_EB_mat(B_index, :); % 作为x值
        dos_slice = dos_EB_mat(B_index, :); % 作为y值
        % 做插值
        dos_fit_list = interp1(density_slice, dos_slice, density_list);
        dos_nB_mat(B_index, :) = dos_fit_list;
    end
end