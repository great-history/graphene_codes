function plot_cnp_line(current_fig, eigvals_LL_xxx, xxx_list, cnp_LL_index, marker_size)
    % 做CNP 中心线（由于K valley和Kp valley各有dims，所以总共有2*dims(偶数)条LLs，因此cnp应为其中两条相邻的LL之间的中心）
    % 我们如果确定了某个Delta1下的cnp对应的两个LL的index，那么所有Delta1下的cnp都对应这两个index下的LL，虽然有可能不是同一条LL
    % cnp_LL_index = dims;
    xxx_steps = length(xxx_list);
    cnp_energy_list = zeros(xxx_steps, 1);
    
    % 首先把K和Kp的LLs合并
    % eigvals_LL_Delta1 = [eigvals_LL_K_Delta1, eigvals_LL_Kp_Delta1];
    % 然后通过cnp_LL_index计算出cnp_LL的位置
    for xxx_index = 1:xxx_steps
        eigvals_LL_list = eigvals_LL_xxx(xxx_index, :);
        eigvals_LL_list = sort(eigvals_LL_list, 'ascend');
        cnp_energy_list(xxx_index) = 1 / 2 * (eigvals_LL_list(cnp_LL_index) + eigvals_LL_list(cnp_LL_index + 1));
    end
    
    figure(current_fig) % 切换当前的figure
    hold on
    % plot(Delta1_list, cnp_energy_list)
    plot(xxx_list, cnp_energy_list, 'r*', "MarkerSize", marker_size)
end