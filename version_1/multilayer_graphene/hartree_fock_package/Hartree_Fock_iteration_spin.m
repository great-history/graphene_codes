%% 进行一次Hartree-Fock迭代
function [error, density_matrix_temp] = Hartree_Fock_iteration_spin(H_hf, density_matrix_last, filling_factor)
    % 对角化自旋向上块（check:应该没问题）
    [eigvecs_temp, eigvals_temp] = eig(H_hf); % MATLAB eig默认是从大到小进行排序
    
    % 根据filling_factor选择其中能量最低的filling_factor个态（check:？？？）
    density_matrix_temp = zeros(3);
    
    [~, index_list] = sort(real(diag(eigvals_temp)), 'ascend'); % 加real是因为厄密矩阵的本征值必为实数，但在这里如果不用real的话,计算机会按照模长的大小进行排序
    for ii = 1:filling_factor
        index = index_list(ii);
        
        % 因为同一个本征态，乘上一个相因子后还跟原来的本征态一样，但是density_matrix却大不相同，因此我们需要规定本征态的第一个分量必须限定为正实数
        % 但是在这里又不需要添加相因子，因为 \rho = |><| 把相因子都抵消掉了
        density_matrix_temp = density_matrix_temp + eigvecs_temp(:, index) * eigvecs_temp(:, index)';
        
    end
    
    % 计算新的density matrix（check:？？？）
    % density_matrix_up_temp = density_matrix_up_temp * transpose(eigvecs_up_temp);
    % density_matrix_up_temp = conj(eigvecs_up_temp) * density_matrix_up_temp;
    % density_matrix_down_temp = density_matrix_down_temp * transpose(eigvecs_down_temp);
    % density_matrix_down_temp = conj(eigvecs_down_temp) * density_matrix_down_temp;
    
    % 计算误差（check:应该没问题）
    error_matrix = abs(density_matrix_temp - density_matrix_last);
    b = sum(sum(abs(density_matrix_last).^2 + abs(density_matrix_last).^2));
    a = sum(sum(abs(error_matrix).^2));
    error = a / b;
end