%% 误差函数
function [error, density_matrix_up_temp, density_matrix_down_temp] = Hartree_Fock_iteration(H_hf_up, H_hf_down, density_matrix_up_last, density_matrix_down_last, filling_factor, dims)
    % 对角化自旋向上块（check:应该没问题）
    [eigvecs_up_temp, eigvals_up_temp] = eig(H_hf_up); % MATLAB eig默认是从大到小进行排序
    % 对角化自旋向下块
    [eigvecs_down_temp, eigvals_down_temp] = eig(H_hf_down); % MATLAB eig默认是从大到小进行排序
    
    % 根据filling_factor选择其中能量最低的filling_factor个态（check:？？？）
    num_up = 0;
    num_down = 0;
    density_matrix_up_temp = zeros(dims);
    density_matrix_down_temp = zeros(dims);
    
    eigvals_temp = [real(diag(eigvals_up_temp)); real(diag(eigvals_down_temp))]; % 加real是因为厄密矩阵的本征值必为实数，但在这里如果不用real的话,计算机会按照模长的大小进行排序
    [~, index_list] = sort(eigvals_temp, 'ascend'); 
    for ii = 1:filling_factor
        index = index_list(ii);
        
        % 因为同一个本征态，乘上一个相因子后还跟原来的本征态一样，但是density_matrix却大不相同，因此我们需要规定本征态的第一个分量必须限定为正实数
        % 但是在这里又不需要添加相因子，因为 \rho = |><| 把相因子都抵消掉了
        if index > dims
            num_down = num_down + 1;
            density_matrix_down_temp = density_matrix_down_temp + eigvecs_down_temp(:, index - dims) * eigvecs_down_temp(:, index - dims)';
        else
            num_up = num_up + 1;
            density_matrix_up_temp = density_matrix_up_temp + eigvecs_up_temp(:, index) * eigvecs_up_temp(:, index)';
        end
        
    end
    
    % 计算新的density matrix（不建议使用这个方法，含义不明确）
    % density_matrix_up_temp = density_matrix_up_temp * transpose(eigvecs_up_temp);
    % density_matrix_up_temp = conj(eigvecs_up_temp) * density_matrix_up_temp;
    % density_matrix_down_temp = density_matrix_down_temp * transpose(eigvecs_down_temp);
    % density_matrix_down_temp = conj(eigvecs_down_temp) * density_matrix_down_temp;
    
    % 计算误差（check:应该没问题）
    error_matrix = abs(density_matrix_up_temp - density_matrix_up_last) + abs(density_matrix_down_temp - density_matrix_down_last);
    b = sum(sum(abs(density_matrix_up_last).^2 + abs(density_matrix_down_last).^2));
    a = sum(sum(abs(error_matrix).^2));
    error = a / b;
end