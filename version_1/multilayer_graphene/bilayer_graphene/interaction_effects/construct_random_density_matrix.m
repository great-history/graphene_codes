function density_matrix_random = construct_random_density_matrix(filling_factor, dim, type)
    % 根据filling_factor，利用随机数产生一个满足一定条件的density matrixs
    density_matrix_random = zeros(dim);
    eigvecs_array = zeros(dim, filling_factor); % 存放随机的正交归一化矢量
    
    if type == "complex"  % 复数
        for ii = 1:filling_factor
            eigvec_temp = rand(dim,1) + 1j * rand(dim,1);
            eigvec_temp = normalize(eigvec_temp,'norm',2);

            % 进行schmit正交归一化
            for jj = 1:(ii - 1)
                eigvec_last = eigvecs_array(:,jj);
                inner_dot = dot(eigvec_last, eigvec_temp);
                eigvec_temp = eigvec_temp - inner_dot * eigvec_last;
                eigvec_temp = normalize(eigvec_temp,'norm',2);
            end

            eigvec_temp = normalize(eigvec_temp,'norm',2);
            eigvecs_array(:, ii) = eigvec_temp;
        end

        for ii = 1:filling_factor
            eigvec_temp = eigvecs_array(:, ii);
            density_matrix_random = density_matrix_random + eigvec_temp * eigvec_temp';
        end
    else  % 实数
        for ii = 1:filling_factor
            eigvec_temp = rand(dim,1);
            eigvec_temp = normalize(eigvec_temp,'norm',2);

            % 进行schmit正交归一化
            for jj = 1:(ii - 1)
                eigvec_last = eigvecs_array(:,jj);
                inner_dot = dot(eigvec_last, eigvec_temp);
                eigvec_temp = eigvec_temp - inner_dot * eigvec_last;
                eigvec_temp = normalize(eigvec_temp,'norm',2);
            end

            eigvec_temp = normalize(eigvec_temp,'norm',2);
            eigvecs_array(:, ii) = eigvec_temp;
        end

        for ii = 1:filling_factor
            eigvec_temp = eigvecs_array(:, ii);
            density_matrix_random = density_matrix_random + eigvec_temp * eigvec_temp';
        end
    end
end