function [eigvecs_last, vel_dirs, vel_array] = get_bloch_velocity(Vel_k, eigvecs_last, degeneracy, eig_num, first_idxs, last_idxs, dims)
    vel_dirs = zeros(dims, 1);
    vel_array = zeros(dims, 1);
    
    for num = 1:eig_num
        if degeneracy(num) == 1
            first = first_idxs(num);
            vec_now = eigvecs_last(:,first);
            vel = real(dot(vec_now, Vel_k * vec_now));
            % vel = dot(vec_now, Vel_k * vec_now); % 这一行代码用来检查vel是否真的为实数
            vel_array(first) = vel;

            if abs(vel) < 1e-8  % 速度为零，一般对应平带
                vel_dirs(first) = 0;
            elseif vel < 0
                vel_dirs(first) = -1;
            elseif vel > 0
                vel_dirs(first) = 1;
            end

        else

            first = first_idxs(num);
            last = last_idxs(num);
            deg_num = (last - first) + 1;

            vec_degenerate = eigvecs_last(:, first:last);
            vel_matrix = vec_degenerate' * Vel_k * vec_degenerate;  % 速度矩阵
            [eigvecs, eigvals] = eig(vel_matrix);
            
            eigvals = real(diag(eigvals));
            % eigvals = diag(eigvals);  % 这一行代码用来检查vel是否真的为实数
            vel_array(first:last) = eigvals;

            for kk = 1:deg_num
                vel = eigvals(kk);
                if abs(vel) < 1e-8  % 速度为零，一般对应平带
                    vel_dirs(kk + first - 1) = 0;
                elseif vel < 0
                    vel_dirs(kk + first - 1) = -1;
                elseif vel > 0
                    vel_dirs(kk + first - 1) = 1;
                end
            end

            eigvecs_last(:, first:last) = vec_degenerate * eigvecs;

        end
    end
end