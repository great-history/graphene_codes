function helper_find_low_ene_bands(flag_start, num_low_band)
   % 对低能能带进行排列
    if flag_start % 第一个k点
        [~, new_order] = sort(abs(eig_vals_Kp_low_ene_temp));
        count = 0;
        for jj = 1:dims
            index = new_order(jj);
            if abs(eig_vals_Kp_low_ene_temp(index)) < low_ene_bound
                count = count + 1;
            end
        end

        % 数低能能带的数目
        num_low_band = count;
        if num_low_band == 0
            % 如果低能能带的数目等于0则无法进行排序
            return
        end

        % 定义存储low_ene的几个变量
        eig_vals_Kp_low_ene = zeros(num_k, num_low_band);
        eig_vecs_Kp_low_ene = zeros(num_k, dims, num_low_band);  % 知道本征态可以用来计算Berry curvature / orbital magnetism / Chern number
        eig_vecs_Kp_low_ene_last = zeros(1, dims, num_low_band);
        eig_vals_Kp_low_ene_last = zeros(1, num_low_band);

        % 更新low_ene_last
        eig_vecs_Kp_low_ene_last(1,:,:) = eigvecs_Kp(:, new_order(1:num_low_band));
        eig_vals_Kp_low_ene_last(1,:) = eig_vals_Kp_low_ene_temp(new_order(1:num_low_band));

        % 存入low_ene
        eig_vals_Kp_low_ene(i, :) = eig_vals_Kp_low_ene_last;
        eig_vecs_Kp_low_ene(i, :, :) = eig_vecs_Kp_low_ene_last;
    else
        if num_low_band == 0
            continue
        end

        % 寻找low_ene
        eps = sqrt(0.5); %

        order_now = [];
        for jj = 1:num_low_band
            vec_last = eig_vecs_Kp_low_ene_last(1, :, jj);
            val_last = eig_vals_Kp_low_ene_last(1,jj);

            degeneracy = 0;
            for ii = 1:dims
                vec_now = eigvecs_Kp(:, ii);
                val_now = eig_vals_Kp_low_ene_temp(ii);

                if abs(val_now - val_last) < 0.01
                    inner_dot = abs(dot(vec_now, vec_last));

                    if inner_dot >= eps
                        if ismember(ii, order_now)
                            continue
                        else
                            order_now = [order_now, ii];
                        end
                    end
                end
            end
            % disp(degeneracy)
        end

        % 更新low_ene_last
        for jj = 1:num_low_band  
            % 有时候order_now中的数目可能比num_low_band多，但我们只要保证能量最低的那些态是连续的就行
            eig_vecs_Kp_low_ene_last(1,:,jj) = eigvecs_Kp(:, order_now(jj));
            eig_vals_Kp_low_ene_last(1,jj) = eig_vals_Kp_low_ene_temp(order_now(jj));
        end

        % 存入low_ene
        eig_vals_Kp_low_ene(i, :) = eig_vals_Kp_low_ene_last;
        eig_vecs_Kp_low_ene(i, :, :) = eig_vecs_Kp_low_ene_last;
    end 
end