function [tbg_K_vx_array, tbg_K_vy_array] = get_tbg_K_velocity(HmK_layer1_vx, HmK_layer1_vy, HmK_layer2_vx, HmK_layer2_vy, eig_vecs_K_low_ene_last)
    %% 单层石墨烯的得到速度算符:HmK_layer1_vx, HmK_layer1_vy, HmK_layer2_vx, HmK_layer2_vy
    
    %% 转角石墨烯的群速度期望
    num_K_low_band = size(eig_vecs_K_low_ene_last, 3);
    tbg_K_vx_array = zeros(1, num_K_low_band);
    tbg_K_vy_array = zeros(1, num_K_low_band);
    
    for i = 1:num_K_low_band
        tbg_K_vx = 0.0;
        tbg_K_vy = 0.0;

        vec_K_temp = eig_vecs_K_low_ene_last(1, :, i);  % 它的尺寸是(1, dims)，是一个行向量
        vec_K_temp = transpose(vec_K_temp);  % 将它转换为一个列向量，转置一下

        for j = 1:num_q_couple
            % layer1
            vec_temp = vec_K_temp(2*j-1:2*j);
            tbg_K_vx = tbg_K_vx + vec_temp' * HmK_layer1_vx * vec_temp;
            tbg_K_vy = tbg_K_vy + vec_temp' * HmK_layer1_vy * vec_temp;
        end

        for j = (num_q_couple + 1):(2 * num_q_couple)
            % layer2
            vec_temp = vec_K_temp(2*j-1:2*j);
            tbg_K_vx = tbg_K_vx + vec_temp' * HmK_layer2_vx * vec_temp;
            tbg_K_vy = tbg_K_vy + vec_temp' * HmK_layer2_vy * vec_temp;
        end

        tbg_K_vx_array(1,i) = real(tbg_K_vx);
        tbg_K_vy_array(1,i) = real(tbg_K_vy);
    end
end