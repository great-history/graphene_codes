%% 得到速度算符(仅单层石墨烯的即可)
[HmK_layer1_vx, HmK_layer1_vy] = construct_HmK_velocity_op_with_rotate(gamma0, 0, 0, -theta / 2); % 第一层顺时针转theta/2
[HmK_layer2_vx, HmK_layer2_vy] = construct_HmK_velocity_op_with_rotate(gamma0, 0, 0, theta / 2); % 第二层逆时针转theta/2
[HmKp_layer2_vx, HmKp_layer2_vy] = construct_HmKp_velocity_op_with_rotate(gamma0, 0, 0, theta / 2); % 第二层逆时针转theta/2
[HmKp_layer1_vx, HmKp_layer1_vy] = construct_HmKp_velocity_op_with_rotate(gamma0, 0, 0, -theta / 2); % 第一层顺时针转theta/2

%% 转角石墨烯的群速度期望
tbg_K_vx = zeros(1, num_K_low_band);
tbg_K_vy = zeros(1, num_K_low_band);
tbg_Kp_vx = zeros(1, num_Kp_low_band);
tbg_Kp_vy = zeros(1, num_Kp_low_band);

for i = 1:num_K_low_band
    tbg_K_vx = 0.0;
    tbg_K_vy = 0.0;
    tbg_Kp_vx = 0.0;
    tbg_Kp_vy = 0.0;
    
    vec_K_temp = eig_vecs_K_low_ene_last(1, :, i);  % 它的尺寸是(1, dims)，是一个行向量
    vec_K_temp = transpose(vec_K_temp);  % 将它转换为一个列向量，转置一下
    vec_Kp_temp = eig_vecs_Kp_low_ene_last(1, :, i);  % 它的尺寸是(1, dims)，是一个行向量
    vec_Kp_temp = transpose(vec_Kp_temp);  % 将它转换为一个列向量，转置一下
    
    for j = 1:num_q_couple
        % layer1
        vec_temp = vec_K_temp(2*j-1:2*j);
        tbg_K_vx = tbg_K_vx + vec_temp' * HmK_layer1_vx * vec_temp;
        tbg_K_vy = tbg_K_vy + vec_temp' * HmK_layer1_vy * vec_temp;
        
        vec_temp = vec_Kp_temp(2*j-1:2*j);
        tbg_Kp_vx = tbg_Kp_vx + vec_temp' * HmKp_layer1_vx * vec_temp;
        tbg_Kp_vy = tbg_Kp_vy + vec_temp' * HmKp_layer1_vy * vec_temp;
    end
    
    for j = (num_q_couple + 1):(2 * num_q_couple)
        % layer2
        vec_temp = vec_K_temp(2*j-1:2*j);
        tbg_K_vx = tbg_K_vx + vec_temp' * HmK_layer2_vx * vec_temp;
        tbg_K_vy = tbg_K_vy + vec_temp' * HmK_layer2_vy * vec_temp;
        
        vec_temp = vec_Kp_temp(2*j-1:2*j);
        tbg_Kp_vx = tbg_Kp_vx + vec_temp' * HmKp_layer2_vx * vec_temp;
        tbg_Kp_vy = tbg_Kp_vy + vec_temp' * HmKp_layer2_vy * vec_temp;
    end
    
    tbg_K_vx_array(1,i) = real(tbg_K_vx);
    tbg_K_vy_array(1,i) = real(tbg_K_vy);
    
    tbg_Kp_vx_array(1,i) = real(tbg_Kp_vx);
    tbg_Kp_vy_array(1,i) = real(tbg_Kp_vy);
end

%% 发生简并时计算非对角元
num_deg = 2;  % 简并度
tbg_K_vx_deg = zeros(num_deg);
tbg_K_vy_deg = zeros(num_deg);
for ii = 1:num_deg
    vec_K_temp1 = eig_vecs_K_low_ene_last(1, :, ii);  % 它的尺寸是(1, dims)，是一个行向量
    vec_K_temp1 = transpose(vec_K_temp1);  % 将它转换为一个列向量，转置一下
    
    for jj = 1:num_deg
        tbg_K_vx = 0.0;
        tbg_K_vy = 0.0;
    
        vec_K_temp2 = eig_vecs_K_low_ene_last(1, :, jj);  % 它的尺寸是(1, dims)，是一个行向量
        vec_K_temp2 = transpose(vec_K_temp2);  % 将它转换为一个列向量，转置一下
        
        % 计算(ii,jj)矩阵元
        for j = 1:num_q_couple
            % layer1
            vec_temp1 = vec_K_temp1(2*j-1:2*j);
            vec_temp2 = vec_K_temp2(2*j-1:2*j);
            tbg_K_vx = tbg_K_vx + vec_temp1' * HmK_layer1_vx * vec_temp2;
            tbg_K_vy = tbg_K_vy + vec_temp1' * HmK_layer1_vy * vec_temp2;
        end

        for j = (num_q_couple + 1):(2 * num_q_couple)
            % layer2
            vec_temp1 = vec_K_temp1(2*j-1:2*j);
            vec_temp2 = vec_K_temp2(2*j-1:2*j);
            tbg_K_vx = tbg_K_vx + vec_temp1' * HmK_layer2_vx * vec_temp2;
            tbg_K_vy = tbg_K_vy + vec_temp1' * HmK_layer2_vy * vec_temp2;
        end
        
        tbg_K_vx_deg(ii, jj) = tbg_K_vx;
        tbg_K_vy_deg(ii, jj) = tbg_K_vy;
    end
    
end

% 进行对角化
eig_vecs_K_low_ene_deg = zeros(dims, num_deg);
[eigvecs_vx, eigvals_vx] = eig(tbg_K_vx_deg);
[eigvecs_vy, eigvals_vy] = eig(tbg_K_vy_deg);


%% 预测下一个k点的能量值
delta_akx = k_vertice_list(1,2) - k_vertice_list(1,1);
delta_aky = k_vertice_list(2,2) - k_vertice_list(2,1);

%%
% eps = 0.99;
% 
% order_now = [];
% for jj = 1:num_K_low_band
%     vec_last = eig_vecs_K_low_ene_last(1, :, jj);
%     val_last = eig_vals_K_low_ene_last(1,jj);
%     
%     degeneracy = 0;
%     for ii = 1:dims
%         vec_now = eigvecs_K(:, ii);
%         val_now = eig_vals_K_temp(ii);
%         
%         if abs(val_now - val_last) < 0.01
%             inner_dot = abs(dot(vec_now, vec_last));
% 
%             if inner_dot >= eps  % 不存在简并            
%                 % 更新low_ene_last
%                 
%                 degeneracy = degeneracy + 1;
%                 order_now = [order_now, ii];
%                 disp(degeneracy)
%                 break
%             elseif inner_dot >= 0.1 && inner_dot < eps  % 存在简并
%                 degeneracy = degeneracy + 1;
%                 if ismember(ii, order_now)
%                     continue
%                 else
%                     order_now = [order_now, ii];
%                 end
%                 disp(degeneracy)
%             end
%         end
%     end
% end

% % 更新low_ene_last
% for jj = 1:num_low_band
%     % 有时候order_now中的数目可能比num_low_band多，但我们只要保证能量最低的那些态是连续的就行
%     eig_vecs_K_low_ene_last(1,:,jj) = eigvecs_K(:, order_now(jj));
%     eig_vals_K_low_ene_last(1,jj) = eig_vals_K_low_ene_temp(order_now(jj));
% end
%     
% % 存入low_ene
% eig_vals_K_low_ene(i, :) = eig_vals_K_low_ene_last;
% eig_vecs_K_low_ene(i, :, :) = eig_vecs_K_low_ene_last;