% 注：其实在计算时没必要把H_tbg_K_vx, H_tbg_K_vy, H_tbg_Kp_vx,
% H_tbg_Kp_vy这些速度算符具体表示出来，用单层石墨烯的速度算符其实就足够了。
function [H_tbg_K_vx, H_tbg_K_vy, H_tbg_Kp_vx, H_tbg_Kp_vy] = construct_tbg_velocity_op(onsite_a, onsite_b, gamma0, w_aa, w_ab, akx, aky, theta, num_q_couple, bm_sup_vecs, bm_sup_mn)
   % dims = 2 * 2 * num_q_couple : 2是layer, 2是sublattice, N是耦合在一起的波矢数目
   dims = 2 * 2 * num_q_couple;
   
   % transfer momentum
   q_norm = 8 * pi * abs(sin(theta / 2)) / 3; % 乘以晶格常数a，保证波矢无量纲
   q_b = q_norm * [0; -1];
   
   %% construct the 1st layer ham @ valley K
   H_K_layer1_vx = zeros(dims / 2);
   H_K_layer1_vy = zeros(dims / 2);
   [HmK_vx, HmK_vy] = construct_HmK_velocity_op_with_rotate(gamma0, 0, 0, -theta / 2); % 第一层顺时针转theta/2
   for i = 1:num_q_couple
       % 第i个参与耦合的波矢对应的速度
       H_K_layer1_vx((2 * i - 1):2*i, (2 * i - 1):2*i) = HmK_vx;
       H_K_layer1_vy((2 * i - 1):2*i, (2 * i - 1):2*i) = HmK_vy;
   end
   
   %% construct the 2nd layer ham @ valley K
   H_K_layer2_vx = zeros(dims / 2);
   H_K_layer2_vy = zeros(dims / 2);
   [HmK_vx, HmK_vy] = construct_HmK_velocity_op_with_rotate(gamma0, 0, 0, theta / 2); % 第一层顺时针转theta/2
   for i = 1:num_q_couple
       % 第i个参与耦合的波矢对应的速度
       H_K_layer2_vx((2 * i - 1):2*i, (2 * i - 1):2*i) = HmK_vx;
       H_K_layer2_vy((2 * i - 1):2*i, (2 * i - 1):2*i) = HmK_vy;
   end
   
   %% construct the 2nd layer ham @ valley Kp
   H_Kp_layer2_vx = zeros(dims / 2);
   H_Kp_layer2_vy = zeros(dims / 2);
   [HmKp_vx, HmKp_vy] = construct_HmKp_velocity_op_with_rotate(gamma0, 0, 0, theta / 2); % 第二层逆时针转theta/2
   for i = 1:num_q_couple
       % 第i个参与耦合的波矢对应的速度
       H_Kp_layer2_vx((2 * i - 1):2*i, (2 * i - 1):2*i) = HmKp_vx;
       H_Kp_layer2_vy((2 * i - 1):2*i, (2 * i - 1):2*i) = HmKp_vy;
   end
   
   %% construct the 1st layer ham @ valley Kp
   H_Kp_layer1_vx = zeros(dims / 2);
   H_Kp_layer1_vy = zeros(dims / 2);
   [HmKp_vx, HmKp_vy] = construct_HmKp_velocity_op_with_rotate(gamma0, 0, 0, -theta / 2); % 第一层顺时针转theta/2
   for i = 1:num_q_couple
       % 第i个参与耦合的波矢对应的速度
       H_Kp_layer1_vx((2 * i - 1):2*i, (2 * i - 1):2*i) = HmKp_vx;
       H_Kp_layer1_vy((2 * i - 1):2*i, (2 * i - 1):2*i) = HmKp_vy;
   end
   
   %% get the whole ham
   H_tbg_K_vx = [H_K_layer1_vx, zeros(dims / 2); zeros(dims / 2), H_K_layer2_vx];
   H_tbg_K_vy = [H_K_layer1_vy, zeros(dims / 2); zeros(dims / 2), H_K_layer2_vy];
   H_tbg_Kp_vx = [H_Kp_layer1_vx, zeros(dims / 2); zeros(dims / 2), H_Kp_layer2_vx];
   H_tbg_Kp_vy = [H_Kp_layer1_vy, zeros(dims / 2); zeros(dims / 2), H_Kp_layer2_vy];
end