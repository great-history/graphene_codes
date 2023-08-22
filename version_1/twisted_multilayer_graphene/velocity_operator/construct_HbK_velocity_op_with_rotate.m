function [HbK_vx, HbK_vy] = construct_HbK_velocity_op_with_rotate(gamma0, gamma1, gamma3, gamma4, delta_dimer, akx, aky, angle, flag_chiral)
    % vf = 2.1354; % eV为单位
    % flag_chiral : 描述了手性:是AB还是BA
    v3 = gamma3 * sqrt(3) / 2;
    v4 = gamma4 * sqrt(3) / 2;
    
    phase = exp(-1j*angle); % angle已经转化为弧度制了
    
    %% 单层
    % 由于速度算符要对k求导，而delta_dimer / gamma1 不是k的函数，所以delta_dimer / gamma1 不会出现在速度算符中
    [HmK_layer1_vx, HmK_layer1_vy] = construct_HmK_velocity_op_with_rotate(gamma0, akx, aky, angle);
    [HmK_layer2_vx, HmK_layer2_vy] = construct_HmK_velocity_op_with_rotate(gamma0, akx, aky, angle);
    
    H_inter_vx = zeros(2,2);
    H_inter_vy = zeros(2,2);
    
    H_inter_vx(1,1) = - v4 * conj(phase);
    H_inter_vy(1,1) = v4 * 1j * conj(phase);
    
    H_inter_vx(2,1) = - v3 * phase;
    H_inter_vy(2,1) = - v3 * 1j * phase;
    
    H_inter_vx(2,2) = - v4 * conj(phase);
    H_inter_vy(2,2) = v4 * 1j * conj(phase);
    
    %% 得到哈密顿量
    HbK_vx = zeros(4,4);
    HbK_vy = zeros(4,4);
    if flag_chiral % AB-stacking
        HbK_vx(1:2,1:2) = HmK_layer1_vx;
        HbK_vx(3:4,3:4) = HmK_layer2_vx;
        HbK_vx(3:4,1:2) = H_inter_vx;
        HbK_vx(1:2,3:4) = H_inter_vx';
        
        HbK_vy(1:2,1:2) = HmK_layer1_vy;
        HbK_vy(3:4,3:4) = HmK_layer2_vy;
        HbK_vy(3:4,1:2) = H_inter_vy;
        HbK_vy(1:2,3:4) = H_inter_vy';
    else % BA-stacking
        HbK_vx(1:2,1:2) = HmK_layer2_vx;
        HbK_vx(3:4,3:4) = HmK_layer1_vx;
        HbK_vx(3:4,1:2) = H_inter_vx';
        HbK_vx(1:2,3:4) = H_inter_vx;
        
        HbK_vy(1:2,1:2) = HmK_layer2_vy;
        HbK_vy(3:4,3:4) = HmK_layer1_vy;
        HbK_vy(3:4,1:2) = H_inter_vy';
        HbK_vy(1:2,3:4) = H_inter_vy;
    end
end