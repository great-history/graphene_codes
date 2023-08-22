function [HbKp_vx, HbKp_vy] = construct_HbKp_velocity_op_with_rotate(gamma0, gamma1, gamma3, gamma4, delta_dimer, akx, aky, angle, flag_chiral)
    % vf = 2.1354; % eV为单位
    % flag_chiral : 描述了手性:是AB还是BA
    v3 = gamma3 * sqrt(3) / 2;
    v4 = gamma4 * sqrt(3) / 2;
    
    phase = exp(1j*angle); % angle已经转化为弧度制了
    
    %% 单层
    [HmKp_layer1_vx, HmKp_layer1_vy] = construct_HmKp_velocity_op_with_rotate(gamma0, akx, aky, angle);
    [HmKp_layer2_vx, HmKp_layer2_vy] = construct_HmKp_velocity_op_with_rotate(gamma0, akx, aky, angle);
    
    H_inter_vx = zeros(2,2);
    H_inter_vy = zeros(2,2);
    
    H_inter_vx(1,1) = v4 * conj(phase);
    H_inter_vy(1,1) = v4 * 1j * conj(phase);
    
    H_inter_vx(2,1) = v3 * phase;
    H_inter_vy(2,1) = - v3 * 1j * phase;
    
    H_inter_vx(2,2) = v4 * conj(phase);
    H_inter_vy(2,2) = v4 * 1j * conj(phase);
    
    %% 得到哈密顿量
    HbKp_vx = zeros(4,4);
    HbKp_vy = zeros(4,4);
    
    if flag_chiral % AB-stacking
        HbKp_vx(1:2,1:2) = HmKp_layer1_vx;
        HbKp_vx(3:4,3:4) = HmKp_layer2_vx;
        HbKp_vx(3:4,1:2) = H_inter_vx;
        HbKp_vx(1:2,3:4) = H_inter_vx';
        
        HbKp_vy(1:2,1:2) = HmKp_layer1_vy;
        HbKp_vy(3:4,3:4) = HmKp_layer2_vy;
        HbKp_vy(3:4,1:2) = H_inter_vy;
        HbKp_vy(1:2,3:4) = H_inter_vy';
        
    else % BA-stacking
        HbKp_vx(1:2,1:2) = HmKp_layer2_vx;
        HbKp_vx(3:4,3:4) = HmKp_layer1_vx;
        HbKp_vx(3:4,1:2) = H_inter_vx';
        HbKp_vx(1:2,3:4) = H_inter_vx;
        
        HbKp_vy(1:2,1:2) = HmKp_layer2_vy;
        HbKp_vy(3:4,3:4) = HmKp_layer1_vy;
        HbKp_vy(3:4,1:2) = H_inter_vy';
        HbKp_vy(1:2,3:4) = H_inter_vy;
    end
end