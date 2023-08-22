function [HmK_vx, HmK_vy] = construct_HmK_velocity_op_with_rotate(gamma0, akx, aky, angle)
    % 得到速度算符，H_m_K关于(akx, aky)的偏导,不是严格的速度算符，少乘一个晶格常数a
    % 对于石墨烯的低能模型，HmK_vx, HmK_vy 是与 (akx, aky)无关的矩阵
    phase = exp(-1j*angle); % angle已经转化为弧度制了
    vf = gamma0 * sqrt(3) / 2;
    
    HmK_vx = zeros(2,2);
    HmK_vy = zeros(2,2);
    
    HmK_vx(1,2) = vf * phase;
    HmK_vx(2,1) = conj(HmK_vx(1,2));
    
    HmK_vy(1,2) = vf * 1j * phase;
    HmK_vy(2,1) = conj(HmK_vy(1,2));
    
end