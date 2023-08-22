% function H_m_Kp = construct_HmKp_with_rotate(vf, onsite_a, onsite_b, akx, aky, angle)
%     % vf = 2.1354; % eV为单位
%     phase = exp(1j*angle); % angle已经转化为弧度制了
%     % vf = gamma0 * sqrt(3) / 2;
%     
%     H_m_Kp = zeros(2,2);
%     H_m_Kp(1,1) = onsite_a;
%     H_m_Kp(2,2) = onsite_b;
%     H_m_Kp(1,2) = - vf * (akx - 1j * aky) * phase;
%     H_m_Kp(2,1) = conj(H_m_Kp(1,2));
% end


function H_m_Kp = construct_HmKp_with_rotate(gamma0, akx, aky, angle)
    % vf = 2.1354; % eV为单位
    phase = exp(1j*angle); % angle已经转化为弧度制了
    vf = gamma0 * sqrt(3) / 2;
    
    H_m_Kp = zeros(2,2);
    H_m_Kp(1,2) = - vf * (akx - 1j * aky) * phase;
    H_m_Kp(2,1) = conj(H_m_Kp(1,2));
end