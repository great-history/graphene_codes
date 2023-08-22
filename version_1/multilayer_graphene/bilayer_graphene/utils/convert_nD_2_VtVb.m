%% convert nD to two gate
function [Vt, Vb] = convert_nD_2_VtVb(carrier_density, disp_field, Ct, Cb)
    % Ct top gate的电容，单位是10^(12) /cm^2 , 已经除去了电荷量
    % Cb bottom gate的电容，单位是10^(12) /cm^2 , 已经除去了电荷量
    % 基本物理量
    single_electron = 1.602176634 * 10^(-19); % C
    epsilon_0 = 8.85 * 10^(-12); % F/m

    a = 2 * epsilon_0 / single_electron * 10^(-7) * disp_field;
    b = 10^(-12) * carrier_density;
    
    Vt = 1 / (2 * Ct) * (b + a);
    Vb = 1 / (2 * Cb) * (b - a);
end