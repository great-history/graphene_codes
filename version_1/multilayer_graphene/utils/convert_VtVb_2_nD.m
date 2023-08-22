%% convert two gate to nD
function [carrier_density, disp_field] = convert_VtVb_2_nD(Vt, Vb, Ct, Cb)
    % Ct top gate的电容，单位是10^(12) /cm^2 , 已经除去了电荷量
    % Cb bottom gate的电容，单位是10^(12) /cm^2 , 已经除去了电荷量
    % 基本物理量
    single_electron = 1.602176634 * 10^(-19); % C
    epsilon_0 = 8.85 * 10^(-12); % F/m

    carrier_density = (Ct * Vt + Cb * Vb) * 10^(12);
    disp_field = (Ct * Vt - Cb * Vb) * 10^(12) * 10^(4);
    disp_field = disp_field * single_electron / (2 * epsilon_0);
    disp_field = disp_field * 10^(-9);
end