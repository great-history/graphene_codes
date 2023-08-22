function thickness = calc_bn_thickness(capacitance, dielectric_const)
    % 计算电容公式：C = epsilon * S / d，但在这里我们使用的是单位面积诱导出多少载流子(而非电荷量)的电容
    % 即 C / (e * S) = epsilon / (d * e)
    one_electron_coulomb = 1.602176634 * 10^(-19); % C
    epsilon_0 = 8.85 * 10^(-14); % F/cm
    epsilon_now = epsilon_0 * dielectric_const;
    
    thickness = epsilon_now / (capacitance * one_electron_coulomb);
end

% example
% calc_bn_thickness(0.4480 * 10^(12), 3)  bn的相对介电常数设为3