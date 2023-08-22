function gap = get_interlayer_asymmetry_gap(D_field, B_field, interlayer_d, rho_layer1, rho_layer2)
    % Self-consistently calculated inter-layer asymmetry gap as a function
    % of electric field , only hold at integer fillings
    electron_charge = 0.0; % 单电子的电荷量
    epsilon_0 = 1.0; % 真空介电常数
    h_plank = 1.0;
    
    macro_deg = 2 * electron_charge * B_field / h_plank; % 宏观简并度
    gap = electron_charge * interlayer_d / epsilon_0 * ( D_field + electron_charge * (rho_layer1 - rho_layer2) );
end