function filling_factor = filling_factor_func_linear(ene, w_bandwidth)
    % only use for linear dos
    if ene <= - w_bandwidth
        filling_factor = - 1;
    elseif (- w_bandwidth < ene) && (ene <= 0)
        filling_factor = - (ene / w_bandwidth)^2;
    elseif (ene < w_bandwidth) && (ene > 0)
        filling_factor = (ene / w_bandwidth)^2;
    else
        filling_factor = 1;
    end
end