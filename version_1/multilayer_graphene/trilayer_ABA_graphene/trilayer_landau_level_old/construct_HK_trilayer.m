function [varargout] = construct_HK_trilayer(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, B_field, N_LL)
    if Delta1 == 0
        [HK_b, HK_m, HKp_b, HKp_m] = construct_HK_without_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, B_field, N_LL);
        varargout{1} = HK_b;
        varargout{2} = HK_m;
        varargout{3} = HKp_b;
        varargout{4} = HKp_m;
    else
        [HK,HKp] = construct_HK_with_D(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta1, Delta2, B_field, N_LL);
        varargout{1} = HK;
        varargout{2} = HKp;
    end
end