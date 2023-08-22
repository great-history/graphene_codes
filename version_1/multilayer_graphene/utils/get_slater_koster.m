function t = get_slater_koster(d_inter, d_intra, theta)
    % a0 is the distance of neighboring A and B sites on graphene
    % d0 is the interlayer spacing of graphite
    % r0 is the decay length of the transfer integral
    a0 = 0.142;
    d0 = 0.335;
    r0 = 0.184 * sqrt(3) * a0;
    
    R = sqrt(d_inter^2 + d_intra^2);
    
    t_pi = 2.7;
    t_sigma = - 0.48;
    
    t = sin(theta) * t_pi * exp(-(R - a0)/r0) + cos(theta) * t_sigma * exp(-(R - d0)/r0);

 end