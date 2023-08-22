function input_value_list = get_random_input_list(ub_list, lb_list, num_input)
    % ub_list上边界  lb_list下边界
    input_value_list = zeros(num_input, 1);
    rand_coeff_list = rand(num_input, 1);
    
    for ii = 1:num_input
        input_value_list(ii) = lb_list(ii) + rand_coeff_list(ii) * (ub_list(ii) - lb_list(ii));
    end
end