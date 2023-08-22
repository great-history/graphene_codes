function decimal_num = num_system_mechine(num_list, weight)
    % 将 num_list 按照一定的位权(weight)转换成一个十进制数
    decimal_num = 0.0;
    max_pow = length(num_list);
    for ii = 1:length(num_list)
        decimal_num = decimal_num + num_list(ii) * (weight)^(max_pow - ii);
    end
end