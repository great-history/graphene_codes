function ind = find_closest_eigval(val_list, target_val)
    % 找到val_list中与target最接近的那个值
    diff_list = abs(val_list - target_val);
    [~, ind] = min(diff_list);
end