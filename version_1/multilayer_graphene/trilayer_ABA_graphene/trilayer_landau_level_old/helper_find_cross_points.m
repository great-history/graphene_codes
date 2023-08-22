function [B_cross, cross_index] = helper_find_cross_points(LL_b, LL_m, LLb2_index, LLm0_index, B_fields, B_steps)
    gaps = LL_b(:, LLb2_index) - LL_m(:, LLm0_index);
    cross_index = 0;
    for i = 1:(B_steps - 1)
        if (gaps(i) <= 0 && gaps(i+1) >= 0) || (gaps(i) >= 0 && gaps(i+1) <= 0)  % 变号的点处就是交叉点处
            cross_index = i;
        end
    end
    
    if cross_index == 0 % 此时交叉点不在当前的磁场范围内
        % 因为LLb2是随着磁场线性变化的，而LLm0是不随磁场发生变化的，因此我们还是可以extropolate出B_cross，而cross_index只能记为-1，表示超出范围
        % 利用gaps(1)和gaps(2)
        cross_index = -1;
        B_cross = B_fields(end) - gaps(end) * (B_fields(end) - B_fields(end-1)) / (gaps(end) - gaps(end-1)); 
    else
       % 找到变号位置处后再用一个线性拟合得到交叉点的“严格”位置
        if gaps(cross_index) == 0
            B_cross = B_fields(cross_index);
        elseif gaps(cross_index+1) == 0
            B_cross = B_fields(cross_index+1);
        else
            B_cross = B_fields(cross_index) - gaps(cross_index) * (B_fields(cross_index+1) - B_fields(cross_index)) / (gaps(cross_index+1) - gaps(cross_index));
        end 
    end
end