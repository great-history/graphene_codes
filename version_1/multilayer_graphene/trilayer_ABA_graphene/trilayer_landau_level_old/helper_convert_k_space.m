function [kxs, kys] = convert_k_space(N1, N2, k1s, k2s)
    % 转换到直角坐标系:即(k1s, k2s) → (kxs, kys)
    vec_b1 = [1, 1 / sqrt(3)];  % 用直角坐标系表示，省略了系数 2 * pi / a
    vec_b2 = [1, - 1 / sqrt(3)];

    kxs = zeros(N1, N2);
    kys = zeros(N1, N2);
    for k2_ind = 1:N2
        for k1_ind = 1:N1
            k1 = k1s(k1_ind);
            k2 = k2s(k2_ind);

            k_xy = k1 * vec_b1 + k2 * vec_b2;
            kxs(k1_ind, k2_ind) = k_xy(1);
            kys(k1_ind, k2_ind) = k_xy(2);
        end
    end
end