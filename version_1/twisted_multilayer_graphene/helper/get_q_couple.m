function [b_moire_sup_mn, b_moire_sup_vecs] = get_q_couple(theta, q_trunc)
    % 任何一个倒格矢都可以表示为b_moire_l和b_moire_r的整数线性组合
    % b_moire_l和b_moire_r之间的夹角为60°
    b_moire_norm = 8 * pi * abs(sin(theta / 2)) / sqrt(3); % 无量纲
    b_moire_r = [ 1/2; sqrt(3)/2];  % b moire reciprocal vector at the right side
    b_moire_l = [-1/2; sqrt(3)/2];  % b moire reciprocal vector at the left side
    
    % 首先计算1/6的区域，之后再用对称性得剩余的5/6
    b_moire_sup_mn = []; % sup refers to superposition
    b_moire_sup_vecs = []; % sup refers to superposition
    
    for i = 0:q_trunc
        for j = 0:q_trunc
            b_moire_sup = i * b_moire_r + j * b_moire_l;
            b_moire_sup_norm = norm(b_moire_sup);

            if b_moire_sup_norm <= q_trunc
                b_moire_sup_mn = [b_moire_sup_mn, [i; j]];  % 第一行是b_moire_r的系数，第二行是b_moire_l的系数
                b_moire_sup_vecs = [b_moire_sup_vecs, b_moire_sup];  % 第一行是x方向的分量，第二行是y方向的分量
            end
        end
    end
    
    % 将第2到第q_trunc+1个矢量删掉，因为之后做对称化会有重复
    % 删的话要从后往前删，因为指标会变化
    for i = (q_trunc+1):-1:2
        b_moire_sup_mn(:, i) = [];
        b_moire_sup_vecs(:, i) = [];
    end
    
    num_q = size(b_moire_sup_vecs, 2);  % 截断之后需要考虑的(m,n)的个数，实际作图时只需考虑低能的十几条能带即可

    % 然后做对称化：因为具有六重对称性
    Rotation_C6 = zeros(2);
    Rotation_C6(1,1) = 1/2;
    Rotation_C6(1,2) = -sqrt(3)/2;
    Rotation_C6(2,1) = sqrt(3)/2;
    Rotation_C6(2,2) = 1/2;

    b_mat = zeros(2);
    b_mat(:, 1) = b_moire_r;
    b_mat(:, 2) = b_moire_l;
    b_mat_inv = inv(b_mat);
    Rotation_C6_mn = b_mat_inv * Rotation_C6 * b_mat; % Rotation_C6_mn矩阵描述的是(m,n)系数的变换方式，Rotation_C6矩阵描述的是xy方向的变量的变换方式 

    for j = 0:4  % 共有5/6需要计算
        b_moire_sup_mn_temp = [0;0];
        b_moire_sup_vec_temp = [0.0;0.0];

        for i = (2+j*(num_q-1)):(num_q+j*(num_q-1))
            b_moire_sup_mn_temp = round(Rotation_C6_mn * b_moire_sup_mn(:,i));
            % b_moire_sup_mn_temp = Rotation_C6_mn * b_moire_sup_mn(:,i);
            b_moire_sup_mn = [b_moire_sup_mn, b_moire_sup_mn_temp];
            b_moire_sup_vec_temp = Rotation_C6 * b_moire_sup_vecs(:,i);
            b_moire_sup_vecs = [b_moire_sup_vecs, b_moire_sup_vec_temp];
        end
    end
    % 至此，b_moire_sup_vecs包含了所有在截断范围内的q,
    % b_moire_sup_mn则包含了这些q在b_moire_l和b_moire_r下的系数(都是整数)
    b_moire_sup_vecs = b_moire_norm * b_moire_sup_vecs;
end