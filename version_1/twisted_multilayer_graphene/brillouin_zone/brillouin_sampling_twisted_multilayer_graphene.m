function [bz_bound, bz_inner, k_K_list] = brillouin_sampling_twisted_multilayer_graphene(np_side, ak_norm)
    point_list1 = [];
    point_list2 = [];
    vector1 = [sqrt(3)/2;-1/2];
    vector2 = [0;1];
    for i = 0:(np_side-1)
        vec_temp = i / (np_side-1) * vector1;
        point_list1 = [point_list1, vec_temp];

        vec_temp = i / (np_side-1) * vector2;
        point_list2 = [point_list2, vec_temp];
    end

    %% 得到布里渊区的内部的点，这是可以通过旋转120°/240°得到其它的点
    bz_inner = [];
    for i = 2:(np_side - 1)
        bz_inner = [bz_inner, point_list2(:, i) + point_list1(:, 1:np_side - 1)];
    end

    %% 不放心的话，可以画一下
    % scatter(bz_bound(1,:), bz_bound(2,:),50,'r','p','filled');
    % hold on
    % scatter(bz_inner(1,:), bz_inner(2,:),30,'b','s');

    %% 旋转120°和240°得到其它的内部点
    Rotation_C3 = zeros(2);
    Rotation_C3(1,1) = -1/2;
    Rotation_C3(1,2) = -sqrt(3)/2;
    Rotation_C3(2,1) = sqrt(3)/2;
    Rotation_C3(2,2) = -1/2;

    np_inner = size(bz_inner, 2);
    % 逆时针旋转120°
    for i = 1:np_inner
        vec_temp = Rotation_C3 * bz_inner(:, i);
        bz_inner = [bz_inner, vec_temp];
    end

    % 逆时针旋转240°
    for i = 1:np_inner
        vec_temp = Rotation_C3 * bz_inner(:, i + np_inner);
        bz_inner = [bz_inner, vec_temp];
    end
    
    %% 别忘了Gamma点（0，0）
    bz_inner = [bz_inner, [0;0]];

    %% 找出K和K'点
    k_K_list = [];
    k_K = [sqrt(3)/2; -1/2];
    k_Kp = [sqrt(3)/2; 1/2];
    k_K_list = [k_K_list, k_K, k_Kp, Rotation_C3 * k_K, Rotation_C3 * k_Kp, Rotation_C3 * Rotation_C3 * k_K, Rotation_C3 * Rotation_C3 * k_Kp];
    
    %% 得到布里渊区的边界，我们只需要三边即可，其它三边是重复的
    bz_bound = [];
    for i = 1:6
        if i == 6
            vec_x_list = linspace(k_K_list(1,i), k_K_list(1,1), np_side);
            vec_y_list = linspace(k_K_list(2,i), k_K_list(2,1), np_side);
        else
            vec_x_list = linspace(k_K_list(1,i), k_K_list(1,i+1), np_side);
            vec_y_list = linspace(k_K_list(2,i), k_K_list(2,i+1), np_side);
        end
        bz_bound = [bz_bound, [vec_x_list(1:end-1); vec_y_list(1:end-1)]];
    end
    
    %% 再乘以ak_norm
    k_K_list = k_K_list * ak_norm;
    bz_inner = bz_inner * ak_norm;
    bz_bound = bz_bound * ak_norm;
end

%% test code
% Gamma
% k_Gamma = [0; 0];
% np_side = 10;
% 
% point_list1 = [];
% point_list2 = [];
% vector1 = [-sqrt(3)/2; -1/2];
% vector2 = [0; 1];
% for i = 0:(np_side-1)
%     vec_temp = i / (np_side-1) * vector1;
%     point_list1 = [point_list1, vec_temp];
%     
%     vec_temp = i / (np_side-1) * vector2;
%     point_list2 = [point_list2, vec_temp];
% end
% 
% % 得到布里渊区的边界，我们只需要两边即可，其它四边是重复的
% bz_bound = [];
% bz_bound = [bz_bound, point_list1(:, end) + point_list2];
% for i = (np_side - 1):-1:2
%     bz_bound = [bz_bound, point_list1(:, i) + point_list2(:, end)];
% end
% 
% % 得到布里渊区的内部的点，这是可以通过旋转120°/240°得到其它的点
% bz_inner = [];
% for i = 2:(np_side - 1)
%     bz_inner = [bz_inner, point_list1(:, i) + point_list2(:, 1:np_side - 1)];
% end
% 
% % 不放心的话，可以画一下
% scatter(bz_bound(1,:), bz_bound(2,:),50,'r','p','filled');
% hold on
% scatter(bz_inner(1,:), bz_inner(2,:),30,'b','s');
% 
% % 旋转120°和240°得到其它的内部点
% Rotation_C3 = zeros(2);
% Rotation_C3(1,1) = -1/2;
% Rotation_C3(1,2) = -sqrt(3)/2;
% Rotation_C3(2,1) = sqrt(3)/2;
% Rotation_C3(2,2) = -1/2;
% 
% np_inner = size(bz_inner, 2);
% 逆时针旋转120°
% for i = 1:np_inner
%     vec_temp = Rotation_C3 * bz_inner(:, i);
%     bz_inner = [bz_inner, vec_temp];
% end
% 
% 逆时针旋转240°
% for i = 1:np_inner
%     vec_temp = Rotation_C3 * bz_inner(:, i + np_inner);
%     bz_inner = [bz_inner, vec_temp];
% end
% 
% % 别忘了Gamma点（0，0）
% bz_inner = [bz_inner, [0;0]];
% 
% % 找出K和K'点
% k_K_list = [];
% k_K = [-sqrt(3)/2; 1/2];
% k_Kp = [-sqrt(3)/2; -1/2];
% k_K_list = [k_K_list, k_K, k_Kp, Rotation_C3 * k_K, Rotation_C3 * k_Kp, Rotation_C3 * Rotation_C3 * k_K, Rotation_C3 * Rotation_C3 * k_Kp];
% 
% % 但是对于转角石墨烯，最后还要再平移一个矢量
% k_K_list = k_K_list + [sqrt(3)/2;1/2];
% bz_inner = bz_inner + [sqrt(3)/2;1/2];
% bz_bound = bz_bound + [sqrt(3)/2;1/2];
% 
% % 不放心的话，可以画一下
% scatter(bz_bound(1,:), bz_bound(2,:),50,'r','p','filled');
% hold on
% scatter(bz_inner(1,1:end-1), bz_inner(2,1:end-1),30,'b','s');
% hold on
% scatter(bz_inner(1,end), bz_inner(2,end),30,'k','+');
% hold on
% plot(k_K_list(1,:), k_K_list(2,:), 'k--');
% plot([k_K_list(1,end), k_K_list(1,1)], [k_K_list(2,end), k_K_list(2,1)], 'k--');