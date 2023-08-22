function [ak_len_array, akx_array, aky_array] = brillouin_k_line(k_vertice_list, np_side_list, ak_norm)
    % k_vertices存放的是高对称折线上的顶点
    % np_sides存放的是两个相邻的顶点之间平均分为几等分
    % ak_norm的计算公式为 ：ak_norm = 8 * pi * abs(sin(theta / 2)) / 3; 无量纲
    
    ak_len = 0.0;
    %% 生成沿着高对称线上的所有k点
    % akx_array = []; % 无量纲
    % aky_array = []; % 无量纲
    % ak_len_array = [];
    
    akx_array = [k_vertice_list(1,1)]; % 无量纲
    aky_array = [k_vertice_list(2,1)]; % 无量纲
    ak_len_array = [0.0]; % 无量纲
    
    n_vertice = size(np_side_list, 2);
    for i = 1:n_vertice
        num_k = np_side_list(i);
        akx_array_temp = linspace(k_vertice_list(1,i), k_vertice_list(1,i+1), num_k);
        aky_array_temp = linspace(k_vertice_list(2,i), k_vertice_list(2,i+1), num_k);
        
        delta_kx = k_vertice_list(1,i+1) - k_vertice_list(1,i);
        delta_ky = k_vertice_list(2,i+1) - k_vertice_list(2,i);
        
        ak_len_temp = ak_len + ak_norm * sqrt(delta_kx^2 + delta_ky^2);
        ak_len_array_temp = linspace(ak_len, ak_len_temp, num_k);
        
        ak_len_array = [ak_len_array, ak_len_array_temp(2:end)];
        akx_array = [akx_array, akx_array_temp(2:end)];
        aky_array = [aky_array, aky_array_temp(2:end)];
        
        ak_len = ak_len_temp;
    end
    
    akx_array = akx_array * ak_norm;
    aky_array = aky_array * ak_norm;
end

%% check
% k_vertice_list = [[0; 0],[sqrt(3)/2; 1/2], [0; 1/2], [0; 1]];
% np_side_list = [101, 101, 101];

%% test code
% % 从K到Gamma:(0,0) \rightarrow (sqrt(3)/2, 1/2)
% num_k1 = 101;
% akx1 = linspace(0, sqrt(3)/2, num_k1);
% aky1 = linspace(0, 1/2, num_k1);
% ak1_len = ak_norm;
% 
% akx_array = [akx_array, akx1];
% aky_array = [aky_array, aky1];
% ak_len_temp = linspace(0, ak1_len, num_k1);
% ak_len_array = [ak_len_array, ak_len_temp];
% 
% % 从Gamma到M:(sqrt(3)/2, 1/2) \rightarrow (0, 1/2)
% num_k2 = 101;
% akx2 = linspace(sqrt(3)/2, 0, num_k2);
% aky2 = linspace(1/2, 1/2, num_k2);
% ak2_len = ak1_len + sqrt(3)/2 * ak_norm;
% 
% akx_array = [akx_array, akx2(2:end)];
% aky_array = [aky_array, aky2(2:end)];
% ak_len_temp = linspace(ak1_len, ak2_len, num_k2);
% ak_len_array = [ak_len_array, ak_len_temp(2:end)];
% 
% % 从M到Kp:(0, 1/2) \rightarrow (0, 1)
% num_k3 = 101;
% akx3 = linspace(0, 0, num_k3);
% aky3 = linspace(1/2, 1, num_k3);
% ak3_len = ak2_len + 1/2 * ak_norm;
% 
% akx_array = [akx_array, akx3(2:end)];
% aky_array = [aky_array, aky3(2:end)];
% ak_len_temp = linspace(ak2_len, ak3_len, num_k3);
% ak_len_array = [ak_len_array, ak_len_temp(2:end)];
% 
% akx_array = ak_norm * akx_array;
% aky_array = ak_norm * aky_array;
% % 得到k点的数目
% num_k = size(akx_array, 2); 