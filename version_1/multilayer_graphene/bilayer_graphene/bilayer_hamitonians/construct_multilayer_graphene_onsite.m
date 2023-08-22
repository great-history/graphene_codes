function H_ham = construct_multilayer_graphene_onsite(U, gamma_n, akx, aky, num_layer, H_ham, type)
    % U是加的电场，gamma_n是次紧邻跃迁, type是哈密顿量的类型——是紧束缚还是连续模型
    u_list = linspace(U, -U, num_layer);
    
    if type == "continuum_K"  % 连续模型 @ valley K        
        gk = - 3;
    elseif type == "continuum_Kp"  % 连续模型 @ valley Kp
        gk = - 3;
    else % 紧束缚模型
        a1_vec = [sqrt(3)/2; 3/2];
        a2_vec = [sqrt(3)/2;-3/2];
        ak_vec = [akx; aky];
        gk = 2 * cos(dot(ak_vec, a1_vec)) + 2 * cos(dot(ak_vec, a2_vec)) + 2 * cos(dot(ak_vec, a1_vec + a2_vec));
    end
    
    
    for layer_index = 1:num_layer
        u_temp = u_list(layer_index);
        
        H_ham((layer_index - 1) * 2 + 1, (layer_index - 1) * 2 + 1) = u_temp + gamma_n * gk;
        H_ham(layer_index * 2, layer_index * 2) = u_temp + gamma_n * gk;
        
    end
end