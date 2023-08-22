function [H0, H1, H_1, y_coords_layer1, y_coords_layer2, y_coords_layer3] = construct_trilayer_GNR_zigzag(type, gamma_params, diagonal_terms, N_rows, yshifts, BS_BS0)
    % 使用两种方法进行构造:一种是实空间，另一种是通过mirror symmetry
    % onsite_A1 = Delta1 + Delta2; onsite_B1 = Delta1 + Delta2 + delta;
    % onsite_A3 = - Delta1 + Delta2; onsite_B3 = -Delta1 + Delta2 + delta;
    % onsite_A2 = -2 * Delta2 + delta; onsite_B2 = -2 * Delta2
    if type == "lattice"
        [H0, H1, H_1, y_coords_layer1, y_coords_layer2, y_coords_layer3] = construct_trilayer_GNR_zigzag_lattice(gamma_params, diagonal_terms, N_rows, yshifts, BS_BS0);
    elseif type == "symmetry"
        [H0, H1, H_1, y_coords_layer1, y_coords_layer2, y_coords_layer3] = construct_trilayer_GNR_zigzag_symmetry(gamma_params, diagonal_terms, N_rows, yshifts, BS_BS0);
    end
end

function [H0, H1, H_1, y_coords_layer1, y_coords_layer2, y_coords_layer3] = construct_trilayer_GNR_zigzag_lattice(gamma_params, onsite_params, N_rows, yshifts, BS_BS0)
    
    % 第一种方法
    gamma0 = gamma_params(1);
    gamma1 = gamma_params(2);
    gamma2 = gamma_params(3);
    gamma3 = gamma_params(4);
    gamma4 = gamma_params(5);
    gamma5 = gamma_params(6);
    
    onsite_A1 = onsite_params(1);
    onsite_B1 = onsite_params(2);
    onsite_A2 = onsite_params(3);
    onsite_B2 = onsite_params(4);
    onsite_A3 = onsite_params(5);
    onsite_B3 = onsite_params(6);
    
    dims = N_rows * 3;
    H0 = zeros(dims, dims);
    H1 = zeros(dims, dims);
    H_1 = zeros(dims, dims);
    
    Im = 1j;
    % layer1 与 layer2 组成的bilayer graphene
    gamma_params_bilayer = [gamma0, gamma1, gamma3, gamma4];
    onsite_params_bilayer = [onsite_A1, onsite_B1, onsite_A2, onsite_B2];
    y_shifts_bilayer = [yshifts(1), yshifts(2)];
    [H0_layer12, H1_layer12, ~, y_coords_layer1, y_coords_layer2] = construct_bilayer_GNR_zigzag(gamma_params_bilayer, onsite_params_bilayer, N_rows, y_shifts_bilayer, BS_BS0);
    % 1:N_rows是layer_1 subspace, (N_rows+1):2*N_rows是layer_1 subspace
    
    gamma_params_bilayer = [gamma0, gamma1, gamma3, gamma4];
    onsite_params_bilayer = [onsite_A3, onsite_B3, onsite_A2, onsite_B2];
    y_shifts_bilayer = [yshifts(3), yshifts(2)];
    [H0_layer32, H1_layer32, ~, y_coords_layer3, y_coords_layer2] = construct_bilayer_GNR_zigzag(gamma_params_bilayer, onsite_params_bilayer, N_rows, y_shifts_bilayer, BS_BS0);
    % 1:N_rows是layer_3 subspace, (N_rows+1):2*N_rows是layer_2 subspace
    
    % % construct H0
    H0(1:2*N_rows, 1:2*N_rows) = H0_layer12;
    H0((2*N_rows + 1):3*N_rows, (2*N_rows + 1):3*N_rows) = H0_layer32(1:N_rows, 1:N_rows);
    % from layer3 to layer2
    H0((N_rows + 1):2*N_rows, (2*N_rows + 1):3*N_rows) = H0_layer32((N_rows + 1):2*N_rows, 1:N_rows);
    % from layer2 to layer3
    H0((2*N_rows + 1):3*N_rows, (N_rows + 1):2*N_rows) = H0_layer32(1:N_rows, (N_rows + 1):2*N_rows);
    
    % from layer1 to layer3: dimer to dimer gamma5; non-dimer to non-dimer gamma2
    for n = 1:N_block
        H0((n-1)*4 + 4, 2*N_rows + (n-1)*4 + 4) = gamma5;
        H0(2*N_rows + (n-1)*4 + 4, (n-1)*4 + 4) = gamma5;
        
        H0((n-1)*4 + 2, 2*N_rows + (n-1)*4 + 2) = gamma5;
        H0(2*N_rows + (n-1)*4 + 2, (n-1)*4 + 2) = gamma5;
        
        H0((n-1)*4 + 3, 2*N_rows + (n-1)*4 + 3) = gamma2;
        H0(2*N_rows + (n-1)*4 + 3, (n-1)*4 + 3) = gamma2;
        
        H0((n-1)*4 + 1, 2*N_rows + (n-1)*4 + 1) = gamma2;
        H0(2*N_rows + (n-1)*4 + 1, (n-1)*4 + 1) = gamma2;
    end
end

function [H0, H1, H_1, y_coords_layer_b1, y_coords_layer_b2, y_coords_layer_m] = construct_trilayer_GNR_zigzag_symmetry(gamma_params, delta_params, N_rows, yshifts, BS_BS0)
    gamma0 = gamma_params(1);
    gamma1 = gamma_params(2);
    gamma2 = gamma_params(3);
    gamma3 = gamma_params(4);
    gamma4 = gamma_params(5);
    gamma5 = gamma_params(6);
    
    % 第二种方法
    Delta1 = delta_params(1);
    Delta2 = delta_params(2);
    delta = delta_params(3);
    
    dims = N_rows * 3;
    H0 = zeros(dims, dims);
    H1 = zeros(dims, dims);
    H_1 = zeros(dims, dims);
    
    Im = 1j;
    
    % basis: 1 / sqrt(2)(|A1> + |A3>); 1 / sqrt(2)(|B1> + |B3>); |A2>; |B2>
    gamma_params_bilayer = [gamma0, sqrt(2) * gamma1, sqrt(2) * gamma3, sqrt(2) * gamma4];
    onsite_params_bilayer = [Delta2 + gamma2 / 2, Delta2 + delta + gamma5 / 2, -2 * Delta2 + delta, -2 * Delta2];
    yshifts_bilayer = [yshifts(1), yshifts(2)];
    [H0_b, H1_b, ~, y_coords_layer_b1, y_coords_layer_b2] = construct_bilayer_GNR_zigzag(gamma_params_bilayer, onsite_params_bilayer, N_rows, yshifts_bilayer, BS_BS0);
    
    % basis: 1 / sqrt(2)(|A1> - |A3>); 1 / sqrt(2)(|B1> - |B3>)
    onsite_params_monolayer = [Delta2 - gamma2 / 2, Delta2 + delta - gamma5 / 2];
    yshifts_monolayer = 0.0;
    [H0_m, H1_m, ~, y_coords_layer_m] = construct_monolayer_GNR_zigzag(gamma0, onsite_params_monolayer, N_rows, yshifts_monolayer, BS_BS0);
    
    H0(1:N_rows, 1:N_rows) = H0_m;
    H0((N_rows + 1):3*N_rows, (N_rows + 1):3*N_rows) = H0_b;
    H1(1:N_rows, 1:N_rows) = H1_m;
    H1((N_rows + 1):3*N_rows, (N_rows + 1):3*N_rows) = H1_b;
    
    % monolayer_A 与 bilayer_A1存在on_site的跃迁;  monolayer_B 与 bilayer_B1存在on_site的跃迁
    H0(1:N_rows, (N_rows + 1):2 * N_rows) = Delta1 * eye(N_rows);
    H0((N_rows + 1):2 * N_rows, 1:N_rows) = Delta1 * eye(N_rows);
    
    H_1 = H1';
end


% 并行计算
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     disp('启动并行运算，核心数：8');
%     % Perform a basic check by entering this code, where "local" is one kind of cluster profile.
%     parpool('local', 8);
% else
%     disp(['并行运o算已启动，核心数：' num2str(poolobj.NumWorkers)]);
% end