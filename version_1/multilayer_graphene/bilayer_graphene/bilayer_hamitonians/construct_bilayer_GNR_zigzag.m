function [H0, H1, H_1, y_coords_layer1, y_coords_layer2] = construct_bilayer_GNR_zigzag(gamma_params, onsite_params, N_rows, y_shifts, BS_BS0)
    gamma0 = gamma_params(1);
    gamma1 = gamma_params(2); % from dimer site to dimer site
    gamma3 = gamma_params(3); % from non-dimer site to non-dimer site
    gamma4 = gamma_params(4); % from dimer site to non-dimer site
    
    onsite_A1 = onsite_params(1);
    onsite_B1 = onsite_params(2);
    onsite_A2 = onsite_params(3);
    onsite_B2 = onsite_params(4);
    
    % 一般 y_shift 写成 1 / 3
    % 不考虑on-site potential
    [H0_layer1, H1_layer1, ~, y_coords_layer1] = construct_monolayer_GNR_zigzag(gamma0, [onsite_A1, onsite_B1], N_rows, y_shifts(1), BS_BS0);
    [H0_layer2, H1_layer2, ~, y_coords_layer2] = construct_monolayer_GNR_zigzag(gamma0, [onsite_A2, onsite_B2], N_rows, y_shifts(2), BS_BS0);
    % helper_check_hermite(H0_layer1, 1e-8)
    % helper_check_hermite(H0_layer2, 1e-8)
    
    N_block = floor(N_rows / 4);
    half_block = mod(N_rows, 4);
    
    H0 = zeros(N_rows * 2, N_rows * 2);  % 1：dims是layer1, (dims+1):2*dims是layer2
    H1 = zeros(N_rows * 2, N_rows * 2);
    H_1 = zeros(N_rows * 2, N_rows * 2);
    Im = 1j;
    
    % % construct the full H0
    % layer subspace
    H0(1:N_rows, 1:N_rows) = H0_layer1;
    H0((N_rows+1):2*N_rows, (N_rows+1):2*N_rows) = H0_layer2;
    
    % layer1_2_layer2
    % from layer1_(n,j) to layer2_(n,j)
    for n = 1:N_block
        phi44 = BS_BS0 * (y_coords_layer2(4*(n-1) + 4) + y_coords_layer1(4*(n-1) + 4)) / 2; %第一个指标是layer2的，第一个指标是layer1的
        phi33 = BS_BS0 * (y_coords_layer2(4*(n-1) + 3) + y_coords_layer1(4*(n-1) + 3)) / 2;
        phi23 = BS_BS0 * (y_coords_layer2(4*(n-1) + 2) + y_coords_layer1(4*(n-1) + 3)) / 2;
        phi22 = BS_BS0 * (y_coords_layer2(4*(n-1) + 2) + y_coords_layer1(4*(n-1) + 2)) / 2;
        phi11 = BS_BS0 * (y_coords_layer2(4*(n-1) + 1) + y_coords_layer1(4*(n-1) + 1)) / 2;
        
        % from layer1_4 to layer2_4: dimer site to non-dimer site
        H0(N_rows + 4*(n-1) + 4, 4*(n-1) + 4) = gamma4 * exp(-Im * phi44);
        H0(4*(n-1) + 4, N_rows + 4*(n-1) + 4) = conj(H0(N_rows + 4*(n-1) + 4, 4*(n-1) + 4));
        
        % from layer1_4 to layer2_3: dimer site to dimer site
        H0(N_rows + 4*(n-1) + 3, 4*(n-1) + 4) = gamma1;
        H0(4*(n-1) + 4, N_rows + 4*(n-1) + 3) = gamma1;
        
        % from layer1_4 to layer2_2: dimer site to non-dimer site
        H0(N_rows + 4*(n-1) + 3, 4*(n-1) + 4) = gamma4;
        H0(4*(n-1) + 4, N_rows + 4*(n-1) + 3) = gamma4;
        
        % from layer1_3 to layer2_3: non-dimer site to dimer site
        H0(N_rows + 4*(n-1) + 3, 4*(n-1) + 3) = gamma4 * exp(-Im * phi33);
        H0(4*(n-1) + 3, N_rows + 4*(n-1) + 3) = conj(H0(N_rows + 4*(n-1) + 3, 4*(n-1) + 3));
        
        % from layer1_3 to layer2_2: non-dimer site to non-dimer site
        H0(N_rows + 4*(n-1) + 2, 4*(n-1) + 3) = gamma3 * exp(-Im * phi23);
        H0(4*(n-1) + 3, N_rows + 4*(n-1) + 2) = conj(H0(N_rows + 4*(n-1) + 2, 4*(n-1) + 3));
        
        % from layer1_2 to layer2_2: dimer site to non-dimer site
        H0(N_rows + 4*(n-1) + 2, 4*(n-1) + 2) = gamma4 * exp(-Im * phi22);
        H0(4*(n-1) + 2, N_rows + 4*(n-1) + 2) = conj(H0(N_rows + 4*(n-1) + 2, 4*(n-1) + 2));
        
        % from layer1_1 to layer2_2: non-dimer site to non-dimer site
        H0(N_rows + 4*(n-1) + 2, 4*(n-1) + 1) = gamma3;
        H0(4*(n-1) + 1, N_rows + 4*(n-1) + 2) = conj(H0(N_rows + 4*(n-1) + 2, 4*(n-1) + 1));
        
        % from layer1_1 to layer2_1: non-dimer site to dimer site
        H0(N_rows + 4*(n-1) + 1, 4*(n-1) + 1) = gamma4 * exp(-Im * phi11);
        H0(4*(n-1) + 1, N_rows + 4*(n-1) + 1) = conj(H0(N_rows + 4*(n-1) + 1, 4*(n-1) + 1));
    end
    
    % from layer1_(n,j+1) to layer2_(n,j)
    for n = 1:(N_block - 1)
        phi41 = BS_BS0 * (y_coords_layer2(4*(n-1) + 4) + y_coords_layer1(4*n + 1)) / 2;
        
        % from layer1_1 to layer2_4: non-dimer site to non-dimer site
        H0(N_rows + 4*(n-1) + 4, 4*n + 1) = gamma3 * exp(-Im * phi41);
        H0(4*n + 1, N_rows + 4*(n-1) + 4) = conj(H0(N_rows + 4*(n-1) + 4, 4*n + 1));
        
        % from layer1_1 to layer2_3: non-dimer site to dimer-site
        H0(N_rows + 4*(n-1) + 3, 4*n + 1) = gamma4;
        H0(4*n + 1, N_rows + 4*(n-1) + 3) = conj(H0(N_rows + 4*(n-1) + 3, 4*n + 1));
    end
    % helper_check_hermite(H0, 1e-8)
    
    % % construct the full H1
    % single layer subspace
    H1(1:N_rows, 1:N_rows) = H1_layer1; % layer1 subspace
    H1((N_rows + 1):2*N_rows, (N_rows + 1):2*N_rows) = H1_layer2; % layer2 subspace
    
    % from layer1_(n+1,j) to layer2_(n,j)
    for n = 1:N_block
        phi44 = BS_BS0 * (y_coords_layer2(4*(n-1) + 4) + y_coords_layer1(4*(n - 1) + 4)) / 2;
        % from layer1_4 to layer2_4: dimer site to non-dimer site
        H1(N_rows + 4*(n-1) + 4, 4*(n-1) + 4) = gamma4 * exp(Im * phi44);
        
        % from layer1_3 to layer2_4: non-dimer site to non-dimer site
        H1(N_rows + 4*(n-1) + 4, 4*(n-1) + 3) = gamma3;
        
        % from layer1_3 to layer2_3: non-dimer site to dimer site
        phi33 = BS_BS0 * (y_coords_layer2(4*(n-1) + 3) + y_coords_layer1(4*(n - 1) + 3)) / 2;
        H1(N_rows + 4*(n-1) + 3, 4*(n-1) + 3) = gamma4 * exp(Im * phi33);
        
        % from layer1_3 to layer2_2: non-dimer site to non-dimer site
        phi23 = BS_BS0 * (y_coords_layer2(4*(n-1) + 3) + y_coords_layer1(4*(n - 1) + 2)) / 2;
        H1(N_rows + 4*(n-1) + 2, 4*(n-1) + 3) = gamma3 * exp(Im * phi23);
        
        % from layer1_3 to layer2_1: non-dimer site to dimer site
        H1(N_rows + 4*(n-1) + 1, 4*(n-1) + 3) = gamma4;
        
        % from layer1_2 to layer2_2: dimer site to non-dimer site
        phi22 = BS_BS0 * (y_coords_layer2(4*(n-1) + 3) + y_coords_layer1(4*(n - 1) + 2)) / 2;
        H1(N_rows + 4*(n-1) + 2, 4*(n-1) + 2) = gamma4 * exp(Im * phi22);
        
        % from layer1_2 to layer2_1: dimer site to dimer site
        H1(N_rows + 4*(n-1) + 1, 4*(n-1) + 2) = gamma1;
        
        % from layer1_1 to layer2_1: non-dimer site to dimer site
        phi11 = BS_BS0 * (y_coords_layer2(4*(n-1) + 1) + y_coords_layer1(4*(n - 1) + 1)) / 2;
        H1(N_rows + 4*(n-1)+1, 4*(n-1)+1) = gamma4 * exp(Im * phi11);
    end
    
    % from layer1_(n+1,j+1) to layer2_(n,j)
    for n = 1:(N_block - 1)
        phi41 = BS_BS0 * (y_coords_layer2(4*(n-1) + 4) + y_coords_layer1(4*n + 1)) / 2;
        % from layer1_1 to layer2_4: non-dimer site to non-dimer site
        H1(N_rows + 4*(n-1) + 4, 4*n+1) = gamma3 * exp(Im * phi41);
        
        % from layer1_2 to layer2_4: dimer site to non-dimer site
        H1(N_rows + 4*(n-1) + 4, 4*n+2) = gamma4;
    end
    
    H_1 = H1';
end