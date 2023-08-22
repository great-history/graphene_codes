function [H0, H1, H_1, y_coords] = construct_monolayer_GNR_zigzag(gamma0, mass, N_rows, y_shift, BS_BS0)
    %% add magnetic field on zigzag
    % a = 0.246; % nm
    % length_ratio = (a / 25.66)^2;
    
    dims = N_rows;
    N_block = floor(N_rows / 4);
    half_block = mod(N_rows, 4);
    
    % 生成y坐标(沿着ribbon宽度方向)
    y_coords = zeros(dims, 1);
    for n = 1:N_block
        y_coords(4*(n-1) + 1) = (2 * n - 5 / 3) + y_shift;
        y_coords(4*(n-1) + 2) = (2 * n - 4 / 3) + y_shift;
        y_coords(4*(n-1) + 3) = (2 * n - 2 / 3) + y_shift;
        y_coords(4*(n-1) + 4) = (2 * n - 1 / 3) + y_shift;
    end

    % construct the tight-binding model under a magnetic field
    Im = 1j;
    H0 = zeros(dims, dims);
    H1 = zeros(dims, dims);
    for n = 1:N_block
        phi1 = BS_BS0 * (y_coords(4*(n-1) + 4) + y_coords(4*(n-1) + 3)) / 2;
        phi2 = BS_BS0 * (y_coords(4*(n-1) + 2) + y_coords(4*(n-1) + 1)) / 2;

        H0(4*(n-1)+4, 4*(n-1)+4) = mass(2);
        H0(4*(n-1)+3, 4*(n-1)+3) = mass(1);
        H0(4*(n-1)+2, 4*(n-1)+2) = mass(2);
        H0(4*(n-1)+1, 4*(n-1)+1) = mass(1);

        H0(4*(n-1)+3, 4*(n-1)+4) = gamma0 * exp(Im * phi1);
        H0(4*(n-1)+4, 4*(n-1)+3) = conj(H0(4*(n-1)+3, 4*(n-1)+4));

        H0(4*(n-1)+2, 4*(n-1)+3) = gamma0;
        H0(4*(n-1)+3, 4*(n-1)+2) = conj(H0(4*(n-1)+2, 4*(n-1)+3));

        H0(4*(n-1)+1, 4*(n-1)+2) = gamma0 * exp(-Im * phi2);
        H0(4*(n-1)+2, 4*(n-1)+1) = conj(H0(4*(n-1)+1, 4*(n-1)+2));

        if ~(n == N_block)
            H0(4*n + 1, 4*(n-1) + 4) = gamma0;
            H0(4*(n-1) + 4, 4*n + 1) = conj(H0(4*n + 1, 4*(n-1)+4));
        end

        % H1不是厄密的矩阵
        H1(4*(n-1)+4, 4*(n-1)+3) = gamma0 * exp(Im * phi1); 
        H1(4*(n-1)+1, 4*(n-1)+2) = gamma0 * exp(Im * phi2);

        % H_{-1}是H1的厄密，故不用创建
    end
    
    H_1 = H1'; % H_{-1}是H1的厄密
end