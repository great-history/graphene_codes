function [rho_layer1, rho_layer2] = get_layer_density_four_bands(eig_vecs, valley, LL_index_cutoff)
    % LL_index_cutoff = 10 说明最高LL的指标为10，即只考虑LL_index = 0,1,2,3,4,5,6,7,8,,9,10，维数为dims = 2 * LL_index_cutoff
    % valley : + 1 for K valley ; -1 for Kp valley
    % in K valley, the Landau levels in the order : [|0,B2>, |1,B2>, |0,B1>, |0,A2> // |2,B2>, |0,A1>, |1,B1>, |1,A2> // |3,B2>, |1,A1>, |2,B1>, |2,A2> // ... ...]
    % in Kp valley, the Landau levels in the order : [|0,A1>, |1,A1>, |0,A2>, |0,B1> // |2,A1>, |0,B2>, |1,A2>, |1,B1> // |3,A1>, |1,B2>, |2,A2>, |2,B1> // ... ...]
    
    rho_layer1 = 0.0;
    rho_layer2 = 0.0;
    
    if valley == 1 % at K valley
        for ii = 1:size(eig_vecs, 2)
            eig_vec = eig_vecs(:, ii);
            rho_layer1 = rho_layer1 + (abs(eig_vec(4)))^2;
            rho_layer2 = rho_layer2 + (abs(eig_vec(1)))^2 + (abs(eig_vec(2)))^2 + (abs(eig_vec(3)))^2;

            for n = 2:LL_index_cutoff
                start_index = 4 * (n - 1);
                rho_layer1 = rho_layer1 + (abs(eig_vec(start_index + 2)))^2 + (abs(eig_vec(start_index + 4)))^2;
                rho_layer2 = rho_layer2 + (abs(eig_vec(start_index + 1)))^2 + (abs(eig_vec(start_index + 3)))^2;
            end
        end
        
    else % at Kp valley
        for ii = 1:size(eig_vecs, 2)
            eig_vec = eig_vecs(:, ii);
            rho_layer1 = rho_layer1 + (abs(eig_vec(1)))^2 + (abs(eig_vec(2)))^2 + (abs(eig_vec(3)))^2;
            rho_layer2 = rho_layer2 + (abs(eig_vec(4)))^2;

            for n = 2:LL_index_cutoff
                start_index = 4 * (n - 1);
                rho_layer1 = rho_layer1 + (abs(eig_vec(start_index + 1)))^2 + (abs(eig_vec(start_index + 3)))^2;
                rho_layer2 = rho_layer2 + (abs(eig_vec(start_index + 2)))^2 + (abs(eig_vec(start_index + 4)))^2;
            end
        end
    end
end