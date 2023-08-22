addpath('D:\matlab\graphene-package\multilayer_graphene\utils')
addpath('D:\matlab\graphene-package\multilayer_graphene\trilayer_ABA_graphene\trilayer_ABA_hamiltonian\')

% effective model of trilayer graphene
% hopping_params = output_value_array(1,:);
% gamma0 = hopping_params(1);
% gamma1 = hopping_params(2);
% gamma2 = hopping_params(3);
% gamma3 = hopping_params(4);  
% gamma4 = hopping_params(5);
% gamma5 = hopping_params(6);
% delta = hopping_params(7);
% Delta2 = hopping_params(8);

gamma0 = 3.1;
gamma1 = 0.390;
gamma2 = -0.028;
gamma3 = 0.315;  
gamma4 = 0.041;
gamma5 = 0.050;
delta = 0.046;

% Delta1 = 0.025; % 0 / 25 / 50 / 150 / 180 / 250
Delta2 = 0;

Delta1_list = [0.025, 0.025, 0.025];
ene_upper_bound = 0.03;
ene_lower_bound = - 0.03;
sigma = 0.0006;
num_ene_point = 201;
ene_list = linspace(ene_lower_bound, ene_upper_bound, num_ene_point);
    
figure;
for ii = 1:3
    Delta1 = Delta1_list(ii); % Delta1 = 0.0000 / 0.0172 / 0.0344

    %% band calculation
    tic
    dim_x = 201;
    dim_y = 201;
    akx_list = linspace(-0.075, 0.075, dim_x);
    aky_list = linspace(-0.075, 0.075, dim_y);
    [akx_mesh, aky_mesh] = meshgrid(akx_list, aky_list);

    eig_enes_K = zeros(dim_y, dim_x, 6);
    eig_enes_Kp = zeros(dim_y, dim_x, 6);

    for i = 1:dim_x
        for j = 1:dim_y
            akx = akx_mesh(j, i);
            aky = aky_mesh(j, i);

            [HK_m_ham, HKp_m_ham, HK_b_ham, HKp_b_ham] = construct_trilayer_ABA_LL_six_bands_continuum_model(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, delta, Delta2, akx, aky);
            if Delta1 == 0
                [~, D] = eig(HK_m_ham);
                eig_enes_K(j, i, 1:2) = diag(D);

                [~, D] = eig(HK_b_ham);
                eig_enes_K(j, i, 3:6) = diag(D);

                [~, D] = eig(HKp_m_ham);
                eig_enes_Kp(j, i, 1:2) = diag(D);

                [~, D] = eig(HKp_b_ham);
                eig_enes_Kp(j, i, 3:6) = diag(D);

            else
                Delta1_mat = zeros(2, 4);
                Delta1_mat(1, 1) = Delta1;
                Delta1_mat(2, 2) = Delta1;

                HK_ham = [HK_m_ham, Delta1_mat; Delta1_mat', HK_b_ham];
                HKp_ham = [HKp_m_ham, Delta1_mat; Delta1_mat', HKp_b_ham];

                [~, D] = eig(HK_ham);
                eig_enes_K(j, i, :) = diag(D);

                [~, D] = eig(HKp_ham);
                eig_enes_Kp(j, i, :) = diag(D);
            end

        end
    end
    toc

    %% 作图
    center_x = round((1 + dim_x) / 2);
    center_y = round((1 + dim_y) / 2);
    subplot(3, 2, 2*(ii - 1) + 1);
    for i = 1:6
        plot(akx_list, eig_enes_K(:,center_y,i), 'b', 'LineWidth', 2)
        hold on;
        plot(akx_list, eig_enes_Kp(:,center_y,i),'r--', 'LineWidth', 2)
        % plot(aky_mesh(:, center_x), eig_enes_K(center_x,:,i),'--','LineWidth', 2)
        hold on
    end
    grid on;
    xlim([-0.10,0.10]) % 无量纲量ka在-0.1到0.1之间
    ylim([ene_lower_bound, ene_upper_bound]) % 30meV
    xlabel('ka \rightarrow');
    ylabel('E(\bf{k}) in eV \rightarrow');
    xticks([-1 -1/2 0 1/2 1]);
    xticklabels({'-1' '-1/2' '0' '1/2' '1'})
    % title(['\Delta_1 = ', num2str(Delta1), ': for k_y = 0 at K valley'])
    hold on

    subplot(3, 2, 2*(ii - 1) + 2);
    dos_K_list = get_dos_list(eig_enes_K, ene_list, sigma);
    dos_Kp_list = get_dos_list(eig_enes_Kp, ene_list, sigma);
    dos_list = dos_K_list + dos_Kp_list;
    plot(ene_list, dos_list, 'b', 'LineWidth', 2)
    grid on;
    title(['DOS at \Delta_1 = ', num2str(Delta1)])
    view(90,270)
    
%     for i = 1:6
%         plot(aky_list, eig_enes_K(center_x,:,i), 'b', 'LineWidth', 2)
%         hold on;
%         plot(aky_list, eig_enes_Kp(center_x,:,i),'--','LineWidth', 2)
%         hold on
%     end
%     grid on;
%     xlim([-0.10,0.10]) % 无量纲量ka在-0.1到0.1之间
%     ylim([ene_lower_bound, ene_upper_bound]) % 30meV
%     xlabel('ka \rightarrow');
%     ylabel('E(\bf{k}) in eV \rightarrow');
%     xticks([-1 -1/2 0 1/2 1]);
%     xticklabels({'-1' '-1/2' '0' '1/2' '1'})
%     title(['\Delta_1 = ', num2str(Delta1), ': for k_y = 0 at Kp valley'])
%     hold on
    
end

sigma = 0.0008;
num_ene_point = 201;
ene_list = linspace(ene_lower_bound, ene_upper_bound, num_ene_point);
dos_K_list = get_dos_list(eig_enes_K, ene_list, sigma);
dos_Kp_list = get_dos_list(eig_enes_Kp, ene_list, sigma);
dos_list = dos_K_list + dos_Kp_list;
plot(ene_list, dos_list, 'b', 'LineWidth', 2)
grid on;
title(['DOS at \Delta_1 = ', num2str(Delta1)])
view(90,270)