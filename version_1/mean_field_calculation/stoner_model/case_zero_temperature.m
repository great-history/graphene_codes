% Dirac Revivals at zero temperature
% common parameters
w_bandwidth = 1.0;
type = "linear";
n_flavor = 4;
Hubbard_U = 1.2 * w_bandwidth;

% global chemical potential
num_chem_pot_points = 801;
chem_pot_list = linspace(-4*w_bandwidth, 4*w_bandwidth, num_chem_pot_points);

% set up fmincon
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     disp('启动并行运算，核心数：6');
%     % Perform a basic check by entering this code, where "local" is one kind of cluster profile.
%     parpool('local', 6);
% else
%     disp(['并行运算已启动，核心数：' num2str(poolobj.NumWorkers)]);
% end

% for each global chemical potential
for ii = 1:num_chem_pot_points
    chem_pot = chem_pot_list(ii); % global chemical potential
    
    kinetic_flavor_list = kinetic_energy_func(type, ene_list, dos_list, ene_num_points, w_bandwidth);
    % mf_grand_pot = mean_field_grand_potential_func(kinetic_flavor_list, chem_pot_flavor_list, Hubbard_U, chem_pot, w_bandwidth, type);
    
    %% chemical_pot_flavor is the object for optimization
    % options = optimoptions('fmincon', 'UseParallel', true); % 使用并行
    % [output_value_list, fval,exitflag,output,lambda,grad,hessian] = fmincon(@mf_grand_pot, initial_value_list, A, b, Aeq, beq, lb_list, ub_list);
end