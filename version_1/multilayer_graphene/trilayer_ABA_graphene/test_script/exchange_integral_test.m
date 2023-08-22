LL_index_max = 12;
num_LL = LL_index_max + 1;

B_field = 1.25;
d_interlayer = 0.335;
d_interlayer_list = [0, d_interlayer, 2 * d_interlayer]; % 以nm为单位
kd_interlayer_list = 2 * d_interlayer_list / 25.6 * sqrt(B_field); % 无量纲值

flag_calc = true; % 是否重新计算并保存数据,否则load原来计算的数据
if flag_calc % 对于LL_index_max = 12, 使用并行计算历时 66.482162 秒
    % 开启多核
    flag_parfor = true;
    if flag_parfor
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            disp('启动并行运算，核心数：6');
            % Perform a basic check by entering this code, where "local" is one kind of cluster profile.
            parpool('local', 6);
        else
            disp(['并行运算已启动，核心数：' num2str(poolobj.NumWorkers)]);
        end
    end

    [exchange_integrals_list, indice_list, indice_cell] = get_exchange_integrals_list(LL_index_max, kd_interlayer_list, flag_parfor);
    file_path = ['.\trilayer_ABA_data\exchange_integrals_B_', num2str(roundn(B_field, -2)), 'T_LL_index_max_', num2str(LL_index_max), '.mat'];
    save(file_path, 'exchange_integrals_list', 'indice_list', 'indice_cell');
else
    file_path = ['.\trilayer_ABA_data\exchange_integrals_B_', num2str(roundn(B_field, -2)), 'T_LL_index_max_', num2str(LL_index_max), '.mat'];
    load(file_path);
end

