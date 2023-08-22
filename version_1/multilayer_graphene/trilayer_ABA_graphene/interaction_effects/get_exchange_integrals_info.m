% 得到与相互作用有关的系数
function [exchange_integrals_list, indice_list, indice_cell] = get_exchange_integrals_info(flag_load, flag_parfor, save_path)
    % flag_load 导入数据还是计算数据    flag_parfor是否开启并行
    if ~(flag_load)
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
        save(save_path, 'exchange_integrals_list', 'indice_list', 'indice_cell');
    else
        [exchange_integrals_list, indice_list, indice_cell] = load(save_path);
    end

end