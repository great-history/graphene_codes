function [exchange_integrals_list, indice_list, indice_cell] = calc_exchange_integrals(LL_index_max, kd_interlayer_list, flag_parfor)
    % LL_index_max : 计算交换积分时最高的LL index
    % kd_interlayer_list中分别是intralayer, nearest interlayer, next-nearest interlayer
    % flag_parfor : 是否进行并行计算
    num_LL = LL_index_max + 1;
    num_interlayer = length(kd_interlayer_list);
    
    format long
    % 构造exchange integral
    indice_list = []; % 存放(a, b, c), 计算结果(a,b,c)和(a,c,b)得到的结果是一样的，因此只需要算其中一个即可，我采取的notation是b>=c
    indice_cell = {};
    % exchange_integrals_cell = {};
    count = 0;
    for a = 0:LL_index_max
        for b = 0:LL_index_max
            if (b + a) > LL_index_max
                break
            end

            for c = 0:b
                % 对于确定的(a,b,c)，一般共有4种可能的(n1,n2,m1,m2)
                % n1 = c + a;
                % n2 = c;
                % m1 = b;
                % m2 = b + a;
                
                count = count + 1;
                indice_list = [indice_list; a, b, c];
                indice_cell{a + 1}{b + 1}{c + 1} = count;
                % exchange_integrals_cell{a + 1}{b + 1}{c + 1} = zeros(num_interlayer, 1);
            end
        end
    end
    exchange_integrals_list = zeros(count, num_interlayer);

    % 开启多核
    % poolobj = gcp('nocreate');
    % if isempty(poolobj)
    %     disp('启动并行运算，核心数：6');
    %     % Perform a basic check by entering this code, where "local" is one kind of cluster profile.
    %     parpool('local', 6);
    % else
    %     disp(['并行运算已启动，核心数：' num2str(poolobj.NumWorkers)]);
    % end

    tic
    if flag_parfor
        parfor ii = 1:count
            a = indice_list(ii, 1);
            b = indice_list(ii, 2);
            c = indice_list(ii, 3);
            
            % 算一种情况就足够了
            n1 = c + a;
            n2 = c;
            m1 = b;
            m2 = b + a;
            
            % exchange_integrals = zeros(num_interlayer, 1);
            % intralayer (Any B)
            % X_nm = get_exchange_integral(n1, n2, m1, m2);
            % exchange_integrals_list(ii, 1) = X_nm;
            
            % Interlayer @ fixed B
            for jj = 1:num_interlayer
                X_nm = get_exchange_integral(n1, n2, m1, m2, kd_interlayer_list(jj));
                exchange_integrals_list(ii, jj) = X_nm;
                % exchange_integrals(jj) = X_nm;
            end
            % exchange_integrals_cell{a + 1}{b + 1}{c + 1} = exchange_integrals;
        end
    else
        for ii = 1:count
            a = indice_list(ii, 1);
            b = indice_list(ii, 2);
            c = indice_list(ii, 3);
            
            % 算一种情况就足够了
            n1 = c + a;
            n2 = c;
            m1 = b;
            m2 = b + a;
            
            % intralayer (Any B)
            % X_nm = get_exchange_integral(n1, n2, m1, m2);
            % exchange_integrals_list(ii, 1) = X_nm;
            
            % Interlayer @ fixed B
            for jj = 1:num_interlayer
                X_nm = get_exchange_integral(n1, n2, m1, m2, kd_interlayer_list(jj));
                exchange_integrals_list(ii, jj) = X_nm;
            end
        end
    end
    toc

end