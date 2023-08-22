%% 计算涉及到的所有exchange integral
% 该函数用来计算不同LL之间的exchange integrals
% 对于ABC三层石墨烯而言，其LLL具有特殊的性质：每个谷都发生了层极化，并且LL_index n = 0,1,2 ———— 总共有 2 * 3 * 3 = 12个态，即存在十二重简并，这个子空间也被称为duodectet subspace
% 因此构造S矩阵只需考虑n = 0,1,2这三个LL之间的层内或层间的交换相互作用：在PRB 85，165139（2021）TABLE I中已经给出了在任意磁场下的层内交换相互作用和B=20T/30T下的层间交换相互作用
% 对于ABA三层石墨烯而言，其LLL子空间也被称为duodectet subspace，(单层1*2 + 双层2*2) * 2 = 12维
% 因此构造S矩阵只需考虑n = 0,1,2这三个LL之间的层内或层间的交换相互作用：在PRB 85，165139（2021）TABLE I中已经给出了在任意磁场下的层内交换相互作用和B=20T/30T下的层间交换相互作用    
function [exchange_integrals_intralayer_cell, exchange_integrals_interlayer_cell] = get_all_exchange_integrals(LL_index_min, LL_index_max, d_interlayer, mag_length)
    % LL_index_min = 0 ; LL_index_max = 2
    % 下面就是构造出这样的S矩阵
    format long
    % 构造exchange integral
    tic
    num_LL = LL_index_max - LL_index_min + 1;
    %% 对(n1,n2)进行编码
    indice_list = zeros(num_LL * num_LL, 2); % 将(n1,n2)按照三进制进行编码
    for ii = 0:(num_LL - 1)
        num1 = ii * num_LL;
        for jj = 0:(num_LL - 1)
            num2 = jj * 1;
            indice = num1 + num2 + 1;
            indice_list(indice, 1) = ii; % 三进制数的第一位(位权为3)
            indice_list(indice, 2) = jj; % 三进制数的第二位(位权为1)
        end
    end
    
    %% 计算层间和层内的交换相互作用
    exchange_integrals_intralayer_cell = cell(num_LL * num_LL, 1); % 其实只需要三个指标即可，因为 X_{n1 n2, m1 m2}不为零需要满足条件：
    exchange_integrals_interlayer_cell = cell(num_LL * num_LL, 1); % 其实只需要三个指标即可，因为 X_{n1 n2, m1 m2}不为零需要满足条件：
    
    % 一种可以减少计算量的方法就是对指标进行编码，然后利用exchange integral的对称性关系(计算量可以减半)
    dims = size(indice_list, 1);
    kd_interlayer = d_interlayer / mag_length; % 计算层间交换相互作用时需要用到
    for ii = dims:-1:1 %% 对应max([n1 n2]_3, [m1 m2]_3)，其中[n1 n2]_3和[m1 m2]_3是三进制数
        exchange_integrals_intralayer_list = zeros(ii, 1); %% 存放min([n1 n2]_3, [m1 m2]_3)
        exchange_integrals_interlayer_list = zeros(ii, 1); %% 存放min([n1 n2]_3, [m1 m2]_3)

        for jj = 1:ii
            n1 = indice_list(ii, 1);
            n2 = indice_list(ii, 2);
            m1 = indice_list(jj, 1);
            m2 = indice_list(jj, 2);
            
            % intralayer (Any B)
            X_nm_intralayer = get_exchange_integral(n1, n2, m1, m2);

            % Interlayer(B = 20T)
            X_nm_interlayer = get_exchange_integral(n1, n2, m1, m2, kd_interlayer);

            % % Interlayer(B = 30T)
            % kd_interlayer = 2 * d_interlayer / 25.6 * sqrt(30); % 乘以2是因为这里计算的是第一层和第三层之间的exchange
            % X_nm_interlayer = get_exchange_integral(n1, n2, m1, m2, d_interlayer);

            % 将交换积分存入list之中
            exchange_integrals_intralayer_list(jj) = X_nm_intralayer;
            exchange_integrals_interlayer_list(jj) = X_nm_interlayer;
        end

        % 将list存入cell之中
        exchange_integrals_intralayer_cell{ii} = exchange_integrals_intralayer_list; %% 存放max([n1 n2]_3, [m1 m2]_3)
        exchange_integrals_interlayer_cell{ii} = exchange_integrals_interlayer_list; %% 存放max([n1 n2]_3, [m1 m2]_3)
    end
    toc
end