 function X_nm = get_exchange_integral(n1, n2, m1, m2, varargin)
    % 计算交换作用的积分 : X^{ab}_{n1n2, m1m2}
    % lb : magnetic length
    X_nm = 0;
    if ~((n1 + m1 - n2 - m2) == 0)
        return
    end
    
    n_min = min(n1, n2);
    n_max = max(n1, n2);
    nd = abs(n1 - n2);
    
    m_min = min(m1, m2);
    m_max = max(m1, m2);
    md = abs(m1 - m2);
    
    coeff = (nd + md - 1) / 2;
    
    format long
    if nargin >= 5 % 计算interlayer exchange integral
        % disp(varargin)
        d_interlayer = varargin{1}; % 层间距离，已经被无量纲化了(除以magnetic length)，d / lb
        func = @(y) exp(-sqrt(2) * d_interlayer .* sqrt(y)) .* exp(-y) .* (y.^coeff) .* laguerreL(n_min, nd, y) .* laguerreL(m_min, md, y);
        % 当LL index > 10时, RelTol不能设为1e-10,而应该设到1e-7左右才行
        % 当LL index == 11时, RelTol应该设到1e-9左右才行
        % 当LL index == 12时, RelTol应该设到1e-8左右才行
        % 当LL index == 13时, RelTol应该设到1e-7左右才行
        % 当LL index == 14时, RelTol应该设到1e-6左右才行
        % 当LL index == 15时, RelTol应该设到1e-5左右才行
        X_nm = integral(func, 0, Inf, 'RelTol', 1e-8);
        % X_nm = integral(func, 0, Inf);
    else % 计算intralayer exchange integral
        % func = @(y) exp(-y) .* (y.^coeff) .* laguerreL(n_min, nd, y) .* laguerreL(m_min, md, y);
        % X_nm = integral(func, 0, Inf, 'RelTol', 1e-6);
        X_nm = I_func(n_min, m_min, coeff, nd, md);
    end
    
    X_nm = X_nm * (-1)^md * sqrt((factorial(n_min) * factorial(m_min)) / (factorial(n_max) * factorial(m_max))) * (1j)^(nd + md) / sqrt(2);
    
    % 对于层内的积分，如果直接做数值积分要花0.08s左右，而如果套公式则需要花不到0.01s，算出来的结果是一样的
end


function I_val = I_func(n, m, alpha, beta, gm)
    x1 = gamma(alpha + 1) * gamma(m + gm + 1) * gamma(beta - alpha + n);
    x2 = gamma(n + 1) * gamma(m + 1) * gamma(gm + 1) * gamma(beta - alpha);
    x3 = hypergeom([-m, alpha + 1, alpha - beta + 1], [gm + 1, alpha + 1 - beta - n], 1);
    
    I_val = x1 / x2 * x3;
end


%% test : 参考PRB 85，165139（2021）TABLE I
% intralyer at different (a, b, c)
% X_nm = get_exchange_integral(0, 0, 0, 0)

% indice_list = zeros(3 * 3 * 3, 3);
% exchange_integrals_list = zeros(3 * 3 * 3, 3);
% count = 0;
% d = 0.335;
% tic
% for a = 0:2
%     for b = 0:2
%         for c = 0:2
%             count = count + 1;
%             indice_list(count, 1) = a;
%             indice_list(count, 2) = b;
%             indice_list(count, 3) = c;
%             
%             % 约定n2 >= n1, m2 <= m1（算这种情况已经完备了）
%             n2 = c;
%             n1 = c + a;
%             m1 = b;
%             m2 = b + a;
%             
%             % intralayer (Any B)
%             X_nm = get_exchange_integral(n1, n2, m1, m2);
%             exchange_integrals_list(count, 1) = X_nm;
%             
%             % Interlayer(B = 20T)
%             d_interlayer = 2 * d / 25.6 * sqrt(20); % 乘以2是因为这里计算的是第一层和第三层之间的exchange 
%             X_nm = get_exchange_integral(n1, n2, m1, m2, d_interlayer);
%             exchange_integrals_list(count, 2) = X_nm;
%             
%             % Interlayer(B = 30T)
%             d_interlayer = 2 * d / 25.6 * sqrt(30); % 乘以2是因为这里计算的是第一层和第三层之间的exchange
%             X_nm = get_exchange_integral(n1, n2, m1, m2, d_interlayer);
%             exchange_integrals_list(count, 3) = X_nm;
%             
%         end
%     end
% end
% toc
% 
% T = array2table([indice_list, exchange_integrals_list], 'VariableNames',{'a','b','c', 'Intralayer(Any B)','Interlayer(B = 20T)','Interlayer(B = 30T)'});


%% 草稿
% alpha = 0.5;
% m = 2;
% beta = 1;
% n = 2;
% gm = 1;
% 
% func = @(x) exp(-x) .* (x.^alpha) .* laguerreL(m, beta, x) .* laguerreL(n, gm, x)
% q = integral(func,0,Inf)
% 
% gamma(alpha + 1) * gamma(n + gm + 1) * gamma(beta - alpha + m) / (gamma(m + 1) * gamma(n + 1) * gamma(gm + 1) * gamma(beta - alpha))
% 
% hypergeom([-n, alpha + 1, alpha - beta + 1], [gm + 1, alpha + 1 - beta - m], 1)
