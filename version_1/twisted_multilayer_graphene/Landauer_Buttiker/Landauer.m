
n_terminal = 6;

% current probes
current_in = 1;
current_out = 4;
I = zeros(n_terminal, 1);
I(current_in) = 1;
I(current_out) = -1;

% transmatrix:
% % 对于Quantum Spin Hall Effect
% T = zeros(n_terminal);
% for i = 1:n_terminal
%     T(i,i) = 2;
%     
%     left = rem(i-1, n_terminal);
%     right = rem(i+1, n_terminal);
%     
%     if left == 0
%         T(i, left + n_terminal) = -1;
%     else
%         T(i, left) = -1;
%     end
%     
%     if right == 0
%         T(i, right + n_terminal) = -1;
%     else
%         T(i, right) = -1;
%     end
% end

% % 对于Quantum Parity Hall Effect
% T = zeros(n_terminal);
% for i = 1:n_terminal
%     T(i,i) = 4;
%     
%     left = rem(i-1, n_terminal);
%     right = rem(i+1, n_terminal);
%     
%     if left == 0
%         T(i, left + n_terminal) = -2;
%     else
%         T(i, left) = -2;
%     end
%     
%     if right == 0
%         T(i, right + n_terminal) = -2;
%     else
%         T(i, right) = -2;
%     end
% end

% 对于Quantum Hall Effect（两个通道）
T = zeros(n_terminal);
for i = 1:n_terminal
    T(i,i) = 2;
    
    next = rem(5+i, n_terminal);
    
    if next == 0
        T(i, next + n_terminal) = -2;
    else
        T(i, next) = -2;
    end
end


G = [T,I];

%% 做高斯消元法
% % 从第一行到最后一行
vec_temp = zeros(1, n_terminal+1);
for i = 1:(n_terminal-1)
    % 如果G(i,i) == 0则要交换一下行
    if G(i,i) == 0
        for ii = (i+1):n_terminal
            if ~(G(ii, i) == 0)
                vec_temp = G(ii, :);
                G(ii, :) = G(i, :);
                G(i, :) = vec_temp;
                break
            end
        end
    end
    
    % 如果G(i,i)仍然为0则直接跳过进行(i+1)
    if G(i, i) == 0
        continue
    end
    
    % 第i分量的高斯消元
    for j = (i+1):n_terminal
        if ~(G(j,i) == 0)
            ratio = G(j,i) / G(i,i);
            G(j,:) = G(j,:) - ratio * G(i,:);
        end
    end
end

% % 从最后一行到第一行
for i = n_terminal:-1:2
    % 如果G(i,i) == 0则跳过直接下一步
    if abs(G(i,i)) < 1e-8
        continue
    end
    
    % 先对第i行进行归一化
    G(i, :) = G(i, :) / G(i, i);
    
    % 第i分量的高斯消元
    for j = (i-1):-1:1
        if ~(G(j,i) == 0)
            ratio = G(j,i);
            G(j,:) = G(j,:) - ratio * G(i,:);
        end
    end
end
% 最后对第1行进行归一化
G(1, :) = G(1, :) / G(1, 1);


%% 计算电导
% 计算two-terminal 或 four-terminal conductance
% % voltage probes
v_left = 2;
v_right = 3;

ratio = G(v_right, end-1) / G(v_left, end-1);
if abs(abs(ratio)-1) > 1e-8
    disp("some problems")
end

delta_G = G(v_right, end) - G(v_left, end);