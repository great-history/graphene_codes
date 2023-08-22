
n_terminal = 6;

% current probes
current_in = 1;
current_out = 4;
I = zeros(n_terminal, 1);
I(current_in) = 1;
I(current_out) = -1;

% transmatrix:
% % ����Quantum Spin Hall Effect
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

% % ����Quantum Parity Hall Effect
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

% ����Quantum Hall Effect������ͨ����
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

%% ����˹��Ԫ��
% % �ӵ�һ�е����һ��
vec_temp = zeros(1, n_terminal+1);
for i = 1:(n_terminal-1)
    % ���G(i,i) == 0��Ҫ����һ����
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
    
    % ���G(i,i)��ȻΪ0��ֱ����������(i+1)
    if G(i, i) == 0
        continue
    end
    
    % ��i�����ĸ�˹��Ԫ
    for j = (i+1):n_terminal
        if ~(G(j,i) == 0)
            ratio = G(j,i) / G(i,i);
            G(j,:) = G(j,:) - ratio * G(i,:);
        end
    end
end

% % �����һ�е���һ��
for i = n_terminal:-1:2
    % ���G(i,i) == 0������ֱ����һ��
    if abs(G(i,i)) < 1e-8
        continue
    end
    
    % �ȶԵ�i�н��й�һ��
    G(i, :) = G(i, :) / G(i, i);
    
    % ��i�����ĸ�˹��Ԫ
    for j = (i-1):-1:1
        if ~(G(j,i) == 0)
            ratio = G(j,i);
            G(j,:) = G(j,:) - ratio * G(i,:);
        end
    end
end
% ���Ե�1�н��й�һ��
G(1, :) = G(1, :) / G(1, 1);


%% ����絼
% ����two-terminal �� four-terminal conductance
% % voltage probes
v_left = 2;
v_right = 3;

ratio = G(v_right, end-1) / G(v_left, end-1);
if abs(abs(ratio)-1) > 1e-8
    disp("some problems")
end

delta_G = G(v_right, end) - G(v_left, end);