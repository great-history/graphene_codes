[file_cell, path] = uigetfile('D:\pycharm_projects\Nano_Lab_App_new\test_data\sample26\*.dat','MultiSelect','on');

if iscell(file_cell) == 0
    if file_cell == 0
        return;
    end
end

name_cell = cell(length(file_cell),1);
for i = 1:length(file_cell)
    name = file_cell{i};
    name = name(1:end-4); % 最后四个是.dat
    name_cell{i} = name;
end

data_cell = cell(length(file_cell),1);
for i = 1:length(file_cell)
    data = readmatrix([path, file_cell{i}]); %读取数据构成矩阵
    if isnan(data(1))  % data(1)应该是一个index，如果是NaN，那么说明无数据
        data_cell{i} = data;
        continue
    end
    
    if strfind(name_cell{i}, "forward")
        n_row = 1;
        for j = 2:size(data, 1)
            if data(j, 1) == 0
                break
            else
                n_row = n_row + 1;
            end
        end
    elseif strfind(name_cell{i}, "backward")
        n_row = 0;
        for j = 1:size(data, 1)
            n_row = n_row + 1;
            if data(j, 1) == 0
                break
            end
        end
    else
        continue
    end
    data = reshape(data, [n_row, size(data, 1)/n_row, size(data, 2)]);
    data_cell{i} = data;
end

% lockin : index               X               Y               R           Theta
% source : index            read           error            leak