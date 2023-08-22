%% get dos list
% common parameters
type = "linear"; % 选取DOS函数的类型
w_bandwidth = 1; % 将能带宽度取为1，所有能量尺度以其为单位
ene_min = - 3 * w_bandwidth;
ene_max = 3 * w_bandwidth;
ene_points = 1201;
ene_list = linspace(ene_min, ene_max, ene_points);

% different types of DOS
% Linear DOS
dos_list = zeros(1, ene_points);
for ii = 1:ene_points
    ene_current = ene_list(ii);
    if abs(ene_current) < w_bandwidth
        dos_list(ii) = 2 * abs(ene_current) / (w_bandwidth)^2;
    else
        dos_list(ii) = 0;
    end
end

% Effect of Van Hove Singularity

% Effect of asymmetry in the conduction & valence flat band

% DOS fropm the continuum oBistritzer-MacDonald model near magic angle

%% figure plot
current_fig = figure;
hold on

p1 = plot(ene_list, dos_list);
current_legend = type + '-' + 'DOS';
legend(current_legend, 'Location','North');

% saveas(gcf, save_path); %保存当前窗口的图像