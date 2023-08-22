%% 计算电容
%% 法二:沿Vt或Vb一个方向分别取两条平行的截线
fprintf(1,'step1 : get the slope of Rxx peaks at cnp\n')

dim_1 = size(lockin_Rxx, 1);
dim_2 = size(lockin_Rxx, 2);

Rxx_cnp_peak_pos_vecs = zeros(1, dim_1);
max_val = max(max(lockin_Rxx));
for i = 1:dim_1
    index_vecs = find(lockin_Rxx(i,:) == max_val);
    mid_point_index = ceil(mean(index_vecs));
    Rxx_cnp_peak_pos_vecs(i) = gate_2(mid_point_index);
end

% 得到这条线的slope
slope_vecs = zeros(1, dim_1-1);
for i = 1:(dim_1-1)
    slope_vecs(i) = (gate_1(i+1) - gate_1(i)) / (Rxx_cnp_peak_pos_vecs(i+1) - Rxx_cnp_peak_pos_vecs(i));
end
slope = mean(slope_vecs);

slope = 1.8146;  % 通过手动拉线得到

fig1 = figure;
im1 = imagesc(gate_2, gate_1, lockin_Rxx);  % gate_1作为y轴, gate_2作为x轴, size(lockin_Rxx, 1) = length(gate_1), size(lockin_Rxx, 2) = length(gate_2)
axis([v2_min v2_max v1_min v1_max]) % 前面两个是x轴的区间，后面两个是y轴的区间
box on
shading interp
hold on
set(gca, 'YDir', 'normal')
xlabel('$V_2$','interpreter','latex', 'FontSize', 12);
ylabel('$V_1$','interpreter','latex', 'FontSize', 12);
colormap(gca, clc_map)
colorbar
title('$R_{xx}(\frac{h}{e^2})$ raw','interpreter','latex', 'FontSize', 14);
hold on
plot(Rxx_cnp_peak_pos_vecs, gate_1, 'm*', 'LineWidth', 0.025, 'MarkerSize',10)

select_index = input('select index is : ');  % 85
fig0 = figure;
set(gcf,'position',[250 300 1000 300])
% 可以把峰直接在图中画出来
plot(gate_2, lockin_Rxx(select_index,:))
findpeaks(lockin_Rxx(select_index,:))
%
[peaks, locs] = findpeaks(lockin_Rxx(select_index,:));
gate_at_peaks = gate_2(locs);
[~,pos] = max(peaks);
peak_pos1 = locs(pos);
peak_V2 = gate_2(peak_pos1);
diff_gate_at_peaks = diff(gate_at_peaks);
delta_V2 = median(diff_gate_at_peaks);
format long
C2 = carrier_density_oneflux * 4 / delta_V2;  % 4来自于自旋*谷的四重简并度
C2 = C2 / 10^(12);

C1 = C2 / abs(slope);

[gate_2_xx, gate_1_yy] = meshgrid(gate_2, gate_1);
[carrier_density, disp_field] = convert_VtVb_2_nD(gate_1_yy, gate_2_xx, C1, C2);

% convert to nD map
figure
subplot(121)
im = pcolor(carrier_density, disp_field, lockin_Rxx);
grid off
axis([min(min(carrier_density)) max(max(carrier_density)) min(min(disp_field)) max(max(disp_field))])
box on
shading interp
hold on
set(gca, 'YDir', 'normal')
xlabel('$n(cm^{-2})$','interpreter','latex', 'FontSize', 12);
ylabel('$D(V/nm)$','interpreter','latex', 'FontSize', 12);
colormap(gca, clc_map)
colorbar
title('$R_{xx}(n,D)$','interpreter','latex', 'FontSize', 14);

subplot(122)
% log plot
lockin_Rxx_log = log10(abs(lockin_Rxx));
im2 = pcolor(carrier_density, disp_field, lockin_Rxx_log);
grid off
axis([min(min(carrier_density)) max(max(carrier_density)) min(min(disp_field)) max(max(disp_field))])
box on
shading interp
hold on
set(gca, 'YDir', 'normal')
xlabel('$n(cm^{-2})$','interpreter','latex', 'FontSize', 12);
ylabel('$D(V/nm)$','interpreter','latex', 'FontSize', 12);
colormap(gca, jet)
colorbar
title('$R_{xx}(n,D)$ log','interpreter','latex', 'FontSize', 14);