% function save_data_as_txt(save_path, flag_all_bands)
%     
% end

%% test code
tbg_path = "D:\matlab\graphene-package\twisted_multilayer_graphene\save_data\tbg_data";
tbg_dir = tbg_path + "\" + datestr(now,'yyyy-mm-dd_HH-MM-SS');
if ~exist(tbg_dir,'dir')
	mkdir(tbg_dir);
end

%% 输入参数的保存：跃迁参数 / k点数组 / 低能能带数 / 等等
parameter_filename = tbg_dir + "\tbg_parameters.txt";
fid=fopen(parameter_filename,'wt');%写入文件路径

name_cell = foo(onsite_a, onsite_b, gamma0, w_aa, w_ab, theta, q_trunc);
for i = 1:length(name_cell)
    fprintf(fid, "%-16s", name_cell{i});
end
fprintf(fid, "\n");

for i = 1:length(name_cell)
    fprintf(fid, "%-16.4f", eval(name_cell{i}));
end
fprintf(fid, "\n");

fclose(fid);

%% 保存k点数组
k_point_filename = tbg_dir + "\tbg_k_point.mat";
save(k_point_filename, 'k_vertice_list', 'np_side_list', 'akx_array', 'aky_array', 'ak_len_array', 'bm_sup_mn', 'bm_sup_vecs');

%% 本征能量的保存(所有能带)    
flag_all_bands = true;
if flag_all_bands
    all_bands_filename = tbg_dir + "\tbg_all_bands.mat";
    save(all_bands_filename, 'eig_vals_K', 'eig_vals_Kp');
end

%% 本征能量/本征态的保存(低能能带)
low_ene_filename = tbg_dir + "\tbg_low_ene.mat";
save(low_ene_filename, 'eig_vals_K_low_ene', 'eig_vecs_K_low_ene', 'eig_vals_Kp_low_ene', 'eig_vecs_Kp_low_ene');

%% 得到变量名的函数
function name_cell = foo(varargin)
    name_cell = cell(nargin,1);
    for k = 1:1:nargin
        name = inputname(k);
        name_cell{k} = name;
    end
end