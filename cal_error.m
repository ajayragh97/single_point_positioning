clc;
clear all;
close all;
% working directory
path = pwd;
% add functions
addpath([path,'/04_functions']);


eqw = load('03_data/EQW_short.mat');
eqw_rdm = load('03_data/EQW+RDM_short.mat');
elv01 = load('03_data/ELV01_short.mat');
elv01_rdm = load('03_data/ELV01+RDM_short.mat');
cn_light = load('03_data/CN_short_light.mat');
cn_rdm_light = load('03_data/CN+RDM_short_light.mat');
cn_heavy = load('03_data/CN_short_heavy.mat');
cn_rdm_heavy = load('03_data/CN+RDM_short_heavy.mat');
elvcn_50 = load('03_data/ELVCN_short_thresh_50.mat');
elvcn_60 = load('03_data/ELVCN_short_thresh_60.mat');
elvcn_rdm_50 = load('03_data/ELVCN+RDM_short_thresh_50.mat');
elvcn_rdm_60 = load('03_data/ELVCN+RDM_short_thresh_60.mat');

start = 1;
stop = size(eqw.user_pos, 1); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EQW %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_eqw, mean_error_ver_eqw, rmse_hor_eqw, ...
 rmse_ver_eqw, max_error_hor_eqw, max_error_ver_eqw]= ...
            calculate_error(eqw.user_pos(start:stop, 1), eqw.user_pos(start:stop, 2), eqw.user_pos(start:stop, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EQW + RDM %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_eqw_rdm, mean_error_ver_eqw_rdm, rmse_hor_eqw_rdm, ...
 rmse_ver_eqw_rdm, max_error_hor_eqw_rdm, max_error_ver_eqw_rdm]= ...
            calculate_error(eqw_rdm.user_pos(start:stop, 1), eqw_rdm.user_pos(start:stop, 2), eqw_rdm.user_pos(start:stop, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elv01, mean_error_ver_elv01, rmse_hor_elv01, ...
 rmse_ver_elv01, max_error_hor_elv01, max_error_ver_elv01]= ...
            calculate_error(elv01.user_pos(start:stop, 1), elv01.user_pos(start:stop, 2), elv01.user_pos(start:stop, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV + RDM %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elv01_rdm, mean_error_ver_elv01_rdm, rmse_hor_elv01_rdm, ...
 rmse_ver_elv01_rdm, max_error_hor_elv01_rdm, max_error_ver_elv01_rdm]= ...
            calculate_error(elv01_rdm.user_pos(start:stop, 1), elv01_rdm.user_pos(start:stop, 2), elv01_rdm.user_pos(start:stop, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CN light %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_cn_li, mean_error_ver_cn_li, rmse_hor_cn_li, ...
 rmse_ver_cn_li, max_error_hor_cn_li, max_error_ver_cn_li]= ...
            calculate_error(cn_light.user_pos(start:stop, 1), cn_light.user_pos(start:stop, 2), cn_light.user_pos(start:stop, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CN + RDM light%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_cn_rdm_li, mean_error_ver_cn_rdm_li, rmse_hor_cn_rdm_li, ...
 rmse_ver_cn_rdm_li, max_error_hor_cn_rdm_li, max_error_ver_cn_rdm_li]= ...
            calculate_error(cn_rdm_light.user_pos(start:stop, 1), cn_rdm_light.user_pos(start:stop, 2), cn_rdm_light.user_pos(start:stop, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CN heavy%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_cn_hv, mean_error_ver_cn_hv, rmse_hor_cn_hv, ...
 rmse_ver_cn_hv, max_error_hor_cn_hv, max_error_ver_cn_hv]= ...
            calculate_error(cn_heavy.user_pos(start:stop, 1), cn_heavy.user_pos(start:stop, 2), cn_heavy.user_pos(start:stop, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CN + RDM heavy%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_cn_rdm_hv, mean_error_ver_cn_rdm_hv, rmse_hor_cn_rdm_hv, ...
 rmse_ver_cn_rdm_hv, max_error_hor_cn_rdm_hv, max_error_ver_cn_rdm_hv]= ...
            calculate_error(cn_rdm_heavy.user_pos(start:stop, 1), cn_rdm_heavy.user_pos(start:stop, 2), cn_rdm_heavy.user_pos(start:stop, 3));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV + CN 50 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elvcn_50, mean_error_ver_elvcn_50, rmse_hor_elvcn_50, ...
 rmse_ver_elvcn_50, max_error_hor_elvcn_50, max_error_ver_elvcn_50]= ...
            calculate_error(elvcn_50.user_pos(start:stop, 1), elvcn_50.user_pos(start:stop, 2), elvcn_50.user_pos(start:stop, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV + CN 60 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elvcn_60, mean_error_ver_elvcn_60, rmse_hor_elvcn_60, ...
 rmse_ver_elvcn_60, max_error_hor_elvcn_60, max_error_ver_elvcn_60]= ...
            calculate_error(elvcn_60.user_pos(start:stop, 1), elvcn_60.user_pos(start:stop, 2), elvcn_60.user_pos(start:stop, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV + CN + RDM 50%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elvcn_rdm_50, mean_error_ver_elvcn_rdm_50, rmse_hor_elvcn_rdm_50, ...
 rmse_ver_elvcn_rdm_50, max_error_hor_elvcn_rdm_50, max_error_ver_elvcn_rdm_50]= ...
            calculate_error(elvcn_rdm_50.user_pos(start:stop, 1), elvcn_rdm_50.user_pos(start:stop, 2), elvcn_rdm_50.user_pos(start:stop, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV + CN + RDM 60%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elvcn_rdm_60, mean_error_ver_elvcn_rdm_60, rmse_hor_elvcn_rdm_60, ...
 rmse_ver_elvcn_rdm_60, max_error_hor_elvcn_rdm_60, max_error_ver_elvcn_rdm_60]= ...
            calculate_error(elvcn_rdm_60.user_pos(start:stop, 1), elvcn_rdm_60.user_pos(start:stop, 2), elvcn_rdm_60.user_pos(start:stop, 3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Horizontal Max Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [ max_error_hor_elvcn_rdm_50, max_error_hor_elvcn_50, ...
                 max_error_hor_elvcn_rdm_60, max_error_hor_elvcn_60, ...
                 max_error_hor_cn_rdm_li, max_error_hor_cn_li, ...
                 max_error_hor_cn_rdm_hv, max_error_hor_cn_hv, ...
                 max_error_hor_elv01_rdm, max_error_hor_elv01,...
                 max_error_hor_eqw_rdm, max_error_hor_eqw];
                
error_types = {'ELVCN+RDM-50', 'ELVCN-50', ...
               'ELVCN+RDM-60','ELVCN-60', ...
               'CN+RDM-light', 'CN-light', ...
               'CN+RDM-heavy', 'CN-heavy', ...
               'ELV+RDM', 'ELV', ...
               'EQW+RDM', 'EQW'};
         

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
% xticklabels(error_types);
ylabel('Max Error');
title('Max Error - Horizontal (Short Sequence)');

for i = 1:numel(error_values)
    text(i, 0.4, error_types{i}, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vertical Max Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [ max_error_ver_elvcn_rdm_50, max_error_ver_elvcn_50, ...
                 max_error_ver_elvcn_rdm_60, max_error_ver_elvcn_60, ...
                 max_error_ver_cn_rdm_li, max_error_ver_cn_li, ...
                 max_error_ver_cn_rdm_hv, max_error_ver_cn_hv, ...
                 max_error_ver_elv01_rdm, max_error_ver_elv01,...
                 max_error_ver_eqw_rdm, max_error_ver_eqw];
                
error_types = {'ELVCN+RDM-50', 'ELVCN-50', ...
               'ELVCN+RDM-60','ELVCN-60', ...
               'CN+RDM-light', 'CN-light', ...
               'CN+RDM-heavy', 'CN-heavy', ...
               'ELV+RDM', 'ELV', ...
               'EQW+RDM', 'EQW'};

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
% xticklabels(error_types);
ylabel('Max Error');
title('Max Error - Vertical (Short Sequence)');

for i = 1:numel(error_values)
    text(i, 0.5, error_types{i}, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Horizontal Mean Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [ mean_error_hor_elvcn_rdm_50, mean_error_hor_elvcn_50, ...
                 mean_error_hor_elvcn_rdm_60, mean_error_hor_elvcn_60, ...
                 mean_error_hor_cn_rdm_li, mean_error_hor_cn_li, ...
                 mean_error_hor_cn_rdm_hv, mean_error_hor_cn_hv, ...
                 mean_error_hor_elv01_rdm, mean_error_hor_elv01,...
                 mean_error_hor_eqw_rdm, mean_error_hor_eqw];
                
error_types = {'ELVCN+RDM-50', 'ELVCN-50', ...
               'ELVCN+RDM-60','ELVCN-60', ...
               'CN+RDM-light', 'CN-light', ...
               'CN+RDM-heavy', 'CN-heavy', ...
               'ELV+RDM', 'ELV', ...
               'EQW+RDM', 'EQW'};

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
% xticklabels(error_types);
ylabel('Mean Error');
title('Mean Error - Horizontal (Short Sequence)');

for i = 1:numel(error_values)
    text(i, 0.15, error_types{i}, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vertical Mean Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [ mean_error_ver_elvcn_rdm_50, mean_error_ver_elvcn_50, ...
                 mean_error_ver_elvcn_rdm_60, mean_error_ver_elvcn_60, ...
                 mean_error_ver_cn_rdm_li, mean_error_ver_cn_li, ...
                 mean_error_ver_cn_rdm_hv, mean_error_ver_cn_hv, ...
                 mean_error_ver_elv01_rdm, mean_error_ver_elv01,...
                 mean_error_ver_eqw_rdm, mean_error_ver_eqw];
                
error_types = {'ELVCN+RDM-50', 'ELVCN-50', ...
               'ELVCN+RDM-60','ELVCN-60', ...
               'CN+RDM-light', 'CN-light', ...
               'CN+RDM-heavy', 'CN-heavy', ...
               'ELV+RDM', 'ELV', ...
               'EQW+RDM', 'EQW'};

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
% xticklabels(error_types);
ylabel('Mean Error');
title('Mean Error - Vertical (Short Sequence)');

for i = 1:numel(error_values)
    text(i, 0.2, error_types{i}, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Horizontal RMS Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [ rmse_hor_elvcn_rdm_50, rmse_hor_elvcn_50, ...
                 rmse_hor_elvcn_rdm_60, rmse_hor_elvcn_60, ...
                 rmse_hor_cn_rdm_li, rmse_hor_cn_li, ...
                 rmse_hor_cn_rdm_hv, rmse_hor_cn_hv, ...
                 rmse_hor_elv01_rdm, rmse_hor_elv01,...
                 rmse_hor_eqw_rdm, rmse_hor_eqw];
                
error_types = {'ELVCN+RDM-50', 'ELVCN-50', ...
               'ELVCN+RDM-60','ELVCN-60', ...
               'CN+RDM-light', 'CN-light', ...
               'CN+RDM-heavy', 'CN-heavy', ...
               'ELV+RDM', 'ELV', ...
               'EQW+RDM', 'EQW'};

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
% xticklabels(error_types);
ylabel('RMS Error');
title('RMS Error - Horizontal (Short Sequence)');

% Place the error type names vertically within the bars
for i = 1:numel(error_values)
    text(i, 0.2, error_types{i}, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vertical RMS Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [ rmse_ver_elvcn_rdm_50, rmse_ver_elvcn_50, ...
                 rmse_ver_elvcn_rdm_60, rmse_ver_elvcn_60, ...
                 rmse_ver_cn_rdm_li, rmse_ver_cn_li, ...
                 rmse_ver_cn_rdm_hv, rmse_ver_cn_hv, ...
                 rmse_ver_elv01_rdm, rmse_ver_elv01,...
                 rmse_ver_eqw_rdm, rmse_ver_eqw];
                
error_types = {'ELVCN+RDM-50', 'ELVCN-50', ...
               'ELVCN+RDM-60','ELVCN-60', ...
               'CN+RDM-light', 'CN-light', ...
               'CN+RDM-heavy', 'CN-heavy', ...
               'ELV+RDM', 'ELV', ...
               'EQW+RDM', 'EQW'};

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
% xticklabels(error_types);
ylabel('RMS Error');
title('RMS Error - Vertical (Short Sequence)');

for i = 1:numel(error_values)
    text(i, 0.2, error_types{i}, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
hold off;