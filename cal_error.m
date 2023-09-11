clc;
clear all;
close all;
% working directory
path = pwd;
% add functions
addpath([path,'/04_functions']);


eqw = load('03_data/EQW_1500_2500.mat');
eqw_rdm = load('03_data/EQW+RDM_1500_2500.mat');
elv01 = load('03_data/ELV01_1500_2500.mat');
elv01_rdm = load('03_data/ELV01+RDM_1500_2500.mat');
cn = load('03_data/CN_1500_2500.mat');
cn_rdm = load('03_data/CN+RDM_1500_2500.mat');
elvcn_45 = load('03_data/ELVCN_1500_2500_thresh_45.mat');
elvcn_50 = load('03_data/ELVCN_1500_2500_thresh_50.mat');
elvcn_60 = load('03_data/ELVCN_1500_2500_thresh_60.mat');
elvcn_rdm_45 = load('03_data/ELVCN+RDM_1500_2500_thresh_45.mat');
elvcn_rdm_50 = load('03_data/ELVCN+RDM_1500_2500_thresh_50.mat');
elvcn_rdm_60 = load('03_data/ELVCN+RDM_1500_2500_thresh_60.mat');

start = 1500;
stop = 2500; 

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
%% CN %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_cn, mean_error_ver_cn, rmse_hor_cn, ...
 rmse_ver_cn, max_error_hor_cn, max_error_ver_cn]= ...
            calculate_error(cn.user_pos(start:stop, 1), cn.user_pos(start:stop, 2), cn.user_pos(start:stop, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CN + RDM %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_cn_rdm, mean_error_ver_cn_rdm, rmse_hor_cn_rdm, ...
 rmse_ver_cn_rdm, max_error_hor_cn_rdm, max_error_ver_cn_rdm]= ...
            calculate_error(cn_rdm.user_pos(start:stop, 1), cn_rdm.user_pos(start:stop, 2), cn_rdm.user_pos(start:stop, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV + CN 45 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elvcn_45, mean_error_ver_elvcn_45, rmse_hor_elvcn_45, ...
 rmse_ver_elvcn_45, max_error_hor_elvcn_45, max_error_ver_elvcn_45]= ...
            calculate_error(elvcn_45.user_pos(start:stop, 1), elvcn_45.user_pos(start:stop, 2), elvcn_45.user_pos(start:stop, 3));

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
%% ELV + CN + RDM 45%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elvcn_rdm_45, mean_error_ver_elvcn_rdm_45, rmse_hor_elvcn_rdm_45, ...
 rmse_ver_elvcn_rdm_45, max_error_hor_elvcn_rdm_45, max_error_ver_elvcn_rdm_45]= ...
            calculate_error(elvcn_rdm_45.user_pos(start:stop, 1), elvcn_rdm_45.user_pos(start:stop, 2), elvcn_rdm_45.user_pos(start:stop, 3));

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

error_values = [max_error_hor_elvcn_rdm_45, max_error_hor_elvcn_rdm_50, max_error_hor_elvcn_rdm_60, ...
                max_error_hor_elvcn_45, max_error_hor_elvcn_50, max_error_hor_elvcn_60, ...
                max_error_hor_cn_rdm, max_error_hor_cn, max_error_hor_elv01_rdm, ...
                max_error_hor_elv01, max_error_hor_eqw_rdm, max_error_hor_eqw];
                
error_types = {'ELVCN+RDM_45', 'ELVCN+RDM_50', 'ELVCN+RDM_60', ...
               'ELVCN_45', 'ELVCN_50', 'ELVCN_60', 'CN+RDM', ...
               'CN', 'ELV+RDM', 'ELV', 'EQW+RDM', 'EQW'};
         

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
xticklabels(error_types);
ylabel('Max Error');
title('Max Error - Horizontal (High Visibility Sequence)');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vertical Max Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [max_error_ver_elvcn_rdm_45, max_error_ver_elvcn_rdm_50, max_error_ver_elvcn_rdm_60, ...
                max_error_ver_elvcn_45, max_error_ver_elvcn_50, max_error_ver_elvcn_60, ...
                max_error_ver_cn_rdm, max_error_ver_cn, max_error_ver_elv01_rdm, ...
                max_error_ver_elv01, max_error_ver_eqw_rdm, max_error_ver_eqw];
                
error_types = {'ELVCN+RDM_45', 'ELVCN+RDM_50', 'ELVCN+RDM_60', ...
               'ELVCN_45', 'ELVCN_50', 'ELVCN_60', 'CN+RDM', ...
               'CN', 'ELV+RDM', 'ELV', 'EQW+RDM', 'EQW'};

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
xticklabels(error_types);
ylabel('Max Error');
title('Max Error - Vertical (High Visibility Sequence)');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Horizontal Mean Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [mean_error_hor_elvcn_rdm_45, mean_error_hor_elvcn_rdm_50, mean_error_hor_elvcn_rdm_60, ...
                mean_error_hor_elvcn_45, mean_error_hor_elvcn_50, mean_error_hor_elvcn_60, ...
                mean_error_hor_cn_rdm, mean_error_hor_cn, mean_error_hor_elv01_rdm, ...
                mean_error_hor_elv01, mean_error_hor_eqw_rdm, mean_error_hor_eqw];
                
error_types = {'ELVCN+RDM_45', 'ELVCN+RDM_50', 'ELVCN+RDM_60', ...
               'ELVCN_45', 'ELVCN_50', 'ELVCN_60', 'CN+RDM', ...
               'CN', 'ELV+RDM', 'ELV', 'EQW+RDM', 'EQW'};

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
xticklabels(error_types);
ylabel('Mean Error');
title('Mean Error - Horizontal (High Visibility Sequence)');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vertical Mean Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [mean_error_ver_elvcn_rdm_45, mean_error_ver_elvcn_rdm_50, mean_error_ver_elvcn_rdm_60, ...
                mean_error_ver_elvcn_45, mean_error_ver_elvcn_50, mean_error_ver_elvcn_60, ...
                mean_error_ver_cn_rdm, mean_error_ver_cn, mean_error_ver_elv01_rdm, ...
                mean_error_ver_elv01, mean_error_ver_eqw_rdm, mean_error_ver_eqw];
                
error_types = {'ELVCN+RDM_45', 'ELVCN+RDM_50', 'ELVCN+RDM_60', ...
               'ELVCN_45', 'ELVCN_50', 'ELVCN_60', 'CN+RDM', ...
               'CN', 'ELV+RDM', 'ELV', 'EQW+RDM', 'EQW'};

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
xticklabels(error_types);
ylabel('Mean Error');
title('Mean Error - Vertical (High Visibility Sequence)');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Horizontal RMS Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [rmse_hor_elvcn_rdm_45, rmse_hor_elvcn_rdm_50, rmse_hor_elvcn_rdm_60, ...
                rmse_hor_elvcn_45, rmse_hor_elvcn_50, rmse_hor_elvcn_60, ...
                rmse_hor_cn_rdm, rmse_hor_cn, rmse_hor_elv01_rdm, ...
                rmse_hor_elv01, rmse_hor_eqw_rdm, rmse_hor_eqw];
                
error_types = {'ELVCN+RDM_45', 'ELVCN+RDM_50', 'ELVCN+RDM_60', ...
               'ELVCN_45', 'ELVCN_50', 'ELVCN_60', 'CN+RDM', ...
               'CN', 'ELV+RDM', 'ELV', 'EQW+RDM', 'EQW'};

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
xticklabels(error_types);
ylabel('RMS Error');
title('RMS Error - Horizontal (High Visibility Sequence)');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vertical RMS Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;

error_values = [rmse_ver_elvcn_rdm_45, rmse_ver_elvcn_rdm_50, rmse_ver_elvcn_rdm_60, ...
                rmse_ver_elvcn_45, rmse_ver_elvcn_50, rmse_ver_elvcn_60, ...
                rmse_ver_cn_rdm, rmse_ver_cn, rmse_ver_elv01_rdm, ...
                rmse_ver_elv01, rmse_ver_eqw_rdm, rmse_ver_eqw];
                
error_types = {'ELVCN+RDM_45', 'ELVCN+RDM_50', 'ELVCN+RDM_60', ...
               'ELVCN_45', 'ELVCN_50', 'ELVCN_60', 'CN+RDM', ...
               'CN', 'ELV+RDM', 'ELV', 'EQW+RDM', 'EQW'};

bar(1:numel(error_values), error_values);
xticks(1:numel(error_values));
xticklabels(error_types);
ylabel('RMS Error');
title('RMS Error - Vertical (High Visibility Sequence)');
hold off;