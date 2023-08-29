% working directory
path = pwd;
% add functions
addpath([path,'/04_functions']);


eqw = load('03_data/EQW.mat');
eqw_rdm = load('03_data/EQW+RDM.mat');
elv01 = load('03_data/ELV01.mat');
elv01_rdm = load('03_data/ELV01+RDM.mat');
cn = load('03_data/CN.mat');
cn_rdm = load('03_data/CN+RDM.mat');
elvcn = load('03_data/ELVCN.mat');
elvcn_rdm = load('03_data/ELVCN+RDM.mat');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EQW %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_eqw, mean_error_ver_eqw, rmse_hor_eqw, ...
 rmse_ver_eqw, max_error_hor_eqw, max_error_ver_eqw]= ...
            calculate_error(eqw.user_pos(:, 1), eqw.user_pos(:, 2), eqw.user_pos(:, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EQW + RDM %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_eqw_rdm, mean_error_ver_eqw_rdm, rmse_hor_eqw_rdm, ...
 rmse_ver_eqw_rdm, max_error_hor_eqw_rdm, max_error_ver_eqw_rdm]= ...
            calculate_error(eqw_rdm.user_pos(:, 1), eqw_rdm.user_pos(:, 2), eqw_rdm.user_pos(:, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elv01, mean_error_ver_elv01, rmse_hor_elv01, ...
 rmse_ver_elv01, max_error_hor_elv01, max_error_ver_elv01]= ...
            calculate_error(elv01.user_pos(:, 1), elv01.user_pos(:, 2), elv01.user_pos(:, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV + RDM %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elv01_rdm, mean_error_ver_elv01_rdm, rmse_hor_elv01_rdm, ...
 rmse_ver_elv01_rdm, max_error_hor_elv01_rdm, max_error_ver_elv01_rdm]= ...
            calculate_error(elv01_rdm.user_pos(:, 1), elv01_rdm.user_pos(:, 2), elv01_rdm.user_pos(:, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CN %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_cn, mean_error_ver_cn, rmse_hor_cn, ...
 rmse_ver_cn, max_error_hor_cn, max_error_ver_cn]= ...
            calculate_error(cn.user_pos(:, 1), cn.user_pos(:, 2), cn.user_pos(:, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CN + RDM %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_cn_rdm, mean_error_ver_cn_rdm, rmse_hor_cn_rdm, ...
 rmse_ver_cn_rdm, max_error_hor_cn_rdm, max_error_ver_cn_rdm]= ...
            calculate_error(cn_rdm.user_pos(:, 1), cn_rdm.user_pos(:, 2), cn_rdm.user_pos(:, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV + CN %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elvcn, mean_error_ver_elvcn, rmse_hor_elvcn, ...
 rmse_ver_elvcn, max_error_hor_elvcn, max_error_ver_elvcn]= ...
            calculate_error(elvcn.user_pos(:, 1), elvcn.user_pos(:, 2), elvcn.user_pos(:, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELV + CN + RDM %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_elvcn_rdm, mean_error_ver_elvcn_rdm, rmse_hor_elvcn_rdm, ...
 rmse_ver_elvcn_rdm, max_error_hor_elvcn_rdm, max_error_ver_elvcn_rdm]= ...
            calculate_error(elvcn_rdm.user_pos(:, 1), elvcn_rdm.user_pos(:, 2), elvcn_rdm.user_pos(:, 3));


