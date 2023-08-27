eqw = load('03_data/EQW.mat');
rdm = load('03_data/RDM.mat');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equal weight %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_eqw, mean_error_ver_eqw, rmse_hor_eqw, ...
 rmse_ver_eqw, max_error_hor_eqw, max_error_ver_eqw]= ...
            calculate_error(eqw.user_pos(:, 1), eqw.user_pos(:, 2), eqw.user_pos(:, 3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Redundancy Matrix %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mean_error_hor_rdm, mean_error_ver_rdm, rmse_hor_rdm, ...
 rmse_ver_rdm, max_error_hor_rdm, max_error_ver_rdm]= ...
            calculate_error(rdm.user_pos(:, 1), rdm.user_pos(:, 2), rdm.user_pos(:, 3));