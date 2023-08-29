function [mean_error_horizontal,mean_error_vertical, ...
    rmse_horizontal,rmse_vertical, ...
    max_error_horizontal, max_error_vertical] = calculate_error(x,y,z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UTM %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lat_s, lon_s, height_s] = convertECEFtoGRS80(x, y, z);
[e_s, n_s ,h_s] = ell2utm32(lat_s, lon_s, height_s);

mean_error_e = mae(e_s - mean(e_s));
mean_error_n = mae(n_s - mean(n_s));
mean_error_h = mae(h_s - mean(h_s));
mean_error_horizontal = norm([mean_error_e mean_error_n]);
mean_error_vertical = mean_error_h;

rmse_e = rmse(e_s, mean(e_s));
rmse_n = rmse(n_s, mean(n_s));
rmse_h = rmse(h_s, mean(h_s));
rmse_horizontal = norm([rmse_e rmse_n]);
rmse_vertical = rmse_h;

max_error_e = max(abs(e_s - mean(e_s)));
max_error_n = max(abs(n_s - mean(n_s)));
max_error_h = max(abs(h_s - mean(h_s)));
max_error_horizontal = norm([max_error_e max_error_n]);
max_error_vertical = max_error_h;
end