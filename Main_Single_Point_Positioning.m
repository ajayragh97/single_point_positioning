clc
clearvars
close all
format longg

%=========================================================================%
%=====                                                               =====%
%=====          Single Point Positioning with GPS                    =====%
%=====                                                               =====%
%=====        MSR-05 Methods for quality assurance                   =====%
%=====                                                               =====%
%=========================================================================%

%% ----- Definition of constant and directory--------------------------- %%
% speed of light (global variable)
global v_light; 
v_light = 299792458;  % [m/s]
el_mask = 0;

% working directory
path = pwd;
% add functions
addpath([path,'/02_functions']);

%% ----- Import data --------------------------------------------------- %%
% navigation file
nav = '01_data/master_igg.16n';
% observation RINEX file
rnx = '01_data/master_igg_short.16o';

% load *.mat files if they already exist 
if isfile('01_data/eph.mat') && isfile('01_data/observations.mat') && ...
        isfile('01_data/satellites.mat') && isfile('01_data/time.mat') && ...
        isfile('01_data/snr.mat')
    
    load('01_data/eph.mat');
    load('01_data/observations.mat');
    load('01_data/time.mat');
    load('01_data/satellites.mat');
    load('01_data/snr.mat');
else
    % ephemeris
    [~, eph, ephTime, ~] = readEphGps(nav,'prepareData',40);
    % observations
    [errNum, header, obs, date, time] = readRinex(rnx,'prepareData',40,1);
    % reallocate dataset
    n = size(date,1);
    observations = cell(n,1);
    satellites = cell(n,1);
    snr = cell(n,1);
    for i = 1:n
        observations{i} = obs(i).obsG(obs(i).satG,1);
        satellites{i} = obs(i).satG;
        snr{i} = obs(i).obsG(obs(i).satG,3);
    end
    % save mat files
    save('01_data/eph_large.mat','eph')
    save('01_data/time_large.mat','time')
    save('01_data/satellites_large.mat','satellites')
    save('01_data/observations_large.mat','observations')
    save('01_data/snr_large.mat','snr')
end

%% ----- Single Point Positioning -------------------------------------- %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO: initialize vector/matrix for parameter and results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----- positioning
user_pos = zeros(length(time), 3);
PDOP = zeros(length(time), 1);
for i=1:length(time)
% for i=1500:2500
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TODO: write/complete the function 'calcSPP' to calculate the receiver
    %%       position based on measurements of a single epoch at time(i)
                                                                                                                                                                                                                                
    sats = satellites{i};
    obs = observations{i};
    time_ = time(i);
    snr_ = snr{i};

    [rec_pos, pdop] = calcSPP(sats, obs, time_, snr_, eph, v_light, el_mask);

     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %% TODO: save output of 'calcSPP' in the initialized result vector/matrix
    user_pos(i,:) = rec_pos;
    PDOP(i) = pdop;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% Plot each coordinate as a line plot with respect to the time values


%% ----- Plot results -------------------------------------------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% TODO: plot results in cartesian coordinates with corresponding time
% % Calculate the mean of each coordinate dimension
% mean_pos = mean(user_pos);
% % Subtract the mean from each coordinate dimension
% centered_pos = user_pos - mean_pos;
% % Plot the centered coordinates as a time series
% figure;
% 
% xPos = 0.05;
% yPos = 0.9;
% annotation('textbox', [xPos, yPos, 0.1, 0.1], 'String', ['Elevation Mask: ' num2str(el_mask)], 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');
% 
% plot(time, centered_pos(:, 1), 'r', 'LineWidth', 1.5);
% hold on;
% plot(time, centered_pos(:, 2), 'g', 'LineWidth', 1.5);
% plot(time, centered_pos(:, 3), 'b', 'LineWidth', 1.5);
% hold off;
% xlabel('Observation Epoch');
% ylabel('Coordinate Value (Mean-Centered)');
% title('Mean-Centered Coordinates over Time');
% legend('X', 'Y', 'Z');
% grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO: plot PDOP values and satellites with corresponding time
figure;
hold on;
plot(PDOP, 'b');
yyaxis right;
% numSatellites = zeros(size(satellites));
% numSatellites(1500:2500) = cellfun(@numel, satellites(1500:2500));
numSatellites= cellfun(@numel, satellites);
% plot(numSatellites(1500:2500), 'r');
plot(numSatellites, 'r');
ylabel('PDOP');
title('PDOP and Number of Satellites');

xlabel('Time Step');
yyaxis left;
ylabel('satellites');
legend( 'PDOP','Number of Satellites');
% title('Number of Satellites');

xPos = 0.05;
yPos = 0.9;
annotation('textbox', [xPos, yPos, 0.1, 0.1], 'String', ['Elevation Mask: ' num2str(el_mask)], 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra plots %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% user_pos_lla = ecef2lla(user_pos);
% geoscatter(user_pos_lla(:,1),user_pos_lla(:,2),'^','filled');
% latlim = [50+43/60+35/3600 50+43/60+45/3600]; % Latitude limits for Bonn
% lonlim = [7+4/60+50/3600 7+5/60+30/3600];    % Longitude limits for Bonn
% geolimits(latlim, lonlim);
% title('Calculated SPP positions');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save('03_data/EQW_5300_6000.mat','user_pos')
% save('03_data/EQW+RDM_5300_6000.mat','user_pos')
% save('03_data/ELV01_5300_6000.mat','user_pos')
% save('03_data/ELV01+RDM_5300_6000.mat','user_pos')
% save('03_data/CN_large.mat','user_pos')
% save('03_data/CN+RDM_large.mat','user_pos')
% save('03_data/ELVCN_short_thresh_60.mat','user_pos')
% save('03_data/ELVCN+RDM_short_thresh_45.mat','user_pos')




