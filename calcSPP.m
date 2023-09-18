function [rec_pos, PDOP] = calcSPP( sats, obs, time_, snr, eph, v_light, el_mask)
% ----------------------------------------------------------------------- %
% function for calculation of single point positioning for a single epoch
% ----------------------------------------------------------------------- %

%% ----- Single point positioning -------------------------------------- %%
% Initialize variables
iteration = 1;
deltax = 9999;
m = length(sats);
delta_rho = zeros(m,1);
rho_init = zeros(m,1);
designmatrix_A = zeros(m,4);
Q_ll = eye(m);
% weight matrix

rec_0 = zeros(4,1);

%% 3.3. Weighting schemes based on both satellite elevation and C/N0
s0 = 10;
s1_threshold = 60;
A_snr = 30;
B_snr = 30;
a_cn = 0.01;
b_cn = 25;
while max(abs(deltax)) > 1e-3
    % iteration for each observation
    for i = 1:m

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TODO: write the calculation of designmatrix A, cofactor
        %%       matrix Q_ll and delta_l

        travel_time = obs(i) / v_light;

        if iteration > 1
            travel_time = rho_init(i) / v_light;
        end
        sow = time_ - travel_time;

        [sat_pos, clock_corr] = satellitePosition(sow, eph, sats(i));
        
        sat_pos_rot = e_r_corr(travel_time, sat_pos);
        
        [~,el(i),dist] = topocent(rec_0(1:3),sat_pos_rot - rec_0(1:3));
        
        [troposphere] = tropo(sin(el(i)),0,1013,293,50,0,0,0);
        sat_clock_off = clock_corr * v_light;
        A = dist;
        B = rec_0(4);
        C = sat_clock_off;
        D = troposphere;
        rho_init(i) = (A+B-C+D);
        delta_rho(i) = obs(i) - rho_init(i);

        xs_xu = (rec_0(1) - sat_pos_rot(1));
        ys_yu = (rec_0(2) - sat_pos_rot(2));
        zs_zu = (rec_0(3) - sat_pos_rot(3));

        designmatrix_A(i,:) = [xs_xu/dist, ys_yu/dist, zs_zu/dist, 1];

        %% 3.3. Weighting schemes based on both satellite elevation and C/N0
        %% PAPER: Using local redundancy to improve GNSS absolute positioning in harsh scenario.
        
        %% ELV + CN
        % if snr(i) < s1_threshold
        %     term1 = 10^-((snr(i) - s1_threshold)/B_snr);
        %     term2 = (A_snr / 10^-((s0-s1_threshold)/B_snr))-1;
        %     term3 = ((snr(i) - s1_threshold)/(s0-s1_threshold));
        %     r = term1 * (term2 * term3 + 1);
        % 
        %     Q_ll(i, i) = r / sin(el(i))^2;
        % else
        %     Q_ll(i, i) = 1;
        % end

        %% CN
        % Q_ll(i, i) = a_cn + (b_cn * (10^-(snr(i)/10)));

        %% ELV
        % Q_ll(i, i) = 1 / sin(el(i))^2;

        %% EQW
        Q_ll60(i, i) = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TODO: - calculate delta_x in a least square adjustment
    %%       - update position for next iteration
    %%
    % redundancy matrix
    I = eye(m);
    R = I - designmatrix_A * inv(designmatrix_A' * designmatrix_A) * designmatrix_A';
    R_diag = diag(R);
    R_diag(R_diag == 0) = 1;
    
    % wi = R_diag./diag(Q_ll);
    wi = 1./ diag(Q_ll);
    P = diag(wi);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Applying elevation Mask %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     idxDel = el(i) < 0;
%     if el(i) < 0
%         idxDel = i;
%         designmatrix_A(idxDel,:) = [];
%         P(idxDel,:) = [];
%         P(:,idxDel) = [];
%         delta_rho(idxDel) = [];
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Setup normal equations
    N = designmatrix_A' * P * designmatrix_A;
    n = designmatrix_A' * P * delta_rho;



    %approximation
    delta_rec_approx = inv(N) * n;

    rec_approx = rec_0 + delta_rec_approx;
    delta_rho_approx = designmatrix_A * delta_rec_approx;
%     rho_approx = rho_init + delta_rho_approx;

    %residual
%     residual = delta_rho_approx - delta_rho;

    deltax = abs(rec_approx - rec_0);
    rec_0 = rec_approx;
    rec_pos = rec_approx(1:3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % control computation
%     control = designmatrix_A' * P * residual;

    iteration = iteration + 1;
    % disp(iteration);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO: calculate PDOP
pdop = diag(inv(N));
PDOP = sqrt(sum(pdop(1:3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% https://github.com/gnss-sdr/gnss-sdr/blob/main/src/utils/matlab/libs/geoFunctions/tropo.m

