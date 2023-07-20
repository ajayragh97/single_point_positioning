function [rec_pos, PDOP] = calcSPP( sats, obs, time_, eph, v_light, el_mask)
% ----------------------------------------------------------------------- %
% function for calculation of single point positioning for a single epoch
% ----------------------------------------------------------------------- %

%% ----- Single point positioning -------------------------------------- %%
% Initialize variables
iteration = 1;
deltax = 9999;
m = length(sats);
delta_rho = zeros(m,1);
designmatrix_A = zeros(m,4);
Q_ll = eye(m);
rec_0 = zeros(4,1);
while max(abs(deltax)) > 1e-3
    % iteration for each observation
    for i = 1:m
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TODO: write the calculation of designmatrix A, cofactor
        %%       matrix Q_ll and delta_l
        sv = sats(i);
        travel_dist = obs(i);
        obs_time = time_;
        sow = obs_time - travel_dist/v_light;
        if iteration > 1
            sow = obs_time - rho_init(i)/v_light;
        end
        [sat_pos, clock_corr] = satellitePosition(sow, eph, sv);
        travel_time = travel_dist / v_light;
        sat_pos_rot = e_r_corr(travel_time, sat_pos);
        [az,el(i,1),dist] = topocent(rec_0(1:3),sat_pos_rot - rec_0(1:3));

        [troposphere] = tropo(sin(el(i)),0,1013,293,50,0,0,0);
        sat_clock_off = clock_corr * v_light;
        A = norm(sat_pos_rot - rec_0(1:3));
        B = rec_0(4);
        C = sat_clock_off;
        D = troposphere;
        rho_init(i,1) = (A+B-C+D) ;
        delta_rho(i,:) = obs(i) - rho_init(i);
        sat = sat_pos_rot;
        xs_xu = (2*rec_0(1) - 2*sat(1));
        ys_yu = (2*rec_0(2) - 2*sat(2));
        zs_zu = (2*rec_0(3) - 2*sat(3));
        denom = (2*((rec_0(1) - sat(1))^2 + (rec_0(2) - sat(2))^2 + (rec_0(3) - sat(3))^2)^(1/2));
        designmatrix_A(i,:) = [xs_xu/denom, ys_yu/denom, zs_zu/denom, 1];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TODO: - calculate delta_x in a least square adjustment
    %%       - update position for next iteration
    %%
    % redundancy matrix
    I = eye(m);
    R = I - designmatrix_A * inv(designmatrix_A' * designmatrix_A) * designmatrix_A';

    % weight matrix
    P = inv(Q_ll);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Applying elevation Mask %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idxDel = el<0;
    designmatrix_A(idxDel,:) = [];
    P(idxDel,:) = [];
    P(:,idxDel) = [];
    delta_rho(idxDel) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Setup normal equations
    N = designmatrix_A' * P * designmatrix_A;
    n = designmatrix_A' * P * delta_rho;



    %approximation
    delta_rec_approx = inv(N) * n;

    rec_approx = rec_0 + delta_rec_approx;
    delta_rho_approx = designmatrix_A * delta_rec_approx;
    %rho_approx = rho_init + delta_rho_approx;

    %residual
    residual = delta_rho_approx - delta_rho;

    deltax = abs(rec_approx - rec_0);
    rec_0 = rec_approx;
    rec_pos = rec_approx(1:3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % control computation
    control = designmatrix_A' * P * residual;

    iteration = iteration + 1;
    % disp(iteration);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO: calculate PDOP
pdop = diag(inv(N));
PDOP = sqrt(sum(pdop(1:3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

