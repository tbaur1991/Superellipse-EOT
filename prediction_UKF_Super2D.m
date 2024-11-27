function [X,P] = prediction_UKF_Super2D(X,P,wm,wc,lambda,dt,sigma_v,sigma_omega,taper)
    % get dymension of system state
    nx = length(X);

    % Calculate Sigma Points
    A = sqrt(nx + lambda) * chol(P)';
    sig_points = [zeros(size(X)) -A A];
    sig_points = sig_points + repmat(X, 1, size(sig_points, 2));
    
    % predict samples according to transition model
    for i = 1:2*nx+1
        sig_points(:,i) = motion_model2D(sig_points(:,i),dt);
    end

    % covariance matrix of process noise
    sigma_a = 1e-3; sigma_b = 1e-3; sigma_e = 1e-4; sigma_t = 1e-4;
    if taper
        G = zeros(nx,6); G(1,1) = dt^2/2*cos(X(4)); G(2,1) = dt^2/2*sin(X(4)); 
        G(3,1) = dt; G(4,2) = dt^2/2; G(5,2) = dt; G(6,3) = dt; G(7,4) = dt; G(8,5) = dt; G(9,6) = dt;
        Q = G*diag([sigma_v sigma_omega sigma_a sigma_b sigma_e sigma_t].^2)*G';
    else
        G = zeros(nx,5); G(1,1) = dt^2/2*cos(X(4)); G(2,1) = dt^2/2*sin(X(4)); 
        G(3,1) = dt; G(4,2) = dt^2/2; G(5,2) = dt; G(6,3) = dt; G(7,4) = dt; G(8,5) = dt;
        Q = G*diag([sigma_v sigma_omega sigma_a sigma_b sigma_e].^2)*G';
    end
    
    % calculate predicted mean
    X = sum(wm.*sig_points,2);
    P = (wc.*(sig_points - X))*(sig_points - X)';
    P = P + Q;
end

% end of function code
% *************************************************************************
%
%
% *************************************************************************
% end of document "prediction_S2KF.m"
% *************************************************************************
