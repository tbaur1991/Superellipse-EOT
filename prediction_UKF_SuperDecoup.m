% *************************************************************************
%    I   SSSS   DDDD
%    I   S      D   D 
%    I   SSSS   D   D    Institut fr Systemdynamik
%    I      S   D   D
%    I   SSSS   DDDD   
% *************************************************************************
% Author: Tim Baur
% University: University of Applied Sciences "HTWG Konstanz"
% E-Mail: tbaur@htwg-konstanz.de
% Creation date: 11.01.2023
% File name: prediction_S2KF.m
% *************************************************************************
% Content: In this file prediction step of S2KF based on Coordinated Turn 
% Model with Polar Velocity is performed.
%
% *************************************************************************
%% start of function code
function [X_2D,P_2D,X_3DE,P_3DE] = prediction_UKF_SuperDecoup(X_2D,P_2D,X_3DE,P_3DE,dt,sigma_v,sigma_omega,wm_2D,wc_2D,lambda_2D,taper)
    % prediction for 2D state
    % sigma point parameters
    nx_2D = length(X_2D);

    % calculate Sigma Points
    A = sqrt(nx_2D + lambda_2D) * chol(P_2D)';
    sig_points = [zeros(size(X_2D)) -A A];
    sig_points = sig_points + repmat(X_2D, 1, size(sig_points, 2));
    
    % predict samples according to transition model
    for i = 1:2*nx_2D+1
        sig_points(:,i) = motion_model2D(sig_points(:,i),dt);
    end

    % covariance matrix of process noise
    sigma_a = 1e-3; sigma_b = 1e-3; sigma_e = 1e-4; 
    if taper
        sigma_t = 1e-5;
        G = zeros(nx_2D,6); 
    else
        G = zeros(nx_2D,5);
    end
    G(1,1) = dt^2/2*cos(X_2D(4)); G(2,1) = dt^2/2*sin(X_2D(4)); 
    G(3,1) = dt; G(4,2) = dt^2/2; G(5,2) = dt; G(6,3) = dt; G(7,4) = dt; G(8,5) = dt; 
    if taper
        G(9,6) = dt;
        Q = G*diag([sigma_v sigma_omega sigma_a sigma_b sigma_e sigma_t].^2)*G';
    else
        Q = G*diag([sigma_v sigma_omega sigma_a sigma_b sigma_e].^2)*G';
    end
    
    % calculate predicted mean
    X_2D = sum(wm_2D.*sig_points,2);
    P_2D = (wc_2D.*(sig_points - X_2D))*(sig_points - X_2D)';
    P_2D = P_2D + Q;

    % linear prediction for line state
    nx_3DE = length(X_3DE);
    sigma_z = 1e-3; sigma_h = 1e-3;
    G = zeros(nx_3DE); G(1,1) = dt; G(2,2) = dt;
    Q = G*diag([sigma_z sigma_h].^2)*G';
    F = eye(2);
    X_3DE = F*X_3DE;
    P_3DE = F*P_3DE*F' + Q;
end

% end of function code
% *************************************************************************
%
%
% *************************************************************************
% end of document "prediction_S2KF.m"
% *************************************************************************
