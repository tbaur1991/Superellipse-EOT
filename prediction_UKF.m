function [X,P] = prediction_UKF(X,P,dt,wm,wc,lambda,sigma_v,sigma_omega)
    % calculate sigma points
    nx = length(X);
    A = sqrt(nx + lambda) * chol(P)';
    sig_points = [zeros(size(X)) -A A];
    sig_points = sig_points + repmat(X, 1, size(sig_points, 2));
    
    % predict samples according to transition model
    for i = 1:2*nx+1
        sig_points(:,i) = motion_model(sig_points(:,i),dt);
    end

    % covariance matrix of process noise
    G = zeros(nx,2); G(1,1) = dt^2/2*cos(X(4)); G(2,1) = dt^2/2*sin(X(4)); 
    G(3,1) = dt; G(4,2) = dt^2/2; G(5,2) = dt;
    Q = G*diag([sigma_v sigma_omega].^2)*G';
    
    % calculate predicted mean and covariance
    X = sum(wm.*sig_points,2);
    P = (wc.*(sig_points - X))*(sig_points - X)';
    P = P + Q;
end