function [X,P] = prediction_UKF_Super3D(X,P,dt,sigma_v,sigma_omega,wm,wc,lambda,taper)
    % sigma point parameters
    nx = length(X);

    %Calculate Sigma Points
    A = sqrt(nx + lambda) * chol(P)';
    sig_points = [zeros(size(X)) -A A];
    sig_points = sig_points + repmat(X, 1, size(sig_points, 2));
    
    % predict samples according to transition model
    for i = 1:2*nx+1
        sig_points(:,i) = motion_model3D(sig_points(:,i),dt);
    end

    % covariance matrix of process noise
    sigma_a = 1e-3; sigma_b = 1e-3; sigma_h = 1e-3; sigma_e = 1e-4;
    if taper
        sigma_t = 1e-5;
        G = zeros(nx,7); 
    else
        G = zeros(nx,6);
    end
    G(1,1) = dt^2/2*cos(X(5)); G(2,1) = dt^2/2*sin(X(5)); G(3,1) = dt^2/2; 
    G(4,1) = dt; G(5,2) = dt^2/2; G(6,2) = dt; G(7,3) = dt; G(8,4) = dt; G(9,5) = dt; G(10,6) = dt; 
    if taper
        G(11,7) = dt;
        Q = G*diag([sigma_v sigma_omega sigma_a sigma_b sigma_h sigma_e sigma_t].^2)*G';
    else
        Q = G*diag([sigma_v sigma_omega sigma_a sigma_b sigma_h sigma_e].^2)*G';
    end
    
    % calculate predicted mean
    X = sum(wm.*sig_points,2);
    P = (wc.*(sig_points - X))*(sig_points - X)';
    P = P + Q;
end
