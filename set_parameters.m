% set measurement uncertainty
sig_m = 0.1;

% System state allocation for superellipse shape estimation
if strcmp(dim,'2D')
    if taper == true
        nx = 9; % [x,y,v,or,omega,a,b,epsilon,ty]
    elseif taper == false
        nx = 8; % [x,y,v,or,omega,a,b,epsilon]
    else
        error('Wrong definition of variable taper. Choose either true or false.')
    end
    % allocation for UKF system state
    X_UKF = zeros(nx,nSamples);
    P_UKF = zeros(nx,nx,nSamples);
    % allocation for S2KF system state
    X_S2KF = zeros(nx,nSamples);
    P_S2KF = zeros(nx,nx,nSamples);
elseif strcmp(dim,'3D')
    if taper == true
        nx = 11; % [x,y,z,v,or,omega,a,b,h,epsilon,ty]
    elseif taper == false
        nx = 10; % [x,y,z,v,or,omega,a,b,h,epsilon]
    else
        error('Wrong definition of variable taper. Choose either true or false.')
    end
    % allocation for single equation system state
    X_3D = zeros(nx,nSamples);
    P_3D = zeros(nx,nx,nSamples);
    % allocations for decoupled system state
    X_2D = zeros(nx-2,nSamples);
    P_2D = zeros(nx-2,nx-2,nSamples);
    X_3DE = zeros(2,nSamples);
    P_3DE = zeros(2,2,nSamples);
else
    error('Wrong definition of variable taper. Choose either true or false.')
end

% UKF sigma point parameters
alpha = 1e-3;
beta = 2;
kappa = 0;
lambda = alpha^2*(nx + kappa) - nx;

% calculate weights for sigma points
wm(1) = lambda/(nx + lambda);
wc(1) = lambda/(nx + lambda) + (1 - alpha^2 + beta);
wm(2:2*nx + 1) = 1/(2*(nx + lambda));
wc(2:2*nx + 1) = 1/(2*(nx + lambda));

% calculate additional samples for LRKF filters
if strcmp(dim,'2D')
    % get samples of S2KF
    S.name = 'Symmetric LCD';
    S.sampling = GaussianSamplingLCD();
    S.sampling.setNumSamplesByFactor(10);
    [samples,~,numSamples] = S.sampling.getStdNormalSamples(nx);
elseif strcmp(dim,'3D')
    % weights for the 2D superellipse filter
    lambda_2D = alpha^2*(nx-2 + kappa) - (nx-2);
    wm_2D(1) = lambda_2D/(nx-2 + lambda_2D);
    wc_2D(1) = lambda_2D/(nx-2 + lambda_2D) + (1 - alpha^2 + beta);
    wm_2D(2:2*(nx-2) + 1) = 1/(2*(nx-2 + lambda_2D));
    wc_2D(2:2*(nx-2) + 1) = 1/(2*(nx-2 + lambda_2D));
    % weights for the 3D extent state filter
    lambda_3DE = alpha^2*(2+n_upd + kappa) - (2+n_upd);
    wm_3DE(1) = lambda_3DE/(2+n_upd + lambda_3DE);
    wc_3DE(1) = lambda_3DE/(2+n_upd + lambda_3DE) + (1 - alpha^2 + beta);
    wm_3DE(2:2*(2+n_upd) + 1) = 1/(2*(2+n_upd + lambda_3DE));
    wc_3DE(2:2*(2+n_upd) + 1) = 1/(2*(2+n_upd + lambda_3DE));
    % weights for the single 3D superellipse filter update
    lambda_3D = alpha^2*(nx+n_upd + kappa) - (nx+n_upd);
    wm_3D(1) = lambda_3D/(nx+n_upd + lambda_3D);
    wc_3D(1) = lambda_3D/(nx+n_upd + lambda_3D) + (1 - alpha^2 + beta);
    wm_3D(2:2*(nx+n_upd) + 1) = 1/(2*(nx+n_upd + lambda_3D));
    wc_3D(2:2*(nx+n_upd) + 1) = 1/(2*(nx+n_upd + lambda_3D));
end

% initialize inside measurement variance estimation
varInX1 = 0; varInY1 = 0; meanInX1 = 0; meanInY1 = 0;
varInX2 = 0; varInY2 = 0; meanInX2 = 0; meanInY2 = 0;
tao = 200;

% evaluation memory allocation
iou1 = zeros(1,nSamples);
iou2 = zeros(1,nSamples);
t1 = zeros(1,nSamples);
t2 = zeros(1,nSamples);

% set number of basis points for polygonal chain approximation
nPoly = 25;