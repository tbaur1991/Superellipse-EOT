% set measurement uncertainty
sig_m = 0.1;

% System state allocation for superellipse shape estimation
nx = 9; % [x,y,v,or,omega,a,b,epsilon,ty]
X = zeros(nx,nSamples);
P = zeros(nx,nx,nSamples);
% measurements to be updated in a single loop
n_upd = 10;

% UKF sigma point parameters
alpha = 0.1;
beta = 3;
kappa = 0;
lambda = alpha^2*(nx + kappa) - nx;

% calculate weights for sigma points
wm(1) = lambda/(nx + lambda);
wc(1) = lambda/(nx + lambda) + (1 - alpha^2 + beta);
wm(2:2*nx + 1) = 1/(2*(nx + lambda));
wc(2:2*nx + 1) = 1/(2*(nx + lambda));

% initialize inside measurement variance estimation
varIn = 0; tao_varIn = 200;

% evaluation memory allocation
iou = zeros(1,nSamples);