% *************************************************************************
%    I   SSSS   DDDD
%    I   S      D   D 
%    I   SSSS   D   D    Institute of System Dynamics
%    I      S   D   D
%    I   SSSS   DDDD   
% *************************************************************************
% Author: Tim Baur
% University: University of Applied Sciences "HTWG Konstanz"
% E-Mail: tbaur@htwg-konstanz.de
% Creation date: 30.07.2024
% File name: main.m
% *************************************************************************
% Content: In this file, the simulation environment and estimation 
% procedure for a superellipse shape tracking in 2D or 3D using the research 
% vessel 'Solgenia' is handled. Measurements are generated from a simulated
% LiDAR sensor: https://de.mathworks.com/help/uav/ref/uavlidarpointcloudgenerator-system-object.html
% *************************************************************************

%% start of function code
close all; clear; clc

% load reference trajectory and simulated measurement sources
load('ref_trajectory.mat')
load('solgenia_sim_meas.mat')

%% settings to be chosen in this simulation

% choose artificieal measurement noise for processing a mixture of boundary
% and interior measurements to be used or not.
art_noise = true;

% choose if tapering coefficient is used for superellipse model
taper = true;

% choose if '2D' or '3D' models are applied
dim = '2D';

% choose if measurement source is approximated using 'radial' or 'projected'
% approximation for superellipse models
source = 'radial';

% set number of measurements to be processed in a single update step
n_upd = 20;

% choose if estimates and measurements are plotted
do_plots = true;

%% simulation environment and EOT filter

% do memory allocation and parameter settings for superellipse model and
% simulation environment
set_parameters

% start simulation
for k = 1:nSamples
    % display time step
    disp(['timestep: ',num2str(k)])

    % get new measurement set and add noise
    meas{k,1} = meas{k,1} + sig_m*randn(size(meas{k,1}));

    if k == 1
        % do state initialization
        pos = mean(meas{k,1}(1:2,:),2);
        % init orientation
        c = cov(meas{k,1}(1,:),meas{k,1}(2,:));
        [v,d] = eig(c); [~,i] = max(diag(d));
        or = atan2(v(2,i),v(1,i));
        % init shape
        R_rot = [cos(or) -sin(or); sin(or) cos(or)];
        mm = R_rot'*(meas{k,1}(1:2,:) - pos);
        a = max(abs(mm(1,:)));
        b = abs(max(mm(2,:)));
        % init superellipse shape tracking
        X(:,k) = [pos;0;or;0;a-1;b-1;3;0];
        P(:,:,k) = blkdiag(c,2,10*pi/180*eye(2),0.5*eye(2),1,1/3);
    else
        % do prediction step for superellipse shape tracking
        [X(:,k),P(:,:,k)] = prediction_UKF(X(:,k-1),P(:,:,k-1),dt,wm,wc,lambda,sigma_v,sigma_omega);
    end
    
    % generate random permutation of measurement set
    idx = randperm(size(meas{k,1},2));

    % do measurement update for superellipse shape tracking
    [X(:,k),P(:,:,k),varIn] = superellipse2D_UKF_update(X(:,k),P(:,:,k),sig_m,meas{k,1}(1:2,idx),wm,wc,lambda,n_upd,varIn,tao_varIn);

    % plot measurements, reference and estimate
    plot_object(meas{k,1}(1:2,:),X(:,k),convHullRef,X_Ref(:,k))

    % calculate intersection over union using reference convex hull
    iou(k) = iouSuper(X(:,k),X_Ref(:,k),convHullRef);
end

% plot IoU results
figure(2)
plot(iou,'LineWidth',2)

% end of function code
% *************************************************************************
%
%
% *************************************************************************
% end of document "main.m"
% *************************************************************************
