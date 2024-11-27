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
% Creation date: 26.11.2024
% *************************************************************************
% Content: In this file, the simulation environment and estimation 
% procedure for a superellipse shape tracking in 2D or 3D using the research 
% vessel 'Solgenia' is handled. Measurements are generated from a simulated
% LiDAR sensor: https://de.mathworks.com/help/uav/ref/uavlidarpointcloudgenerator-system-object.html
%
% !! The toolbox for the S2KF must be downloaded from
% https://nonlinearestimation.bitbucket.io/ !!
% *************************************************************************

%% start of function code
close all; clear; clc

% load reference trajectory and simulated measurement sources
load('ref_trajectory.mat')
load('solgenia_sim_meas.mat')

%% settings to be chosen in this simulation

% choose if '2D' or '3D' models are applied
dim = '3D';

% choose artificieal measurement noise for processing a mixture of boundary
% and interior measurements to be used or not.
art_noise = true;

% choose if tapering coefficient is used for superellipse model
taper = true;

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
        state_initialization
    else
        % do prediction step for superellipse shape tracking
        prediction_steps
    end

    % do measurement update for superellipse shape tracking
    update_steps

    % plot measurements, reference and estimate
    if do_plots
        if strcmp(dim,'2D')
            plot_scene2D
        elseif strcmp(dim,'3D')
            plot_scene3D
            axo = ax;
        end
    end

    % calculate intersection over union using reference convex hull
    if strcmp(dim,'2D')
        iou1(k) = iouSuper(X_UKF(:,k),X_Ref(:,k),convHullRef,taper);
        iou2(k) = iouSuper(X_S2KF(:,k),X_Ref(:,k),convHullRef,taper);
    elseif strcmp(dim,'3D')
        iou1(k) = iouSuper3D(X_3D(:,k),[X_Ref(:,k);h_ref;h_ref_min],convHullRef,taper);
        iou2(k) = iouSuper3D([X_2D(1:2,k);X_3DE(1,k);X_2D(3:7,k);X_3DE(2,k);X_2D(8:end,k)],[X_Ref(:,k);h_ref;h_ref_min],convHullRef,taper);
    end
end

% plot results
figure(2)
tcl = tiledlayout(2,1);
% plot IoUs
nexttile(tcl)
plot(1:nSamples,iou1,'LineWidth',1.5)
hold on
plot(1:nSamples,iou2,'LineWidth',1.5)
% plot computation times
nexttile(tcl)
semilogy(1:nSamples,t1,'LineWidth',1.5)
hold on
semilogy(1:nSamples,t2,'LineWidth',1.5)
% legend
if strcmp(dim,'2D')
    leg = legend('UKF','S2KF');
elseif strcmp(dim,'3D')
    leg = legend('UKF','UKF Decoupled');
end
leg.Orientation = "horizontal";
leg.Box = "off";
leg.Layout.Tile = "south";
set(gcf,'Color','w')

% end of function code
% *************************************************************************
%
%
% *************************************************************************
% end of document "main.m"
% *************************************************************************
