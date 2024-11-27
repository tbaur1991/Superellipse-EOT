if strcmp(dim,'2D')
    % do prediction steps for superellipse shape tracking
    [X_UKF(:,k),P_UKF(:,:,k)] = prediction_UKF_Super2D(X_UKF(:,k-1),P_UKF(:,:,k-1),wm,wc,lambda,dt,sigma_v,sigma_omega,taper);
    [X_S2KF(:,k),P_S2KF(:,:,k)] = prediction_UKF_Super2D(X_S2KF(:,k-1),P_S2KF(:,:,k-1),wm,wc,lambda,dt,sigma_v,sigma_omega,taper);
elseif strcmp(dim,'3D')
    % do prediction steps for superellipse shape tracking
    [X_3D(:,k),P_3D(:,:,k)] = prediction_UKF_Super3D(X_3D(:,k-1),P_3D(:,:,k-1),dt,sigma_v,sigma_omega,wm,wc,lambda,taper);
    [X_2D(:,k),P_2D(:,:,k),X_3DE(:,k),P_3DE(:,:,k)] = prediction_UKF_SuperDecoup(X_2D(:,k-1),P_2D(:,:,k-1),X_3DE(:,k-1),P_3DE(:,:,k-1),...
                                                        dt,sigma_v,sigma_omega,wm_2D,wc_2D,lambda_2D,taper);
    if do_plots
        % delete old axes
        delete(axo)
    end
end