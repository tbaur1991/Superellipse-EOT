% generate random permutation of measurement set
idx = randperm(size(meas{k,1},2));

if strcmp(dim,'2D')
    % do measurement updates for superellipse shape tracking
    tic
    [X_UKF(:,k),P_UKF(:,:,k),meanInX1,meanInY1,varInX1,varInY1] = superellipse2D_UKF_update(X_UKF(:,k),P_UKF(:,:,k),sig_m,meas{k,1}(1:2,idx),...
                                                      wm,wc,lambda,art_noise,n_upd,meanInX1,meanInY1,varInX1,varInY1,tao,taper,source,nPoly);
    t1(k) = toc; t1(k) = t1(k)/size(meas{k,1},2);
    tic
    [X_S2KF(:,k),P_S2KF(:,:,k),meanInX2,meanInY2,varInX2,varInY2] = superellipse2D_S2KF_update(X_S2KF(:,k),P_S2KF(:,:,k),sig_m,meas{k,1}(1:2,idx),...
                                                     samples,numSamples,art_noise,n_upd,meanInX2,meanInY2,varInX2,varInY2,tao,taper,source,nPoly);
    t2(k) = toc; t2(k) = t2(k)/size(meas{k,1},2);
elseif strcmp(dim,'3D')
    % do measurement updates for superellipse shape tracking
    tic
    [X_3D(:,k),P_3D(:,:,k),meanInX1,meanInY1,varInX1,varInY1] = superellipse3D_UKF_update(X_3D(:,k),P_3D(:,:,k),sig_m,meas{k,1}(:,idx),wm_3D,wc_3D,...
                                              lambda_3D,art_noise,n_upd,meanInX1,meanInY1,varInX1,varInY1,tao,taper,source,nPoly,alpha,beta,kappa);
    t1(k) = toc; t1(k) = t1(k)/size(meas{k,1},2);
    tic
    [X_2D(:,k),P_2D(:,:,k),X_3DE(:,k),P_3DE(:,:,k),meanInX2,meanInY2,varInX2,varInY2] = superellipse3D_Decoup_update(X_2D(:,k),P_2D(:,:,k),X_3DE(:,k),...
                                              P_3DE(:,:,k),sig_m,meas{k,1}(:,idx),wm_2D,wc_2D,lambda_2D,wm_3DE,wc_3DE,lambda_3DE,art_noise,n_upd,...
                                              meanInX2,meanInY2,varInX2,varInY2,tao,taper,source,nPoly,alpha,beta,kappa);
    t2(k) = toc; t2(k) = t2(k)/size(meas{k,1},2);
end
