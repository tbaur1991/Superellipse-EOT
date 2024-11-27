% *************************************************************************
%    I   SSSS   DDDD
%    I   S      D   D 
%    I   SSSS   D   D    Institut fur Systemdynamik
%    I      S   D   D
%    I   SSSS   DDDD   
% *************************************************************************
% Author: 
% University: University of Applied Sciences "HTWG Konstanz"
% E-Mail: 
% Creation date: 
% File name: function_template.m
% *************************************************************************
% Content:
% 
% Input:
%
% Output:
%
% *************************************************************************

%% start of function code
function [X_2D,P_2D,X_3DE,P_3DE,meanInX,meanInY,varInX,varInY] = superellipse3D_Decoup_update(X_2D,P_2D,X_3DE,P_3DE,sig_m,Y,wm_2D,wc_2D,lambda_2D,wm_3DE,wc_3DE,lambda_3DE,art_noise,n_upd,meanInX,meanInY,varInX,varInY,tao,taper,source,nPoly,alpha,beta,kappa)
    % get state dimensions
    nx_2D = length(X_2D); nx_3DE = length(X_3DE);
    nx_3DE = nx_3DE + n_upd;

    % check if number of measurements smaller than n_upd
    nMeas = size(Y,2);
    if nMeas < n_upd
        % change number of measurements to be updated and UKF weights
        n_upd = nMeas;
        nx_3DE = length(X_3DE) + n_upd;
        lambda_3DE = alpha^2*(nx_3DE + kappa) - nx_3DE;
        wm_3DE = zeros(1,2*nx_3DE+1); wc_3DE = wm_3DE;
        wm_3DE(1) = lambda_3DE/(nx_3DE + lambda_3DE);
        wc_3DE(1) = lambda_3DE/(nx_3DE + lambda_3DE) + (1 - alpha^2 + beta);
        wm_3DE(2:2*nx_3DE + 1) = 1/(2*(nx_3DE + lambda_3DE));
        wc_3DE(2:2*nx_3DE + 1) = 1/(2*(nx_3DE + lambda_3DE));
    end

    % iterate j times over n_upd measurements
    for n = 1:ceil(nMeas/n_upd)
        % update 2D state using UKF
        % calculate sigma points
        A = sqrt(nx_2D + lambda_2D) * chol(P_2D)';
        sig_points = [zeros(size(X_2D)) -A A];
        sig_points = sig_points + repmat(X_2D, 1, size(sig_points, 2));
        samplesUKF = size(sig_points,2);

        % get new measurement set
        if n <= floor(nMeas/n_upd)
            meas_2D = Y(1:2,n_upd*n-n_upd+1:n_upd*n);
            meas_3D = Y(3,n_upd*n-n_upd+1:n_upd*n);
        else
            meas_2D = Y(1:2,end-mod(nMeas,n_upd)+1:end);
            meas_3D = Y(3,end-mod(nMeas,n_upd)+1:end);
            % change number of measurements to be updated and UKF weights
            n_upd = mod(nMeas,n_upd);
            nx_3DE = length(X_3DE) + n_upd;
            lambda_3DE = alpha^2*(nx_3DE + kappa) - nx_3DE;
            wm_3DE = zeros(1,2*nx_3DE+1); 
            wc_3DE = wm_3DE;
            wm_3DE(1) = lambda_3DE/(nx_3DE + lambda_3DE);
            wc_3DE(1) = lambda_3DE/(nx_3DE + lambda_3DE) + (1 - alpha^2 + beta);
            wm_3DE(2:2*nx_3DE + 1) = 1/(2*(nx_3DE + lambda_3DE));
            wc_3DE(2:2*nx_3DE + 1) = 1/(2*(nx_3DE + lambda_3DE));
        end

        % predict sigma points for every measurement
        z_predict = zeros(2*n_upd,samplesUKF);
        for i = 1:samplesUKF
            % get sample variables
            pos = sig_points(1:2,i);
            or = sig_points(4,i);
            a = positive_constrained(sig_points(6,i)); 
            b = positive_constrained(sig_points(7,i)); 
            e = positive_constrained(sig_points(8,i)) + 1; 
            % transform measurements to local coordinate system
            R = [cos(or) -sin(or); sin(or) cos(or)];
            meas_loc = R'*(meas_2D - pos);
            % do inverse tapering transformation
            if taper
                ty = 2/pi*atan(sig_points(9,i));
                meas_loc(2,:) = meas_loc(2,:)./(ty*meas_loc(1,:)./a + 1);
            end
            % do sample prediction with ray function or polygonal chain
            % approximation
            if strcmp(source,'radial')
                % calculate intersection point
                xi = sign(meas_loc(1,:)).*power(complex(1./(1/(abs(a)^e) + abs(power(complex((meas_loc(2,:)./(meas_loc(1,:)*b))),e)))),1/e);
                yi = xi.*(meas_loc(2,:)./meas_loc(1,:));
                zs = [xi; yi];
            elseif strcmp(source,'projected')
                % calculate polygonal chain and projection point
                ps = as_polygon(sig_points(:,i),nPoly);
                zs = project(ps, meas_loc(1:2,:)); 
            else
                error('Wrong definition of variable source. Choose either radial or projected')
            end 
            z_predict(:,i) = real(zs(:) - meas_loc(:));
        end

        % get predicted measurements
        z_pred = sum(wm_2D.*z_predict,2);

        % calculate measurement noise covariance matrix
        if art_noise
            % for asymmetric noise
            % first step measurements in local coordinates
            pos = X_2D(1:2); or = X_2D(4); 
            R = [cos(or) -sin(or); sin(or) cos(or)]; 
            a = positive_constrained(X_2D(6)); 
            b = positive_constrained(X_2D(7)); 
            e = positive_constrained(X_2D(8)) + 1; 
            meas_loc = R'*(meas_2D - pos);
            if taper
                ty = 2/pi*atan(X_2D(9));
                meas_loc(2,:) = meas_loc(2,:)./(ty*meas_loc(1,:)./a + 1);
            end
            % second step measurements inside or outside
            inside_outside = abs(meas_loc(1,:)/a).^e + abs(meas_loc(2,:)/b).^e - 1;
            % indices of measurements inside
            idx = find(inside_outside < 0);
            % update estimate of inside measurement mean
            meanInX = 1/(1 + length(idx)/tao)*meanInX + 1/(tao + length(idx))*sum(z_pred(2*idx-1));
            meanInY = 1/(1 + length(idx)/tao)*meanInY + 1/(tao + length(idx))*sum(z_pred(2*idx));
            z_pred(2*idx-1) = z_pred(2*idx-1) + meanInX; z_pred(2*idx) = z_pred(2*idx) + meanInY;
            % update estimate of inside measurement variance
            varInX = 1/(1 + length(idx)/tao)*varInX + 1/(tao + length(idx))*sum((z_pred(2*idx-1) - meanInX).^2);
            varInY = 1/(1 + length(idx)/tao)*varInY + 1/(tao + length(idx))*sum((z_pred(2*idx) - meanInY).^2);
            % build covariance matrix
            R = sig_m^2*ones(1,2*n_upd);
            R(2*idx-1) = varInX; R(2*idx) = varInY; R = diag(R);
        else
            R = sig_m^2*eye(2*n_upd);
        end

        % get update matrices
        S = zeros(2*n_upd);
        Psi = zeros(nx_2D,2*n_upd);
        for i = 1:size(sig_points,2)
            S = S + wc_2D(i)*(z_predict(:,i) - z_pred)*(z_predict(:,i) - z_pred)';
            Psi = Psi + wc_2D(i)*(sig_points(:,i) - X_2D)*(z_predict(:,i) - z_pred)';
        end
        S = S + R;

        % do UKF measurement update
        K = Psi/S;
        X_2D = X_2D + K*(zeros(2*n_upd,1) - z_pred);
        P_2D = P_2D - K*Psi';
        P_2D = 0.5*(P_2D+P_2D');

        % update line state using UKF
        % extend system state by mean of u for sample transformation
        Xu = [X_3DE; zeros(n_upd,1)];
        L = blkdiag(chol(P_3DE)',eye(n_upd));

        % calculate sigma points
        A = sqrt(nx_3DE + lambda_3DE) * L;
        sig_points = [zeros(size(Xu)) -A A];
        sig_points = sig_points + repmat(Xu, 1, size(sig_points, 2));
        sig_points(3:2+n_upd,:) = 0.5*(1 + erf(sig_points(3:2+n_upd,:)/sqrt(2)));
        nSamples = size(sig_points,2);
        
        % predict sigma points for every measurement
        z_predict = zeros(n_upd,nSamples);
        for i = 1:nSamples
            % get sample variables
            z = sig_points(1,i);
            h = positive_constrained(sig_points(2,i));
            us = sig_points(3:end,i);
            z_predict(:,i) = abs(meas_3D' - z) - us*h;
        end

        % get predicted measurements
        z_pred = sum(wm_3DE.*z_predict,2);

        % get update matrices
        S = zeros(n_upd);
        Psi = zeros(2,n_upd);
        for i = 1:size(sig_points,2)
            S = S + wc_3DE(i)*(z_predict(:,i) - z_pred)*(z_predict(:,i) - z_pred)';
            Psi = Psi + wc_3DE(i)*(sig_points(1:2,i) - X_3DE)*(z_predict(:,i) - z_pred)';
        end
        S = S + sig_m^2*eye(n_upd);

        % do UKF measurement update
        K = Psi/S;
        X_3DE = X_3DE + K*(zeros(n_upd,1) - z_pred);
        P_3DE = P_3DE - K*Psi';
        P_3DE = 0.5*(P_3DE+P_3DE');
    end
end
% end of function code
% *************************************************************************
%
%
% *************************************************************************
% end of document "function_template.m"
% *************************************************************************
