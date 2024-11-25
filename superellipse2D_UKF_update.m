function [X,P,varIn] = superellipse2D_UKF_update(X,P,sig_r,Y,wm,wc,lambda,n_upd,varIn,tao)
    fprintf(['\n*** call of function: ',mfilename,' ***\n'])
    
    % get dimension of system state
    nx = length(X);
    
    % iterate n times over n_upd measurements
    for n = 1:size(Y,2)/n_upd
        % calculate sigma points
        A = sqrt(nx + lambda) * chol(P)';
        sig_points = [zeros(size(X)) -A A];
        sig_points = sig_points + repmat(X, 1, size(sig_points, 2));
        numSamples = size(sig_points,2);
    
        % get new measurement set
        meas = Y(:,n_upd*n-n_upd+1:n_upd*n);
        
        % predict sigma points for every measurement
        z_predict = zeros(2*n_upd,numSamples);
        for i = 1:numSamples
            % get sample quantities
            pos = sig_points(1:2,i);
            or = sig_points(4,i);
            a = positive_constrained(sig_points(6,i)); 
            b = positive_constrained(sig_points(7,i)); 
            e = positive_constrained(sig_points(8,i)) + 1;
            ty = 2/pi*atan(sig_points(9,i)); 
            % transform measurements to local coordinate system
            R = [cos(or) -sin(or); sin(or) cos(or)];
            meas_loc = R'*(meas - pos);
            % do inverse tapering transformation
            meas_loc(2,:) = meas_loc(2,:)./(ty*meas_loc(1,:)./a + 1);
            % do sample prediction with ray function
            xi = sign(meas_loc(1,:)).*power(complex(1./(1/(abs(a)^e) + abs(power(complex((meas_loc(2,:)./(meas_loc(1,:)*b))),e)))),1/e);
            yi = xi.*(meas_loc(2,:)./meas_loc(1,:));
            zs = [xi; yi];
            z_predict(:,i) = real(zs(:) - meas_loc(:));
        end
        
        % get predicted measurements
        z_pred = sum(wm.*z_predict,2);
    
        % calculate measurement noise covariance matrix for asymmetric noise
        % first step measurements in local coordinates
        pos = X(1:2); or = X(4); 
        R = [cos(or) -sin(or); sin(or) cos(or)]; 
        a = positive_constrained(X(6)); b = positive_constrained(X(7)); 
        e = positive_constrained(X(8)) + 1; ty = 2/pi*atan(X(9));
        meas_loc = R'*(meas - pos);
        meas_loc(2,:) = meas_loc(2,:)./(ty*meas_loc(1,:)./a + 1);
        % second step measurements inside or outside
        inside_outside = abs(meas_loc(1,:)/a).^e + abs(meas_loc(2,:)/b).^e - 1;
        % indices of measurements inside
        idx = find(inside_outside < 0);
        % update estimate of inside measurement variance
        varIn = 1/(1 + length(idx)/tao)*varIn + 1/(tao + 2*length(idx))*(sum(z_pred(2*idx-1).^2) + sum(z_pred(2*idx).^2));
        % build covariance matrix
        R = sig_r^2*ones(1,2*n_upd);
        R(2*idx-1) = varIn; R(2*idx) = varIn; R = diag(R);
        
        % get update matrices
        S = zeros(2*n_upd);
        Psi = zeros(nx,2*n_upd);
        for i = 1:size(sig_points,2)
            S = S + wc(i)*(z_predict(:,i) - z_pred)*(z_predict(:,i) - z_pred)';
            Psi = Psi + wc(i)*(sig_points(1:nx,i) - X)*(z_predict(:,i) - z_pred)';
        end
        S = S + R;
        
        % do UKF measurement update
        K = Psi/S;
        X = X + K*(zeros(2*n_upd,1) - z_pred);
        P = P - K*Psi';
    end

    fprintf(['*** return of function: ',mfilename,' ***\n'])
end