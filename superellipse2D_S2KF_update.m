function [X,P,meanInX,meanInY,varInX,varInY] = superellipse2D_S2KF_update(X,P,sig_r,Y,samples,numSamples,art_noise,n_upd,meanInX,meanInY,varInX,varInY,tao,taper,source,nPoly)
    % get dymension of system state
    nx = length(X);

    % check if number of measurements smaller than n_upd
    nMeas = size(Y,2);
    if nMeas < n_upd
        n_upd = nMeas;
    end
    
    % iterate j times over n_upd measurements
    for n = 1:ceil(nMeas/n_upd)
        % calculate samples using previous update
        Xu = X; L = chol(P)';
        samples_upd = L*samples + Xu;

        % get new measurement set
        if n <= floor(nMeas/n_upd)
            meas = Y(:,n_upd*n-n_upd+1:n_upd*n);
        else
            meas = Y(:,end-mod(nMeas,n_upd)+1:end);
            n_upd = mod(nMeas,n_upd);
        end
    
        % predict samples for every measurement
        z_predict = zeros(2*n_upd,numSamples);
        for i = 1:numSamples
            % get sample variables
            pos = samples_upd(1:2,i);
            or = samples_upd(4,i);
            a = positive_constrained(samples_upd(6,i)); 
            b = positive_constrained(samples_upd(7,i)); 
            e = positive_constrained(samples_upd(8,i)) + 1;
            % transform measurements to local coordinate system
            R = [cos(or) -sin(or); sin(or) cos(or)];
            meas_loc = R'*(meas - pos);
            if taper
                % do tapering transformation
                ty = 2/pi*atan(samples_upd(9,i));
                meas_loc(2,:) = meas_loc(2,:)./(ty*meas_loc(1,:)/a + 1);
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
                ps = as_polygon(samples_upd(:,i),nPoly);
                zs = project(ps, meas_loc); 
            else
                error('Wrong definition of variable source. Choose either radial or projected')
            end
            z_predict(:,i) = real(zs(:) - meas_loc(:));
        end
        
        % get predicted measurement
        z_pred = mean(z_predict, 2);

        % calculate measurement noise covariance matrix
        if art_noise
            % for artificial noise
            % first step measurements in local coordinates
            pos = X(1:2); or = X(4); 
            R = [cos(or) -sin(or); sin(or) cos(or)]; 
            a = positive_constrained(X(6)); 
            b = positive_constrained(X(7)); 
            e = positive_constrained(X(8)) + 1; 
            meas_loc = R'*(meas - pos);
            if taper
                ty = 2/pi*atan(X(9));
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
            R = sig_r^2*ones(1,2*n_upd);
            R(2*idx-1) = varInX; R(2*idx) = varInY; R = diag(R);
        else
            R = sig_r^2*eye(2*n_upd);
        end

        % get state measurement cross-covariance matrix
        z_p = z_predict - z_pred;
        YY = z_p*z_p'/numSamples + R;
        C = (samples_upd(1:nx,:) - X)*z_p'/numSamples;
    
        % do S2KF measurement update
        K = C/YY;
        X = X + K*(zeros(2*n_upd,1) - z_pred);
        P = P - K*C';
        P = 0.5*(P+P');
    end
end