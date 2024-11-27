% init orientation
c = cov(meas{k,1}(1,:),meas{k,1}(2,:));
[v,d] = eig(c); [~,i] = max(diag(d));
or = atan2(v(2,i),v(1,i));
% init position
pos_2D = mean(meas{k,1}(1:2,:),2);
if strcmp(dim,'2D')
    % do measurement transformation to local coordinate
    R_rot = [cos(or) -sin(or); sin(or) cos(or)];
    mm = R_rot'*(meas{k,1}(1:2,:) - pos_2D);
elseif strcmp(dim,'3D')
     % do position initialization and measurement transformation to local coordinate
    pos = zeros(3,1);
    pos(1:2) = pos_2D; 
    pos(3) = min(meas{k,1}(3,:));
    R_rot = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
    mm = R_rot'*(meas{k,1} - pos);
end
% init 2D extent
a = max(abs(mm(1,:)));
b = abs(max(mm(2,:)));

% init 2D superellipse shape tracking system state
if strcmp(dim,'2D')
    % init UKF state
    if taper
        X_UKF(:,k) = [pos_2D;0;or;0;a-1;b-1;3;0];
        P_UKF(:,:,k) = blkdiag(c,2,10*pi/180*eye(2),0.5*eye(2),1,1/3);
    else
        X_UKF(:,k) = [pos_2D;0;or;0;a-1;b-1;3];
        P_UKF(:,:,k) = blkdiag(c,2,10*pi/180*eye(2),0.5*eye(2),1);
    end
    % init S2KF state
    X_S2KF(:,k) = X_UKF(:,k); P_S2KF(:,:,k) = P_UKF(:,:,k);
% init 3D superellipse shape tracking system state
elseif strcmp(dim,'3D')
    h = max(mm(3,:)) - min(mm(3,:));
    if taper
        X_3D(:,k) = [pos(1:2);pos(3)+h/2;0;or;0;a-1;b-1;(h/2)-1;3;0];
        P_3D(:,:,k) = blkdiag(c,2,2,10*pi/180*eye(2),0.5*eye(3),1,1/3);
        X_2D(:,k) = [pos(1:2);0;or;0;a-1;b-1;3;0];
        P_2D(:,:,k) = blkdiag(c,2,10*pi/180*eye(2),0.5*eye(2),1,1/3);
    else
        X_3D(:,k) = [pos(1:2);pos(3)+h/2;0;or;0;a-1;b-1;(h/2)-1;3];
        P_3D(:,:,k) = blkdiag(c,2,2,10*pi/180*eye(2),0.5*eye(3),1);
        X_2D(:,k) = [pos(1:2);0;or;0;a-1;b-1;3];
        P_2D(:,:,k) = blkdiag(c,2,10*pi/180*eye(2),0.5*eye(2),1);
    end
    X_3DE(:,k) = [pos(3)+h/2;(h/2)-1];
    P_3DE(:,:,k) = diag([2,0.5]);
end