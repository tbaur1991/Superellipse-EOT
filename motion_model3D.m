function y = motion_model3D(x,T)
    % x(1) = x1 Position in x1 direction [m]
    % x(2) = x2 Position in x2 direction [m]
    % x(3) = x3 Position in x3 direction [m]
    % x(4) = v velocity in polar coordinates [m/s]
    % x(5) = h heading angle [rad]
    % x(6) = omega turn rate [rad/s]
    % x(7:end) = shape coefficients
    
    % Memory Allocation
    y = zeros(size(x));
    
    if ~x(6)==0
        y(1) = x(1) + 2*x(4)/x(6)*sin(x(6)*T/2)*cos(x(5) + x(6)*T/2); 
        y(2) = x(2) + 2*x(4)/x(6)*sin(x(6)*T/2)*sin(x(5) + x(6)*T/2);
    else
        y(1) = x(1) + cos(x(5))*x(4)*T;
        y(2) = x(2) + sin(x(5))*x(4)*T;
    end
    y(3:4) = x(3:4);
    y(5) = x(5) + x(6)*T;
    y(6:end) = x(6:end);
end