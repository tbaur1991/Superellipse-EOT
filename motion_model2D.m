function y = motion_model2D(x,T)
    % x(1) = x1 Position in x1 direction [m]
    % x(2) = x2 Position in x2 direction [m]
    % x(3) = v velocity in polar coordinates [m/s]
    % x(4) = h heading angle [rad]
    % x(5) = omega turn rate [rad/s]
    % x(6:end) = shape coefficients
    
    % Memory Allocation
    y = zeros(size(x));
    
    if ~x(5)==0
        y(1) = x(1) + 2*x(3)/x(5)*sin(x(5)*T/2)*cos(x(4) + x(5)*T/2); 
        y(2) = x(2) + 2*x(3)/x(5)*sin(x(5)*T/2)*sin(x(4) + x(5)*T/2);
    else
        y(1) = x(1) + cos(x(4))*x(3)*T;
        y(2) = x(2) + sin(x(4))*x(3)*T;
    end
    y(3) = x(3);
    y(4) = x(4) + x(5)*T;
    y(5:end) = x(5:end);
end