function ps = as_polygon3D(X,nPoly)
    % Returns the polygon chain that corresponds to the given x.
    a = half_constraint(X(7)); 
    b = half_constraint(X(8)); 
    e = half_constraint(X(10)) + 1;
    
    as = (0:360/(nPoly-1):360) * pi / 180;
    ps = [a;b].*sign([cos(as); sin(as)]).*abs([cos(as); sin(as)]).^(2/e);
end