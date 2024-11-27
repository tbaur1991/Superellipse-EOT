function ps = as_polygon(X,nPoly)
    % Returns the polygon chain that corresponds to the given x.
    a = positive_constrained(X(6)); 
    b = positive_constrained(X(7)); 
    e = positive_constrained(X(8)) + 1;
    
    as = (0:360/(nPoly-1):360) * pi / 180;
    ps = [a;b].*sign([cos(as); sin(as)]).*abs([cos(as); sin(as)]).^(2/e);
end