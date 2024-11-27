function iou = iouSuper(X,X_Ref,pts_ref,taper)
    % get shape points for estimate
    pts = getShapePoints2dTop(X,taper);

    % translate and rotate reference points
    pos_ref = X_Ref(1:2); or_ref = X_Ref(5);
    R = [cos(or_ref) -sin(or_ref); sin(or_ref ) cos(or_ref )];
    pts_ref = R*pts_ref + pos_ref;

    % calculate intersection over union using alpha shapes and convex hull
    iou = calculateIoU(pts, pts_ref);
end

function pts = getShapePoints2dTop(X,taper)
    % get pose and shape parameters
    pos = X(1:2); or = X(4); 
    a = positive_constrained(X(6)); 
    b = positive_constrained(X(7));
    e = positive_constrained(X(8)) + 1;

    % calculate shape points
    theta = 0:0.1:2*pi;
    pts(1,:) = a*sign(cos(theta)).*abs(cos(theta)).^(2/e);
    pts(2,:) = b*sign(sin(theta)).*abs(sin(theta)).^(2/e);
    if taper
        ty = 2/pi*atan(X(9));
        pts(2,:) = (ty*pts(1,:)/a + 1).*pts(2,:);
    end

    % translate and rotate points
    R = [cos(or) -sin(or); sin(or) cos(or)];
    pts = R*pts + pos;
end

function iou = calculateIoU(pts, pts_ref)
    % polyshapes
    poly = polyshape(pts(1,:),pts(2,:));
    poly_ref = polyshape(pts_ref(1,:),pts_ref(2,:));

    % intersection area
    poly_int = intersect(poly,poly_ref);
    area_int = area(poly_int);

    % union area
    poly_union = union(poly,poly_ref);
    area_union = area(poly_union);

    % intersection over union
    iou = area_int/area_union;
end