function iou = iouSuper(X,X_Ref,pts_ref)
    % get shape points for estimate
    pts_est = getShapePoints2dTop(X);

    % translate and rotate reference points
    pos_ref = X_Ref(1:2); or_ref = X_Ref(5);
    R = [cos(or_ref) -sin(or_ref); sin(or_ref ) cos(or_ref )];
    pts_ref = R*pts_ref + pos_ref;

    % calculate intersection over union using convex hull reference
    iou = calculateIoU(pts_est, pts_ref);
end

function pts = getShapePoints2dTop(X)
    % get pose and shape parameters
    pos = X(1:2); or = X(4); a = positive_constrained(X(6)); b = positive_constrained(X(7));
    e = positive_constrained(X(8)) + 1; ty = 2/pi*atan(X(9));

    % calculate shape points
    theta = 0:0.1:2*pi;
    pts(1,:) = a*sign(cos(theta)).*abs(cos(theta)).^(2/e);
    pts(2,:) = b*sign(sin(theta)).*abs(sin(theta)).^(2/e);
    pts(2,:) = (ty*pts(1,:)/a + 1).*pts(2,:);

    % translate and rotate points
    R = [cos(or) -sin(or); sin(or) cos(or)];
    pts = R*pts + pos;
end

function iou = calculateIoU(pts_est, pts_ref)
    % polyshapes
    poly_est = polyshape(pts_est(1,:),pts_est(2,:));
    poly_ref = polyshape(pts_ref(1,:),pts_ref(2,:));

    % intersection area
    poly_int = intersect(poly_est,poly_ref);
    area_int = area(poly_int);

    % union area
    poly_union = union(poly_est,poly_ref);
    area_union = area(poly_union);

    % intersection over union
    iou = area_int/area_union;
end