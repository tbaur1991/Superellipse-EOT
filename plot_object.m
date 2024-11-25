function plot_object(meas,X,convHullRef,poseBoat)
    % plot measurements
    m = plot(meas(1,:),meas(2,:),'.r','MarkerSize',5,'MarkerFaceColor','r');
    hold on
    axis equal
    xlabel('x/m')
    ylabel('y/m')

    % plot solgenia reference convex hull
    pos = poseBoat(1:2); or = poseBoat(5);
    R = [cos(or) -sin(or); sin(or) cos(or)];
    ref = R*convHullRef + pos;
    r = plot(ref(1,:),ref(2,:),'k','LineWidth',2);

    % plot estimated shape of superellipse UKF
    thetas = 0:0.1:2*pi+0.1;
    pos = X(1:2);
    or = X(4);
    a = positive_constrained(X(6));
    b = positive_constrained(X(7));
    epsilon = positive_constrained(X(8)) + 1;
    ty = 2/pi*atan(X(9));
    x = a*sign(cos(thetas)).*abs(cos(thetas)).^(2/epsilon);
    y = b*sign(sin(thetas)).*abs(sin(thetas)).^(2/epsilon);
    % tapering y
    y = (ty*x/a + 1).*y;
    % translate and rotate
    R = [cos(or) -sin(or); sin(or) cos(or)];
    pts = R*[x;y] + pos;
    s1 = plot(pts(1,:),pts(2,:),'LineWidth',2);
    
    % legend
    legend([m,r,s1],{'measurement','reference','Superellipse UKF'},'Location','northwest','Orientation','vertical');
    set(gcf,'Color','w')
    drawnow
    hold off
end