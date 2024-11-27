% init figure
tcl = tiledlayout(1,1);
nexttile(tcl)

% plot measurements
m = plot(meas{k,1}(1,:),meas{k,1}(2,:),'.r','MarkerSize',5,'MarkerFaceColor','r');
hold on
axis equal
xlabel('x/m')
ylabel('y/m')

% plot solgenia reference convex hull
pos = X_Ref(1:2,k); or = X_Ref(5,k);
R = [cos(or) -sin(or); sin(or) cos(or)];
ref = R*convHullRef + pos;
r = plot(ref(1,:),ref(2,:),'k','LineWidth',2);

% plot estimated shape of superellipse UKF
thetas = 0:0.1:2*pi+0.1;
pos = X_UKF(1:2,k);
or = X_UKF(4,k);
a = positive_constrained(X_UKF(6,k));
b = positive_constrained(X_UKF(7,k));
e = positive_constrained(X_UKF(8,k)) + 1;
x = a*sign(cos(thetas)).*abs(cos(thetas)).^(2/e);
y = b*sign(sin(thetas)).*abs(sin(thetas)).^(2/e);
% tapering y
if taper
    ty = 2/pi*atan(X_UKF(9,k));
    y = (ty*x/a + 1).*y;
end
% translate and rotate
R = [cos(or) -sin(or); sin(or) cos(or)];
pts = R*[x;y] + pos;
s1 = plot(pts(1,:),pts(2,:),'LineWidth',2);

% plot estimated shape of superellipse S2KF
thetas = 0:0.1:2*pi+0.1;
pos = X_S2KF(1:2,k);
or = X_S2KF(4,k);
a = positive_constrained(X_S2KF(6,k));
b = positive_constrained(X_S2KF(7,k));
e = positive_constrained(X_S2KF(8,k)) + 1;
x = a*sign(cos(thetas)).*abs(cos(thetas)).^(2/e);
y = b*sign(sin(thetas)).*abs(sin(thetas)).^(2/e);
% tapering y
if taper
    ty = 2/pi*atan(X_S2KF(9,k));
    y = (ty*x/a + 1).*y;
end
% translate and rotate
R = [cos(or) -sin(or); sin(or) cos(or)];
pts = R*[x;y] + pos;
s2 = plot(pts(1,:),pts(2,:),'LineWidth',2);

% legend
leg = legend([m,r,s1,s2],{'Measurement','Reference','Superellipse UKF','Superellipse S2KF'},'Location','northwest','Orientation','vertical');
leg.Orientation = "horizontal";
leg.NumColumns = 2;
leg.Box = "off";
leg.Layout.Tile = "south";
set(gcf,'Color','w')
set(gcf,'Position',[680,546,560,332])
drawnow
hold off
