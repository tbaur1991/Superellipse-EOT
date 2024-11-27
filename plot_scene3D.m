% init figure
tcl = tiledlayout(1,1);
nexttile(tcl)
ax = gca;

% plot measurements
mm = plot3(meas{k,1}(1,:),meas{k,1}(2,:),meas{k,1}(3,:),'.r','MarkerSize',5,'MarkerFaceColor','r');
hold on
axis equal
xlabel('x/m')
ylabel('y/m')
zlabel('z/m')

% plot solgenia reference convex hull
pos = X_Ref(1:2,k); or = X_Ref(5,k);
R = [cos(or) -sin(or); sin(or) cos(or)];
ref = R*convHullRef + pos;
X = repmat(ref(1,:),2,1);
Y = repmat(ref(2,:),2,1);
Z = [h_ref_min*ones(1,size(ref,2)); (h_ref_min+h_ref)*ones(1,size(ref,2))];
r = surf(X,Y,Z,'FaceAlpha',0.1,'FaceColor','#000000');

% plot estimate with UKF
syms theta u
pos = X_3D(1:3,k);
or = X_3D(5,k);
a = positive_constrained(X_3D(7,k));
b = positive_constrained(X_3D(8,k));
h = positive_constrained(X_3D(9,k));
e = positive_constrained(X_3D(10,k)) + 1;
% get coordinates
x = a*sign(cos(theta)).*abs(cos(theta)).^(2/e);
y = b*sign(sin(theta)).*abs(sin(theta)).^(2/e);
z = u*h;
% tapering x,y
if taper
    ty = 2/pi*atan(X_3D(11,k));
    y = (ty*x/a + 1).*y;
end
s1 = fsurf(x,y,z,[0 2*pi -1 1],'FaceAlpha',0.6,'FaceColor','#3CB043');
% translate surf plot
t = hgtransform('Parent',ax);
set(s1,'Parent',t);
m = makehgtform('translate',pos);
set(t,'Matrix',m)
% rotate surf plot
t = hgtransform('Parent',t);
set(s1,'Parent',t);
Rz = makehgtform('zrotate',or);
set(t,'Matrix',Rz);

% plot estimate with decoupled state
pos = [X_2D(1:2,k); X_3DE(1,k)];
or = X_2D(4,k);
a = positive_constrained(X_2D(6,k));
b = positive_constrained(X_2D(7,k));
h = positive_constrained(X_3DE(2,k));
e = positive_constrained(X_2D(8,k)) + 1;
% get coordinates
x = a*sign(cos(theta)).*abs(cos(theta)).^(2/e);
y = b*sign(sin(theta)).*abs(sin(theta)).^(2/e);
z = u*h;
% tapering x,y
if taper
    ty = 2/pi*atan(X_2D(9,k));
    y = (ty*x/a + 1).*y;
end
s2 = fsurf(x,y,z,[0 2*pi -1 1],'FaceAlpha',0.6,'FaceColor','#D5B85A');
% translate surf plot
t = hgtransform('Parent',ax);
set(s2,'Parent',t);
m = makehgtform('translate',pos);
set(t,'Matrix',m)
% rotate surf plot
t = hgtransform('Parent',t);
set(s2,'Parent',t);
Rz = makehgtform('zrotate',or);
set(t,'Matrix',Rz);

% legend
leg = legend([mm,r,s1,s2],{'Measurement','Reference','Superellipse 3D','Superellipse 3D Decoupled'},'Location','northwest','Orientation','vertical');
leg.Orientation = "horizontal";
leg.NumColumns = 2;
leg.Box = "off";
leg.Layout.Tile = "south";
set(gcf,'Color','w')
set(gcf,'Position',[680,546,560,332])
drawnow
hold off