%% Investigating various ways to draw an ellipse
%% Using fsurf (slow)
clear;
close all;

P = [2.9087, 2.4783, 1.1200;
     2.4783, 6.9100, 3.0680;
     1.1200, 3.0680, 4.5327];
mu = [1;2;3];
c = 1;
[V,L] = eig(P);
Ri = V*sqrt(L);

fx = @(u,v) sqrt(c) * (Ri(1,1) * cos(u).*sin(v) + Ri(1,2) * sin(u).*sin(v) + Ri(1,3) * cos(v)) + mu(1);
fy = @(u,v) sqrt(c) * (Ri(2,1) * cos(u).*sin(v) + Ri(2,2) * sin(u).*sin(v) + Ri(2,3) * cos(v)) + mu(2);
fz = @(u,v) sqrt(c) * (Ri(3,1) * cos(u).*sin(v) + Ri(3,2) * sin(u).*sin(v) + Ri(3,3) * cos(v)) + mu(3);

fs = fsurf(fx,fy,fz,[0 2*pi 0 pi], 'FaceAlpha', .6, 'LineStyle', 'none');
axis equal;
set(gca, 'box', 'on', 'BoxStyle', 'full');

% P = eye(3)
% mu = [3;-2;-3];
% c = 1;
% [V,L] = eig(P);
% Ri = V*sqrt(L);
% 
% fx = @(u,v) sqrt(c) * (Ri(1,1) * cos(u).*sin(v) + Ri(1,2) * sin(u).*sin(v) + Ri(1,3) * cos(v)) + mu(1);
% fy = @(u,v) sqrt(c) * (Ri(2,1) * cos(u).*sin(v) + Ri(2,2) * sin(u).*sin(v) + Ri(2,3) * cos(v)) + mu(2);
% fz = @(u,v) sqrt(c) * (Ri(3,1) * cos(u).*sin(v) + Ri(3,2) * sin(u).*sin(v) + Ri(3,3) * cos(v)) + mu(3);
% 
% 
% set(fs, 'XFunction', fx, 'YFunction', fy, 'ZFunction', fz);

%% Using Points
figure;
P = [2.9087, 2.4783, 1.1200;
     2.4783, 6.9100, 3.0680;
     1.1200, 3.0680, 4.5327];
mu = [1;2;3];
c = 1;
[V,L] = eig(P);
Ri = V*sqrt(L);

urange = linspace(0,2*pi,10);
urange = urange(1:end-1);
vrange = linspace(0,pi,10);
Z = zeros(3,numel(vrange),numel(urange));
i=1;
for u=urange
    Z(:,:,i) = [cos(u).*sin(vrange);sin(u).*sin(vrange);cos(vrange)];
    i = i+1;
end

ZZ = reshape(Z,3,numel(vrange)*numel(urange));
ZZ(:, [10*(1:8)+1, 10*(2:9)]) = [];
X = (sqrt(c) * Ri * ZZ) + mu;
XX = X.';

plot3(X(1,:), X(2,:), X(3,:), 'bo');
axis equal;
set(gca, 'box', 'on', 'BoxStyle', 'full');
%%
DT = delaunayTriangulation(XX(:,1), XX(:,2), XX(:,3));
figure;
tetramesh(DT,'FaceAlpha',0.3, 'LineStyle', 'none');
figure;
[K,v] = convexHull(DT);
trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3))
axis equal;
set(gca, 'box', 'on', 'BoxStyle', 'full');

%%
figure;
shp = alphaShape(XX(:,1), XX(:,2), XX(:,3),1);
plot(shp)
axis equal
%%
figure;
k = boundary(XX);
hold on
ts = trisurf(k,XX(:,1),XX(:,2),XX(:,3),'Facecolor','red','FaceAlpha',0.1,'LineStyle','none')

axis equal;
set(gca, 'box', 'on', 'BoxStyle', 'full');
view(3);
%%
figure;
k = convhull(XX);
ts = trisurf(k,XX(:,1),XX(:,2),XX(:,3),'Facecolor','red','FaceAlpha',0.1,'LineStyle','none')
axis equal;
set(gca, 'box', 'on', 'BoxStyle', 'full');
view(3);
