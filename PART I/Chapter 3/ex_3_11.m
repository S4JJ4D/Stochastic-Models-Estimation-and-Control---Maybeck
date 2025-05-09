%%
clear;
close all;

mu = [0 0];
% eigenvectors
v1 = [2;1]/vecnorm([2;1]);
v2 = [-1;2]/vecnorm([-1;2]);
V = [v1,v2];
D = diag([4 2]);
Di = inv(D);

D12 = sqrt(D);
D12_i = inv(D12);

Q = V*D12;
Qi = D12_i * V.';
P = V*D*V.';
Pi = V*Di*V.';

f = @(x, y, mu, P, Pi) 1/(2*pi*sqrt(det(P))) * ...
    exp(-1/2 * ( (x-mu(1)).^2 * Pi(1,1) + (y-mu(2)).^2 * Pi(2,2) + 2*(x-mu(1)).*(y-mu(2))*Pi(1,2) ));

f_level = @(c, P) 1/(2*pi*sqrt(det(P))) * exp(-1/2 * c^2);

a = 2;
nPoints = 100;
xrange = linspace(mu(1) - a*P(1), mu(1) + a*P(1), nPoints);
yrange = linspace(mu(2) - a*P(2,2), mu(2) + a*P(2,2), nPoints);

[X,Y] = meshgrid(xrange, yrange);
Z = f(X,Y, mu, P, Pi);

surfPlt = surf(X,Y,Z, 'FaceAlpha', .5, 'LineStyle', 'none');
daspect([1 1 .02])

hold on;

c = 1;
levels = [f_level(c, P), f_level(c, P)];
contour(X,Y,Z,levels)


hold on;

quiver(0,0,v1(1),v1(2),'AutoScale','off');
quiver(0,0,v2(1),v2(2),'AutoScale','off');

quiver(0,0,Q(1,1),Q(2,1),'AutoScale','off');
quiver(0,0,Q(1,2),Q(2,2),'AutoScale','off');


box on;


Points = [];
for x=xrange
    for y=yrange
        p = [x;y];
        if p.'*Pi*p < c^2 
            Points = [Points;p.'];
        end
    end
end

plot(Points(:,1), Points(:,2), 'bo');

% axis equal;

%%
dt = delaunayTriangulation(Points)
IC = incenter(dt);
triplot(dt)
hold on
plot(IC(:,1),IC(:,2),'*r')

cl = dt.ConnectivityList;
triCount = size(cl,1);

areaVec = zeros(triCount, 1);
for i=1:triCount
    areaVec(i) = 1/2 * abs(det(...
        [1,1,1;Points(cl(i,1), 1), Points(cl(i,2), 1), Points(cl(i,3), 1); ...
        Points(cl(i,1), 2), Points(cl(i,2), 2), Points(cl(i,3), 2)]));
end

f2 = @(x,y) 1/(2*pi) * exp(-1/2 *(x.^2 + y.^2));

z = arrayfun(@(x,y) f(x,y, mu, P, Pi), IC(:,1), IC(:,2), 'UniformOutput', true);
integral_val = areaVec.' * z

z1 = arrayfun(@(x,y) f(x,y, mu, P, Pi), dt.Points(:,1), dt.Points(:,2), 'UniformOutput', true);
int1 = integrateTriangulation(dt, z1)

fh = @(c) 1 - exp(-1/2 * c.^2);
fh(c)
%%
function int = integrateTriangulation(trep, z)
P = trep.Points; T = trep.ConnectivityList;
d21 = P(T(:,2),:)-P(T(:,1),:);
d31 = P(T(:,3),:)-P(T(:,1),:);
areas = abs(1/2*(d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1)));
int = areas.' * mean(z(T),2);
end