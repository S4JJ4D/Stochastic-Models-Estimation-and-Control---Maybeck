cla;
clear;
n = 10;

% 1
% P = rand([n 2]);

% 2
% xrange = linspace(0,1,n);
% yrange = linspace(0,1,n);
% [X,Y] = meshgrid(xrange,yrange);
% P = [X(:),Y(:)];

% 3
P = [...
    0.2642    0.7592;
    0.7673    0.7086;
    0.8350    0.8668;
    0.9178    0.8237;
    0.0960    0.7218;
    0.1774    0.5595;
    0.7936    0.9881;
    0.7489    0.9099;
    0.0547    0.9393;
    0.6602    0.4028;
];


f1 = @(x,y) 1*(x.^2 + y.^2);
fsurf(f1, [0, 1, 0, 1], 'FaceAlpha', .5, 'LineStyle', 'none')
hold on;

dt = delaunayTriangulation(P)
IC = incenter(dt);
triplot(dt)
hold on
plot(IC(:,1),IC(:,2),'*r')

for i=1:n
    text(P(i,1), P(i,2), {sprintf('%i', i),''}, ...
        "FontSize", 12, 'FontWeight', 'bold');
end


cl = dt.ConnectivityList;
triCount = size(cl,1);

areaVec = zeros(triCount, 1);
for i=1:triCount
    areaVec(i) = 1/2 * abs(det(...
        [1,1,1;P(cl(i,1), 1), P(cl(i,2), 1), P(cl(i,3), 1); ...
        P(cl(i,1), 2), P(cl(i,2), 2), P(cl(i,3), 2)]));
end

z = arrayfun(f1, IC(:,1), IC(:,2))
integral_val = areaVec.' * z
