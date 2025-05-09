clear;
close all;
%% One-time computation and generation of ZZ matrix
meshN = 20;
urange = linspace(0,2*pi,meshN);
urange = urange(1:end-1);
vrange = linspace(0,pi,meshN);
Z = zeros(3,numel(vrange),numel(urange));
i=1;
for u=urange
    Z(:,:,i) = [cos(u).*sin(vrange);sin(u).*sin(vrange);cos(vrange)];
    i = i+1;
end

ZZ = reshape(Z,3,numel(vrange)*numel(urange));
ZZ(:, [meshN*(1:8)+1, meshN*(2:9)]) = []; % removing duplicate points

%% 1st plot
[X,k,Z] = cellipse([0;0;0], eye(3), 3);

ts = trisurf(k,X(:,1),X(:,2),X(:,3),...
    'Facecolor','red','FaceAlpha',0.5,'LineStyle','-');
axis equal; hold on;
set(gca, 'box', 'on', 'BoxStyle', 'full');
view(3);

%% 2nd plot
P = [2.9087, 2.4783, 1.1200;
     2.4783, 6.9100, 3.0680;
     1.1200, 3.0680, 4.5327];
mu = [1;2;3];
c = 1;

[X,k] = cellipse(mu,P,c,Z);
pause(1);
set(ts, 'Faces', k, 'Vertices', X);

%% Nth plot
n = 10;
for i=1:n
    % generate some random ellipsoids
    [V,~] = qr(randi(10, 3, 3));
    L = diag(5*rand(1,3));
    P = V*L*V.';
    mu = 5*rand(3,1) - 2.5;
    [X,k] = cellipse(mu,P,c,Z);
    set(ts, 'Faces', k, 'Vertices', X);
    pause(.4);

end



