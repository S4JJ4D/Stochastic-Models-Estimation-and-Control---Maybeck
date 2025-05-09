%% performance comparison of bounding objects computation

clear;
close all;
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
%%
P = [2.9087, 2.4783, 1.1200;
     2.4783, 6.9100, 3.0680;
     1.1200, 3.0680, 4.5327];
mu = [1;2;3];
c = 1;
[V,L] = eig(P);
Ri = V*sqrt(L);

X = (sqrt(c) * Ri * ZZ) + mu;
XX = X.';
%%
N = 2e3;
XX = 100*rand(75,3, N);

tvec = zeros(1,numel(1:N));
for i=1:N
    tic 
    k = boundary(XX(:,:,i));
    tvec(i) = toc;
end

% XX = 100*rand(75,3, N);
pause(1);
tvec2 = zeros(1,numel(1:N));
for i=1:N
    tic 
    k = convhull(XX(:,:,i));
    tvec2(i) = toc;
end

plot(tvec); hold on; plot(tvec2);

legend({'boundary(.)', 'convhull(.)'});
xlabel('iteration'); ylabel('computation time');
