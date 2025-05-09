clear;
close all;

P = [1, 1/2;
    1/2, 1];
mu = [0;0];

X = cellipse(mu, P, 2);
patch('XData',X(1,:),'YData',X(2,:), 'FaceColor', 'r', 'FaceAlpha', .2, ...
    'DisplayName', '$2\sigma$ - ellipse');
axis equal; hold on; box on; 
xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 13);

ax = gca;
set(ax, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', ...
    'XTickLabel', {}, 'YTickLabel', {});

legend('Interpreter', 'latex');
%%
figure;

P = [0, 0;
    0, 3/4];
mu = [1;1/2];

X = cellipse(mu, P, 2);
patch('XData',X(1,:),'YData',X(2,:), 'FaceColor', 'r', 'FaceAlpha', .2, ...
    'DisplayName', '$2\sigma$ - ellipse');
axis equal; hold on; box on; 
xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 13);

ax = gca;
set(ax, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin', ...
    'XTickLabel', {}, 'YTickLabel', {});


%%

plot(mu(1),mu(2), 'kx', 'HandleVisibility', 'off');

f = @(z,sigma,mu) 1/(sigma * sqrt(2*pi)) * exp(-1/2 * ((z - mu)/sigma).^2);

sigmaZ = 1;
muZ = 0;
z = linspace(-2*sigmaZ, 2*sigmaZ, 100);
fz = f(z,sigmaZ, muZ);

X = [z;fz;zeros(1,100);ones(1,100)];

T = makehgtform('translate', [mu;0]) * makehgtform('zrotate', pi/2);
Xt = T*X;

% plt = area(z,f(z,sigmaZ, muZ), 'FaceAlpha', .1, 'DisplayName', '$f_{x_2|x_1=z}$');

plt = plot(Xt(1,:),Xt(2,:), 'DisplayName', '$f_{x_2|x_1=z}$');

hold on;

% h = hgtransform;
% plt.Parent = h;
% h.Matrix = T;

legend('Interpreter', 'latex', 'FontSize', 12);





