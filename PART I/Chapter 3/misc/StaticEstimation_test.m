%%
% mu = [0 0];
% P = diag([1 1]);
% 
% R = mvnrnd(mu, P, 1e3);
% 
% cla;
% plot(R(:,1), R(:,2), 'b.');
% hold on;
% xline(0, 'k-');
% yline(0, 'k-');
% axis equal;


%%

% Measurement model
% z = H*x + v
x = [3;1]; % a constant deterministic value
H = [1 2;0 1]; % measuring both components through the sensing device
% H = [1 1];

% Noise Covariance
[s1,s2] = deal(1,5);
T = 1/sqrt(2) * [[1;1], [-1;1]];
R = T*diag([s1, s2])/T;

% R = 4;

x1 = [0;0];
P1 = diag([1 1]);

subplot(2,2,1);
cla;
xplt = PlotMultivariateNormal(x1, P1);
hold on;
plot3(x(1),x(2),0, 'bo', 'MarkerFaceColor', 'b');

% if exist('fg', 'var') && isvalid(fg)
%     figure(fg);
%     ax = gca;
% else
%     fg = figure;
%     ax = axes;
% end
subplot(2,2,2);
cla;
e = x1 - x;
ePlt = plot(0, vecnorm(e), 'k-o');

subplot(2,2,[3 4]);
cla;
sigma = eig(P1);
varsPlt(1) = plot(0, sigma(1), 'r-^');
hold on;
varsPlt(2) = plot(0, sigma(2), 'b-x');



pause(1);
N = 200;
for i=1:N

    % make a measurement
    z = H*x + mvnrnd([0 0], R, 1).';
%     z = H*x + randn();

    % use the measurement to update the estimate of x
    K = P1 * H.' / (H*P1*H.' + R);
    x2 = x1 + K*(z-H*x1);
    P2 = P1 - K*H*P1;

    subplot(2,2,1);
%     PlotMultivariateNormal(x2, P2);
    UpdatePlot(xplt, x2, P2);
    subplot(2,2,2);
    set(ePlt, 'XData', [ePlt.XData, i], 'YData', [ePlt.YData, vecnorm(x2-x)]);


    subplot(2,2,[3 4]);
    sigma = eig(P2);
    set(varsPlt(1), 'XData', [varsPlt(1).XData, i], 'YData', [varsPlt(1).YData, sigma(1)]);
    set(varsPlt(2), 'XData', [varsPlt(2).XData, i], 'YData', [varsPlt(2).YData, sigma(2)]);

    pause(.05);

    x1 = x2;
    P1 = P2;
end

%%
function normPlt = PlotMultivariateNormal(mu, P)

Pi = inv(P);

f = @(x, y, mu, P, Pi) 1/sqrt(2*pi*det(P)) * ...
    exp(-1/2 * ( (x-mu(1)).^2 * Pi(1,1) + (y-mu(2)).^2 * Pi(2,2) + 2*(x-mu(1)).*(y-mu(2))*Pi(1,2) ));

[V,D] = eig(P);
Q = V*sqrt(D);
nn = 2;
min_int = min(nn*Q,[],2);
max_int = max(nn*Q,[],2);

xrange = linspace(min_int(1), max_int(1), 50);
yrange = linspace(min_int(2), max_int(2), 50);

[X,Y] = meshgrid(xrange, yrange);

Z = f(X,Y, mu, P, Pi);
surfPlt = surf(X,Y,Z, 'FaceAlpha', .5, 'LineStyle', 'none');
hold on;
estimatePlt = plot3(mu(1),mu(2),0, 'r', 'Marker','square', 'MarkerFaceColor', 'r');
estimateLinePlt = plot3([mu(1), mu(1)], [mu(2), mu(2)], [0, f(mu(1), mu(2), mu, P, Pi)], 'k-');

normPlt = [surfPlt, estimatePlt, estimateLinePlt];

end

function UpdatePlot(normPlt, mu, P)

Pi = inv(P);

f = @(x, y, mu, P, Pi) 1/sqrt(2*pi*det(P)) * ...
    exp(-1/2 * ( (x-mu(1)).^2 * Pi(1,1) + (y-mu(2)).^2 * Pi(2,2) + 2*(x-mu(1)).*(y-mu(2))*Pi(1,2) ));

[V,D] = eig(P);
Q = V*sqrt(D);
nn = 2;
min_int = min(nn*Q,[],2);
max_int = max(nn*Q,[],2);

xrange = linspace(min_int(1), max_int(1), 50);
yrange = linspace(min_int(2), max_int(2), 50);

[X,Y] = meshgrid(xrange, yrange);

set(normPlt(1), 'XData', X, 'YData', Y, 'ZData', f(X,Y, mu, P, Pi));
set(normPlt(2), 'XData', mu(1), 'YData', mu(2));
set(normPlt(3), 'XData', [mu(1), mu(1)], 'YData', [mu(2), mu(2)],...
    'ZData', [0, f(mu(1), mu(2), mu, P, Pi)]);

end

