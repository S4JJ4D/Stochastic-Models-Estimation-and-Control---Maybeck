%%
% NOTE: the batch model has full col-rank H;
% Measurement Model:
% z = h^T * x + v
% where h is a a column vector, x is the column vector of unknown but
% constant parameters to be estimated. z is a scalar measurement. v is the
% scalar Gaussian random noise with specified covariance.
% each measurement is processed once it has been obtained through
% observation: Recursive estimation is used.

% Specify noise covariance in a single scalar measurement:
R = 3; Ri = inv(R);
% Inital estimate and initial covariance matrix for the parameter
xhat0 = [1;1];
V = [1 -1;1 1]; 
V = V/diag(vecnorm(V)); % normalize
P0 = V*diag([1 2])/V;
% Unknown parameter to be estimated
theta = [1.2, 3.4];
simOut = sim("filtering1.slx", 'StopTime', '15');
%%
time = simOut.tout;
xhat = simOut.logsout.getElement('x-kalman').Values.Data;
P = simOut.logsout.getElement('P-kalman').Values.Data;

close all;
subplot(2,2,1);
xplt = PlotMultivariateNormal(xhat0, P0);
hold on;
plot3(theta(1),theta(2),0, 'bo', 'MarkerFaceColor', 'b');
title('Static Parameter Estimation: Kalman Filtering');

subplot(2,2,[3 4]);
cla;
varsPlt(1) = plot(0, xhat(1,1), 'r', 'Marker', 'square', 'DisplayName', '$$\hat{x}_1$$');
hold on;
varsPlt(2) = plot(0, xhat(1,2), 'b-o', 'DisplayName', '$$\hat{x}_2$$');
grid on;
yline(theta(1),'r-','HandleVisibility','off');
yline(theta(2),'b-','HandleVisibility','off');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location','northwest');
tobj = title('Samples Processed = 0');


pause(1);
N = numel(time);
for i=1:N
    subplot(2,2,1);
    UpdatePlot(xplt, xhat(i,:), P(:,:,i));

    subplot(2,2,[3 4]);
    set(varsPlt(1), 'XData', [varsPlt(1).XData, i], 'YData', [varsPlt(1).YData, xhat(i,1)]);
    set(varsPlt(2), 'XData', [varsPlt(2).XData, i], 'YData', [varsPlt(2).YData, xhat(i,2)]);

    tobj.String = sprintf('Samples Processed = %i', i-1);

    pause(.01);
end

%%
function normPlt = PlotMultivariateNormal(mu, P)

Pi = inv(P);
% 2D normal distribution
f = @(x, y, mu, P, Pi) 1/sqrt(2*pi*det(P)) * ...
    exp(-1/2 * ( (x-mu(1)).^2 * Pi(1,1) + (y-mu(2)).^2 * Pi(2,2) + 2*(x-mu(1)).*(y-mu(2))*Pi(1,2) ));

[V,D] = eig(P);
Q = V*sqrt(D);
nn = 2;
min_int = min(nn*[Q,-Q],[],2);
max_int = max(nn*[Q,-Q],[],2);

xrange = linspace(mu(1) + min_int(1), mu(1) + max_int(1), 50);
yrange = linspace(mu(2) + min_int(2), mu(2) + max_int(2), 50);

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
nn = 3;
min_int = min(nn*[Q,-Q],[],2);
max_int = max(nn*[Q,-Q],[],2);

xrange = linspace(mu(1) + min_int(1), mu(1) + max_int(1), 50);
yrange = linspace(mu(2) + min_int(2), mu(2) + max_int(2), 50);

[X,Y] = meshgrid(xrange, yrange);

set(normPlt(1), 'XData', X, 'YData', Y, 'ZData', f(X,Y, mu, P, Pi));
set(normPlt(2), 'XData', mu(1), 'YData', mu(2));
set(normPlt(3), 'XData', [mu(1), mu(1)], 'YData', [mu(2), mu(2)],...
    'ZData', [0, f(mu(1), mu(2), mu, P, Pi)]);

end
