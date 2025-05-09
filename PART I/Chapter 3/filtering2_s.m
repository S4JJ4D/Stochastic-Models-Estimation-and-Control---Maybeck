%% H is not full col-rank: can't estimate x: only some 'aspects' or 'features' of x are estimated
% non-full column rank observations give only partial information about the
% state.
%
% Measurement Model:
% z = h^T * x + v
% where h is a a column vector, x is the column vector of unknown but
% constant parameters to be estimated. z is a scalar measurement. v is the
% scalar Gaussian random noise with specified covariance.
% each measurement is processed once it has been obtained through
% observation: Recursive estimation is used.

clear;
% z \in R^2
% x \in R^3
% v \in R^2
% H \in R^(2x3)

% Specify noise covariance
V = [1 -1;1 1]; 
D = diag([3, 5]);
Q = V*sqrt(D);
V = V/diag(vecnorm(V)); % normalize
R = V*D/V;
Ri = inv(R);

% Measurement Model: model of an atomic measurement which is a two
% component vector 
H = [1 3 -2;4 4 0]; % row-rank = 2

% Inital estimate and initial covariance matrix for the parameter
xhat0 = [0;0;0];
P0 = 1e1*eye(3); % almost no initial knowledge about the constant variable;
% Unknown parameter to be estimated
theta = [8;-3;15];

%%
simOut = sim("filtering2.slx", 'StopTime', '20');
%%
time = simOut.tout;
xhat = squeeze(simOut.logsout.getElement('x-kalman').Values.Data).';
P = simOut.logsout.getElement('P-kalman').Values.Data;

clf;

tl = tiledlayout(2,3,'TileSpacing','compact', 'Padding','compact');
nexttile;
varsPlt(1) = plot(0, xhat(1,1), 'r', 'Marker', 'square', 'DisplayName', '$$\hat{x}_1$$');
hold on;
yline(theta(1),'r-','HandleVisibility','off');
grid on;
% legend('Interpreter', 'latex', 'FontSize', 14, 'Location','northwest');
title('x1');

nexttile;
varsPlt(2) = plot(0, xhat(1,2), 'b-o', 'DisplayName', '$$\hat{x}_2$$');
hold on;
yline(theta(2),'b-','HandleVisibility','off');
grid on;
title('x2');
% legend('Interpreter', 'latex', 'FontSize', 14, 'Location','northwest');

nexttile;
varsPlt(3) = plot(0, xhat(1,2), 'm-^', 'DisplayName', '$$\hat{x}_3$$');
grid on;
yline(theta(3),'m-','HandleVisibility','off');
title('x3');
% legend('Interpreter', 'latex', 'FontSize', 14, 'Location','northwest');

% --- feature extraction
xx = H * theta;
xxhat = H * xhat.';

nexttile;
varsPlt(4) = plot(0, xxhat(1,1), 'k-', 'marker', 'diamond', 'DisplayName', '$$\hat{x}_3$$');
grid on;
yline(xx(1),'m-','HandleVisibility','off');
title('xx1');

nexttile;
varsPlt(5) = plot(0, xxhat(2,1), 'k-', 'marker', 'hexagram', 'DisplayName', '$$\hat{x}_3$$');
grid on;
yline(xx(2),'m-','HandleVisibility','off');
title('xx2');

% covPlt = plot(0, xhat(1,2), 'k+');

tobj = title(tl, 'Samples Processed = 0');

pause(1);
N = numel(time);
for i=1:N

    set(varsPlt(1), 'XData', [varsPlt(1).XData, i], 'YData', [varsPlt(1).YData, xhat(i,1)]);
    set(varsPlt(2), 'XData', [varsPlt(2).XData, i], 'YData', [varsPlt(2).YData, xhat(i,2)]);
    set(varsPlt(3), 'XData', [varsPlt(3).XData, i], 'YData', [varsPlt(3).YData, xhat(i,3)]);

    set(varsPlt(4), 'XData', [varsPlt(4).XData, i], 'YData', [varsPlt(4).YData, xxhat(1,i)]);
    set(varsPlt(5), 'XData', [varsPlt(5).XData, i], 'YData', [varsPlt(5).YData, xxhat(2,i)]);

%     set(covPlt, 'XData', [covPlt.XData, i], 'YData', [covPlt.YData, norm(P(:,:,i))]);

    tobj.String = sprintf('Samples Processed = %i', i-1);

    pause(.01);
end

%% Batch Processing of Data
