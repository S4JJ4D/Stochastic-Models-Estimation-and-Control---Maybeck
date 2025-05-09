clc; clear; close all;
F = [0 0;1 0];
B = [0;0];
G = [0;0];
H = [0 1];
Q = 0;
R = 2;

Ts = .01;            % observation sample time
tspan = [0,30];      % time vector
x0 = [3 0].';        % deviations from nominal trajectory @ t = t_0
xhat0 = [0 0].';     % initial estimate of deviations from nominal trajectory @ t = t_0
P0 = diag([1 0]);    % confidence over the initial estimate of deviations from nominal trajectory @ t = t_0


FCN(F,B,H,G,Q,R,Ts,tspan,x0,xhat0,P0)

function [] = FCN(F,B,H,G,Q,R,Ts,tspan,x0,xhat0,P0)

% continuous-time process system model
% 'F' dynamics matrix
% 'B' control input matrix
% 'H' observation matrix
% 'G' noise input matrix:
% 'Ts' sampling time of the process for measurement purposes
% 'tspan' the time interval in which the measurements are taken: the first measurement is taken at tspan(1) and the last on at tspan(2)
% 'x0' initial value for the state at time tspan(1) [Generally unkown to the observer]
% 'xhat0' initial estimate of the state
% 'P0'  estimate certainty

% In simulation, we are also concerned with simulating the physical
% continuous-time process. To do that we should use (4-125)
% time-invariance is assumed in the following computations
Phi = expm(F*Ts);
Qd = integral(@(tau) expm(F * (Ts - tau)) * G * Q * G.' * expm(F.' * (Ts - tau)), 0, Ts, ...
    'AbsTol', 1e-12, 'RelTol', 1e-9, 'ArrayValued', true);

[Vw,Dw] = eig(Qd);
Tw = (Vw*sqrt(Dw)); % to be used in SIMULINK

[Vv,Dv] = eig(R);
Tv = (Vv*sqrt(Dv)); % % to be used in SIMULINK

[m,n] = size(H);


tvec = (tspan(1):Ts:tspan(2)).'; % uniform sampling of the process
% tvec = (nugrid(tspan(1), tspan(2), 801)).'; %non-uniform sampling of the process
N = numel(tvec); % number of measurements taken

% when discretizing, wd has dimension [n x 1] and thus the B matrix in the
% discretized system is the identity matrix.
% take both states in the output for subsequent processing. sysd is only
% used to simulate the cont-time process.
% sysd = ss(Phi, eye(n), eye(n), zeros(n), Ts);


% generating discrete-time noise process
rng('default');
w = mvnrnd(zeros(1,n), Qd, numel(tvec)).'; % process noise
v = mvnrnd(zeros(1,m), R, numel(tvec)).';  % measurement nosie

% simulate the system under the action of noise
x = [x0, zeros(n,N-1)];
for i=1:N-1
    x(:,i+1) = expm(F*(tvec(i+1) - tvec(i)))*x(:,i) + w(:,i);
end
% alternatively:
% x = lsim(sysd, w.', tvec.', x0).';

figure('Position',[252.2000 272.2000 1072 544]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

pltNames = {'x_1=a', ...
    'x_2=at',};
axlims = {[0 6.5 -3 4],[0 6. -0.19 25.6]};
for i=1:n
    ax = nexttile(i);
    plot(tvec.', x(i,:), 'b-', 'DisplayName', pltNames{i});
    hold on; box on; xlabel('time'); ylabel(pltNames{i});
    set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'Box', 'on', 'FontName', 'Source Code Pro');
    set([ax.XAxis, ax.YAxis], 'Color', '#737373', 'TickLabelColor', 'k',...
        'LineWidth',.8);
    set([ax.XLabel, ax.YLabel], 'Color', 'k');
    grid on;
    axis(axlims{i});
end

%% kalman filter loop
xhat = zeros(n,N); % corresponding to each measurement there is an estimate
P = zeros(n,n,N); % uncertainty in the estimate

z = H*x + v; % noisy observations of the state process x(.,.) at sample times in tspan

% x_ are the best estimates of x based on all previous observations, just
% before the new measurement is taken.


K = zeros(n,m,N); % filter gain
for i=1:N
    if i==1
        % your knowledge of the process just before the first measurement
        % is taken. If you have information about the process in the other
        % previous time < t0, you can propagate it into the time t0(-) using
        % the time update equations (5-36, 5-37).
        x_ = xhat0; % xhat @ t = t0(-)
        P_ = P0;  % P    @ t = t0(-)
    end
    % 1. process the measurement by updating the estimate and its
    % covariance
    K(:,:,i) = P_ * ((H.') / (H*P_*H.' + R));
    xhat(:,i) = x_ + K(:,:,i)*(z(:,i) - H*x_);
    P(:,:,i) = P_ - K(:,:,i)*H*P_;

    % 2. propagate the estimate until the next measurement is obtained. The
    % estimate moves under the action of the vector field. Interval of
    % integration [t_{i}, t_{i+1}]
    if i~=N
        x_ = Phi * xhat(:,i);
        P_ = Phi * P(:,:,i) * Phi.' + Qd;
    end
end

for i=1:n
    nexttile(i);
    % plot estimated states
    plot(tvec.', xhat(i,:), 'k-', 'DisplayName', [char([120 770]),'_',num2str(i)]);
    % plot 95 confidence region
    pu = xhat(i,:).' + 2*sqrt(squeeze(P(i,i,:)));
    pl = xhat(i,:).' - 2*sqrt(squeeze(P(i,i,:)));
    sp = cshade(tvec, pl, pu);
    set(sp, 'DisplayName', '95% confidence region');
    legend('FontSize',10,'Location','southeast');
end

title(tl, ['Q=', mat2str(Q), ' | R=', mat2str(R), ' | Ts=', num2str(Ts), ...
    newline, ['H=',mat2str(H)]], ...
    'FontName', 'Source Code Pro');

%% Plot Covariance
Pcell = squeeze(mat2cell(P,n,n,ones(1,N)));
Peigs = cellfun(@eig,Pcell,'UniformOutput',false);
PeigsNorm = cellfun(@vecnorm,Peigs); %1

Peigs_ = Peigs.';
PeigsMat = [Peigs_{:}];  %2


% figure('Position',[313 57 984.8000 803.2000]);
% ax = axes;
% plot(tvec, PeigsNorm, 'b-', 'DisplayName', '|\sigma(P)|');
% hold on; box on; xlabel('time'); ylabel('|\sigma(P)|');
% set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
%     'Box', 'on', 'FontName', 'Source Code Pro');
% set([ax.XAxis, ax.YAxis], 'Color', '#737373', 'TickLabelColor', 'k',...
%     'LineWidth',.8);
% set([ax.XLabel, ax.YLabel], 'Color', 'k');
% grid on;
% title(['Estimation Error Covariance Norm |\sigma(P)| Over Time', ...
%     newline, ['H=',mat2str(H)]]);
% legend;


% figure('Position',[313 57 984.8000 803.2000]);
% tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% for i=1:n
%     ax = nexttile(i);
%     plot(tvec.', PeigsMat(i,:), 'b-', 'DisplayName', ['\sigma_', num2str(i)]);
%     hold on; box on; xlabel('time'); ylabel(['\sigma_', num2str(i)]);
%     set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
%         'Box', 'on', 'FontName', 'Source Code Pro');
%     set([ax.XAxis, ax.YAxis], 'Color', '#737373', 'TickLabelColor', 'k',...
%         'LineWidth',.8);
%     set([ax.XLabel, ax.YLabel], 'Color', 'k');
%     grid on;
% end
% title(tl, ['Case: H=' , mat2str(H)]);


figure('Position',[252.2000 272.2000 1072 544]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for i=1:n
    ax = nexttile(i);
    plot(tvec, squeeze(P(i,i,:)), 'b-', 'DisplayName', sprintf('P_{%d%d}',i,i));
    hold on; box on; xlabel('time'); ylabel(sprintf('P_{%d%d}',i,i));
    set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'Box', 'on', 'FontName', 'Source Code Pro');
    set([ax.XAxis, ax.YAxis], 'Color', '#737373', 'TickLabelColor', 'k',...
        'LineWidth',.8);
    set([ax.XLabel, ax.YLabel], 'Color', 'k');
    grid on;
end
title(tl, ['Case: H=', mat2str(H)], 'FontName', 'Source Code Pro');

end
