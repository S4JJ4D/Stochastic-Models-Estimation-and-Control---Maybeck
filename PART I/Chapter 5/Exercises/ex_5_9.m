% ------------- MATHEMATICA STYLE PLOTTING ---------------
% Mathematica Style:
% 'FontName', 'Source Code Pro'
% Area Colors: 'FaceColor', '#DAD9EB', 'FaceColor', '#F6E1BE'
% Edge Colors: 'EdgeColor', '#4A457F', 'EdgeColor', '#E5A73C'
% Axes Style:
% set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
%     'XAxisLocation', 'origin', 'YAxisLocation', 'origin', ...
%     'Box', 'on', ...
%     'FontName', 'Source Code Pro');
% --------------------------------------------------------
clear;
close all;
FCN(.02,10);


% FCN(2*pi,1000);

function [] = FCN(Ts, TF)

% continuous-time process system model
F = [0 1;-1 0]; % dynamics matrix
B = 0; % control input matrix
H = [1 0]; % observation matrix
G = 0; % noise input matrix:

Q = 0;
R = 1;

% Ts = .05; % sampling time of the process for measurement purposes

% In simulation, we are also concerned with simulating the physical
% continuous-time process. To do that we should use (4-125)
Phi = expm(F*Ts);
Qd = integral(@(tau) expm(F * tau) * G * Q * G.' * expm(F.' * tau), 0, Ts, ...
    'AbsTol', 1e-12, 'RelTol', 1e-9, 'ArrayValued', true);

[Vw,Dw] = eig(Qd);
Tw = (Vw*sqrt(Dw));

[Vv,Dv] = eig(R);
Tv = (Vv*sqrt(Dv));

[m,n] = size(H);


tspan = [0, TF]; % the time interval in which the measurements are taken.
% the first measurement is taken at tspan(1) and the last on at tspan(2)
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
x0 = [3;2]; % initial value for the state at time tspan(1) [Generally unkown to the observer]

x = [x0, zeros(n,N-1)];
for i=1:N-1
    x(:,i+1) = expm(F*(tvec(i+1) - tvec(i)))*x(:,i) + w(:,i);
end
% alternatively:
% x = lsim(sysd, w.', tvec.', x0).';


figure;
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile(1);
plot(tvec.', x(1,:), 'b-','DisplayName', '$x_1(t)$');
hold on; box on; xlabel('time'); ylabel('x1');
nexttile(2);
plot(tvec.', x(2,:), 'b-', 'DisplayName', '$x_2(t)$');
hold on; box on; xlabel('time'); ylabel('x2');
%% State-Space Plot
if false
    figure;
    plot(x(1,:), x(2,:), 'k-', 'DisplayName', 'x(t)');
    hold on; axis equal; box on;
    
    mu0 = [1;2]; % initial estimate of the state
    Pxx0 = [4 1;1 2];   % estimate certainty
    
    mx = [mu0, zeros(2,N-1)];
    Pxx = zeros(2,2,N);
    Pxx(:,:,1) = Pxx0;
    
    X = cellipse(mu0, Pxx0, 2);
    cpatch = patch('XData',X(1,:),'YData',X(2,:), 'FaceColor', 'r', 'FaceAlpha', .2);  
    mplt = plot(mu0(1), mu0(2), 'k-');
    
    for i=1:N-1
        mx(:,i+1) = Phi * mx(:,i);
        Pxx(:,:,i+1) = Phi * Pxx(:,:,i) * Phi.';
        
        X = cellipse(mx(:,i+1), Pxx(:,:,i+1), 2);
        set(cpatch,'XData',X(1,:),'YData',X(2,:));
        set(mplt, 'XData', [mplt.XData, mx(1,i+1)], ...
            'YData', [mplt.YData, mx(2,i+1)]);
    
        pause(.1);
    end
end

%% kalman filter loop
xhat = zeros(2,N); % corresponding to each measurement there is an estimate
P = zeros(2,2,N); % uncertainty in the estimate

z = H*x + v; % noisy observations of the state process x(.,.) at sample times in tspan
xhat0 = [1;1]; % initial estimate of the state
P0 = 1*eye(2);   % estimate certainty
P0 = [4 1;1 2];   % estimate certainty

% x_ are the best estimates of x based on all previous observations, just
% before the new measurement is taken.

figure;
plot(x(1,:), x(2,:), 'k-', 'LineWidth', 1, 'DisplayName', 'x(t)');
hold on; axis equal; box on;
ssax = gca;
X = cellipse(xhat0, P0, 2);
cpatch = patch('XData',X(1,:),'YData',X(2,:), 'FaceColor', 'r', 'FaceAlpha', .2);  
mplt = plot(xhat0(1), xhat0(2), 'k-o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b');
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

    X = cellipse(xhat(:,i), P(:,:,i), 2);
    set(cpatch,'XData',X(1,:),'YData',X(2,:));
    set(mplt, 'XData', [mplt.XData, xhat(1,i)], ...
        'YData', [mplt.YData, xhat(2,i)]);
    
    pause(.02);

    % 2. propagate the estimate until the next measurement is obtained. The
    % estimate moves under the action of the vector field. Interval of
    % integration [t_{i}, t_{i+1}]
    if i~=N
        x_ = Phi * xhat(:,i);
        P_ = Phi * P(:,:,i) * Phi.' + Qd;
    end
end

set(ssax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Box', 'on', ...
    'FontName', 'Source Code Pro');
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 13);
title('State Space', 'Interpreter', 'latex', 'FontSize', 14);
xline(0, 'k-', 'HandleVisibility', 'off');
yline(0, 'k-', 'HandleVisibility', 'off');
ssax.XTickLabel = {}; ssax.YTickLabel = {};
legend('$x(t)$','$86$\% Confidence Region','$\widehat{x}(t_i)$', 'Interpreter', 'latex', 'FontSize', 11);

nexttile(tl,1);
plot(tvec.', xhat(1,:), 'k-', 'DisplayName', '$\widehat{x}_1(t_i)$');
pu = xhat(1,:).' + 2*sqrt(squeeze(P(1,1,:)));
pl = xhat(1,:).' - 2*sqrt(squeeze(P(1,1,:)));
sp1 = cshade(tvec, pl, pu);
set(sp1, 'DisplayName', '$95$\% Confidence Region');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Box', 'on', ...
    'FontName', 'Source Code Pro');
xlabel('time', 'FontSize', 10);
ylabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
legend('Interpreter', 'latex', 'FontSize', 12);

nexttile(tl,2);
plot(tvec.', xhat(2,:), 'k-', 'DisplayName', '$\widehat{x}_2(t_i)$');
pu = xhat(2,:).' + 2*sqrt(squeeze(P(2,2,:)));
pl = xhat(2,:).' - 2*sqrt(squeeze(P(2,2,:)));
sp2 = cshade(tvec, pl, pu);
set(sp2, 'DisplayName', '$95$\% Confidence Region');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Box', 'on', ...
    'FontName', 'Source Code Pro');
xlabel('time', 'FontSize', 10);
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
legend('Interpreter', 'latex', 'FontSize', 12);


title(tl, ['Q = ', num2str(Q), ' | R = ', num2str(R), ' | Ts = ', num2str(Ts)]);


figure;
tl = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

nexttile;
plot(tvec, squeeze(P(1,1,:)), 'k-o', 'DisplayName', '\sigma^2_{e1}');
box on; xlabel('time'); ylabel('\sigma^2_{e1}');
legend;

nexttile;
plot(tvec, squeeze(P(2,2,:)), 'k-o', 'DisplayName', '\sigma^2_{e2}');
box on; xlabel('time'); ylabel('\sigma^2_{e2}');
legend;

Kend = sprintf('[%.3f,%.3f]',K(1,1,end),K(2,1,end));
title(tl, ['K(\infty) = ', Kend]);

end