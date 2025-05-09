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
% set([ax.XAxis, ax.YAxis], 'Color', '#737373', 'TickLabelColor', 'k',...
%    'LineWidth',.8);
% set([ax.XLabel, ax.YLabel], 'Color', 'k')
% --------------------------------------------------------

%%
clc;
clear;
close all;
G = tf([.3 .3*.01], [1 .006 .003]);
% system characterstic polynomial
a = G.Denominator{:};
% system order
n = numel(a) - 1;

% arbitrary state-space representation
sys = ss(G);

%% Controllable Canonical Form:
% You can simply do it by inspection. Here we adopt a systematic approach
% to perform the transformation

% See SECTION 3.9 of Elbert Hendricks, Ole Jannerup, Paul Haase Sørensen -
% Linear Systems Control_ Deterministic and Stochastic Methods-Springer (2008)
A = sys.A; B = sys.B;

P = zeros(n);
P(:,1) = B;
for i = 2:n
    P(:,i) = A*P(:,i-1) + a(i)*P(:,1);
end
P = P(:,end:-1:1);
Tcc = inv(P);
sys_cc = ss2ss(sys, Tcc)


damp(sys_cc)
stepinfo(sys_cc)

step(sys_cc)
%%
F = [sys_cc.A, sys_cc.B;0 0 0];
G = [sys_cc.B;0];
B = [0 0 0].';
H = [sys_cc.C,0];
Q = 400;
R = 900;

% for information purposes
% Fd = expm(F);
% Qd = integral(@(tau) expm(F * (Ts - tau)) * G * Q * G.' * expm(F.' * (Ts - tau)), 0, Ts, ...
%     'AbsTol', 1e-12, 'RelTol', 1e-9, 'ArrayValued', true);

% true value of the initial altitude (unknown to observer) 6000
% true value of the commanded altitude (unknown to observer) 12000
% initial estimate of the altitude: 2000 - std: 1000ft -> variance: 1e6 ft^2
% initial estimate of the commanded altitude: 10000 - std: 500ft -> variance: 250000 ft^2

[Ts,tspan,x0,xhat0,P0] = deal(1,[0, 1e3],...
    [10000;20000;12000],[0;6666.7;10000],blkdiag(blkdiag(0,1.111e+07),250000));
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

x = [x0, zeros(n,N-1)];
for i=1:N-1
    x(:,i+1) = expm(F*(tvec(i+1) - tvec(i)))*x(:,i) + w(:,i);
end
% alternatively:
% x = lsim(sysd, w.', tvec.', x0).';


figure;
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(1);
plot(tvec.', x(1,:), 'DisplayName', 'x_1');
hold on; box on; xlabel('time'); ylabel('x_1');
ax.FontName = 'Source Code Pro';
grid on;

ax = nexttile(2);
plot(tvec.', x(2,:), 'b-', 'DisplayName', 'x_2');
hold on; box on; xlabel('time'); ylabel('x_2');
ax.FontName = 'Source Code Pro';
grid on;

ax = nexttile(3);
plot(tvec.', x(3,:), 'b-', 'DisplayName', 'x_3: h_c: Commanded Altitude');
hold on; box on; xlabel('time'); ylabel('x_3: h_c');
set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Box', 'on', 'FontName', 'Source Code Pro');
set([ax.XAxis, ax.YAxis], 'Color', '#737373', 'TickLabelColor', 'k',...
    'LineWidth',.8);
set([ax.XLabel, ax.YLabel], 'Color', 'k');
grid on;
% axis([0	26	9000	12500]);

ax = nexttile(4);
plot(tvec.', H*x, 'b-', 'DisplayName', 'h: Altitude');
hold on; box on; xlabel('time'); ylabel('h: Altitude');
ax.FontName = 'Source Code Pro';
% axis([964         968       13909       14779]);
grid on;

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


nexttile(1);
plot(tvec.', xhat(1,:), 'k-x', 'DisplayName', 'xhat_1');
pu = xhat(1,:).' + 2*sqrt(squeeze(P(1,1,:)));
pl = xhat(1,:).' - 2*sqrt(squeeze(P(1,1,:)));
sp1 = cshade(tvec, pl, pu);
set(sp1, 'DisplayName', '95% confidence region');
legend;

nexttile(2);
plot(tvec.', xhat(2,:), 'k-x', 'DisplayName', 'xhat_2');
pu = xhat(2,:).' + 2*sqrt(squeeze(P(2,2,:)));
pl = xhat(2,:).' - 2*sqrt(squeeze(P(2,2,:)));
sp2 = cshade(tvec, pl, pu);
set(sp2, 'DisplayName', '95% confidence region');
legend;

nexttile(3);
plot(tvec.', xhat(3,:), 'k-x', 'DisplayName', 'xhat_3');
pu = xhat(3,:).' + 2*sqrt(squeeze(P(3,3,:)));
pl = xhat(3,:).' - 2*sqrt(squeeze(P(3,3,:)));
sp3 = cshade(tvec, pl, pu);
set(sp3, 'DisplayName', '95% confidence region');
legend;

nexttile(4);
hhat = H*xhat;
Ph = squeeze(pagemtimes(pagemtimes(H,P),H.'));
plot(tvec.', hhat, 'k-x', 'DisplayName', 'h-hat');
hold on;
plot(tvec.', z, 'wo', 'MarkerEdgeColor', 'k', 'DisplayName', 'z');
pu = hhat.' + 2*sqrt(Ph);
pl = hhat.' - 2*sqrt(Ph);
sp4 = cshade(tvec, pl, pu);
set(sp4, 'DisplayName', '95% confidence region');
legend;


title(tl, ['Q=', num2str(Q), ' | R=', num2str(R), ' | Ts=', num2str(Ts)], ...
    'FontName', 'Source Code Pro');

    
figure;
tl = tiledlayout(2,2,"TileSpacing","compact","Padding","compact");

ax = nexttile(1);
plot(tvec, squeeze(P(1,1,:)), 'k-o', 'DisplayName', '\sigma^2_{e1}');
box on; xlabel('time'); ylabel('\sigma^2_{e1}');
ax.FontName = 'Source Code Pro';
legend;

ax = nexttile(2);
plot(tvec, squeeze(P(2,2,:)), 'k-o', 'DisplayName', '\sigma^2_{e2}');
box on; xlabel('time'); ylabel('\sigma^2_{e2}');
ax.FontName = 'Source Code Pro';
legend;

ax = nexttile(3);
plot(tvec, squeeze(P(3,3,:)), 'k-o', 'DisplayName', '\sigma^2_{e3}');
box on; xlabel('time'); ylabel('\sigma^2_{e3}');
title({'Error Variance in', 'Commanded Altitude Estimate'});
ax.FontName = 'Source Code Pro';
legend;

ax = nexttile(4);
plot(tvec, Ph, 'k-o', 'DisplayName', '\sigma^2_{eh}');
box on; xlabel('time'); ylabel('\sigma^2_{eh}');
title({'Error Variance in', 'Altitude Estimate'});
% title('Error Variance in Altitude Estimate');
ax.FontName = 'Source Code Pro';
legend;

Kend = sprintf('[%.3f,%.3f,%.3f]',K(1,1,end),K(2,1,end),K(3,1,end));
title(tl, ['K(\infty) = ', Kend], 'FontName', 'Source Code Pro');

end


