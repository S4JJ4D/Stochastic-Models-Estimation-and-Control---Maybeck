% (I): NO control inputs are applied. System is driven only by WGN.
% continous-time model (model of the real physical process) is given by
% \dot{x} = A x + G w
% where w is WGN with intensity Q (possibly time-varying)
% In order to "simulate the system", we construct the discrete-time
% "equivalent model" in the sense that the discrete-time model’s 
% values of x(t1), x(t2), . . . , are identical (in probabilistic sense)
% to those of the continuous-time model at these particular times. 
% this discrete-time model is given by
% x(k+1) = Ad x(k) + wd
% where Ad = expm(A * Ts) and Ts is the sampling period
% and wd is the discrete-time GWN with covariance matrix Qd given by
% Equation (4-127b). 

close all;
clear;

[Q,R] = deal(2,0);
FCN(Q,R);

pause(1);

[Q,R] = deal(2,0.2);
FCN(Q,R);


function [] = FCN(Q,R)


[r1,r2,c1,c2] = deal(1,1,1,1);
% continuous-time process system model
F = [(-1/c1 * (1/r1 + 1/r2)), 1/(c1*r2);(1/(c2*r2)), (-1/(c2*r2))]; % dynamics matrix
B = [1/(c1*r1);0]; % control input matrix
H = [0 1]; % observation matrix
G = [1;0]; % noise input matrix:

Ts = .5; % sampling time of the process for measurement purposes

% In simulation, we are also concerned with simulating the physical
% continuous-time process. To do that we should use (4-125)
Phi = expm(F*Ts);
Qd = integral(@(tau) expm(F * (Ts - tau)) * G * Q * G.' * expm(F.' * (Ts - tau)), 0, Ts, ...
    'AbsTol', 1e-12, 'RelTol', 1e-9, 'ArrayValued', true);

[Vw,Dw] = eig(Qd);
Tw = (Vw*sqrt(Dw));

[Vv,Dv] = eig(R);
Tv = (Vv*sqrt(Dv));

[m,n] = size(H);


tspan = [0, 30]; % the time interval in which the measurements are taken.
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
x0 = [0;0]; % initial value for the state at time tspan(1) [Generally unkown to the observer]

x = [x0, zeros(n,N-1)];
for i=1:N-1
    x(:,i+1) = expm(F*(tvec(i+1) - tvec(i)))*x(:,i) + w(:,i);
end
% alternatively:
% x = lsim(sysd, w.', tvec.', x0).';


figure;
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile(1);
plot(tvec.', x(1,:), 'DisplayName', 'x1');
hold on; box on; xlabel('time'); ylabel('x1');
nexttile(2);
plot(tvec.', x(2,:), 'b-o', 'DisplayName', 'x2');
hold on; box on; xlabel('time'); ylabel('x2');

%% kalman filter loop
xhat = zeros(2,N); % corresponding to each measurement there is an estimate
P = zeros(2,2,N); % uncertainty in the estimate

z = H*x + v; % noisy observations of the state process x(.,.) at sample times in tspan
xhat0 = [0;0]; % initial estimate of the state
P0 = 1*eye(2);   % estimate certainty

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