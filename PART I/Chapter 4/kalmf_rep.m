%% Reproducible Behavior for 'kalman' function

% xdot = -a x + g w
%    y =  h x + v
% a = 1/T T: correlation time

clear;
close all;

Ts = .02; % sampling time of the process
% The equivalent discrete-time system
a = 1;
Ad = exp(-a*Ts);
B = 0;
G = 1; % noise-input matrix
C = 1; % observation matrix
D = 0;
H = 0; 

sys = ss(Ad,[B G],C,D,Ts,'InputName',{'u' 'w'},'OutputName','y');

Q = 1; % process nosie covariance (cont-time WGN: noise power)
R = .4; % measurement noise covariance (discrete-time WGN)
Qd = Q * G^2 * 1/(2*a) * (1 - exp(-2*a*Ts));

[kalmf,L,PP,Mx,Z,My] = kalman(sys,Qd,R);

sys.InputName = {'u', 'w'};
sys.OutputName = {'yt'};
vIn = sumblk('y=yt+v');

kalmf.InputName = {'u', 'y'};

simModel = connect(sys,vIn,kalmf,{'u', 'w', 'v'},{'yt', 'y_e', 'x1_e'});

tspan = [0, 20]; % the time interval in which the measurements are taken.
t = (tspan(1):Ts:tspan(2)).'; % uniform sampling of the process
N = numel(t); % number of measurements taken

u = zeros(N,1);

rng('default');
w = sqrt(Qd)* randn(N,1);
v = sqrt(R) * randn(N,1);

x_0 = 2; % initial value for the state at time tspan(1) [Generally unkown to the observer]
xe_0 = 3;

out = lsim(simModel,[u,w,v],t,[x_0;xe_0]);

yt = out(:,1);  % true response
ye = out(:,2);  % filtered response
xe = out(:,3);  % estimated state
y = yt + v;     % measured response

figure;
tl = tiledlayout(2,1,"TileSpacing",'compact','Padding','compact');
nexttile, plot(t,yt,'b',t,ye,'k-x'), 
xlabel('Time'), ylabel('Output')
title('Kalman Filter Response')
legend('True','Filtered: Built-in kalman')
nexttile, plot(t,yt-y,'g',t,yt-ye,'r--'),
xlabel('Number of Samples'), ylabel('Error')
legend('True - measured','True - filtered')


% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

% The equivalent discrete-time system
A = -a;
Ad = exp(-a*Ts);

% the first measurement is taken at tspan(1) and the last on at tspan(2)
tvec = (tspan(1):Ts:tspan(2)).'; % uniform sampling of the process
N = numel(tvec); % number of measurements taken
sysd = ss(Ad,G,C,D,Ts);

x = lsim(sysd, w, tvec, x_0); % samples of the state process

z = C*x + v; % noisy observations of the state process x(.,.) at sample times in tspan
% A quick segue:
% mean and covariance matrix functions of the process (.,.) obey the
% following differential equations:
xhat0 = 0;
P0 = 0;
[tm, mx] = ode45(@(t,x) (A*x), linspace(tspan(1), tspan(2), 100), xhat0);
[tp, px] = ode45(@(t,p) A*p + p*A.' + G*Q*G.', linspace(tspan(1), tspan(2), 100), P0);


nexttile(1);
hold on;
% plot(tvec, x, 'm-^', 'DisplayName', 'samples of the process x(.,.)');
box on; hold on;
plot(tm, mx, 'DisplayName', 'process mean function mx(.)');

pu_2sigma = mx + 2*sqrt(px);
pl_2sigma = mx - 2*sqrt(px);
% cshade(tp, pl_2sigma, pu_2sigma);

% hold on;
% plot(tvec, z, 'DisplayName', 'Noisy Measurements of x(.,.) \cdot');
% box on;

legend;

% kalman filter loop
xhat = zeros(N,1); % corresponding to each measurement there is an estimate
P = zeros(N,1); % uncertainty in the estimate

% x_ are the best estimates of x based on all previous observations, just
% before the new measurement is taken.

tt1vec = [];
xhatvec = [];

tt2vec = [];
phatvec = [];

K = zeros(N,1);
K = Mx * ones(N,1); % use ss gain
% K = 0.25588465615;
for i=1:N
    if i==1
        x_ = xe_0;
        P_ = Z;
        
        % your knowledge of the process just before the first measurement
        % is taken. If you have information about the process in the other
        % previous time < t0, you can propagate it into the time t0 using
        % the propagation equations.
%         x_ = 0;
%         P_ = 0;
    end

    % 1. process the measurement by updating the estimate and its
    % covariance
%     K(i) = P_ * (C / (C*P_*C.' + R));
    xhat(i) = x_ + K(i)*(z(i) - C*x_);
    P(i) = P_ - K(i)*C*P_;

    % 2. propagate the estimate until the next measurement is obtained. The
    % estimate moves under the action of the vector field. Interval of
    % integration [t_{i}, t_{i+1}]
    if i~=N
    [tt1,xx] = ode45(@(t,x) x_odefun(t,x,A),    [tvec(i), tvec(i+1)], xhat(i));
    [tt2,pp] = ode45(@(t,p) p_odefun(t,p,A,Q,G),[tvec(i), tvec(i+1)], P(i));

%     tt1vec = [tt1vec; tt1];
%     xhatvec = [xhatvec; xx];
% 
%     tt2vec = [tt2vec; tt2];
%     phatvec = [phatvec; pp];

    x_ = xx(end);
    P_ = pp(end);
    end
end

nexttile(1);
hold on;
plot(tvec, xhat, 'r-o', 'MarkerSize', 4, 'DisplayName', 'Filtered: Manually Coded kalman');
pu = xhat.' + 2*sqrt(P.');
pl = xhat.' - 2*sqrt(P.');
sp2 = cshade(tvec, pl, pu);
set(sp2, 'DisplayName', '2\sigma confidence region');

% axis([-0.0535    0.7119   -0.2495    1.0385]);

% plot noisy observations:
nexttile(1);
plot(tvec, z, 'bo', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'DisplayName', 'Observations');



%%
function dxdt = x_odefun(t,x,A)
dxdt = A*x;
end

function dpdt = p_odefun(t,p,A,Q,G)
dpdt = A*p + p*A.' + G*Q*G.';
end