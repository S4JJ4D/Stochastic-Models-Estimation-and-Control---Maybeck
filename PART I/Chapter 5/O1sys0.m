%% First Order System: Kalman Filter: Coded
% xdot = -a x + g w
%    y =  h x + v

% a = 1/T T: correlation time

[a,g,h] = deal(1,1,1);

Ts = .02; % sampling time of the process
% The equivalent discrete-time system
A = -a;
Ad = exp(A*Ts);
H = h; % observation matrix
G = g; % noise-input matrix

Q = 1; % process nosie covariance (cont-time WGN: noise power)
R = .5; % measurement noise covariance (discrete-time WGN)
Qd = Q*g^2*1/(2*a) * (1 - exp(-2*a*Ts));


tspan = [0, 20]; % the time interval in which the measurements are taken.

% the first measurement is taken at tspan(1) and the last on at tspan(2)
tvec = (tspan(1):Ts:tspan(2)).'; % uniform sampling of the process
% tvec = (nugrid(tspan(1), tspan(2), 801)).'; %non-uniform sampling of the process

N = numel(tvec); % number of measurements taken

% generating discrete-time process noise signal 
rng('default');
w = sqrt(Qd) * randn(N,1);
v = sqrt(R)  * randn(N,1);

x0 = 3; % initial value for the state at time tspan(1) [Generally unkown to the observer]
sysd = ss(Ad,1,1,0,Ts);
%%
% Generating a particular realization of the process which is unkown to the
% observer
% The same as 'x = lsim(sysd, w, tvec, x0)', only works for non-uniform time grid too
x = [x0; zeros(N-1,1)];
for i=1:N-1
    x(i+1) = expm(A*(tvec(i+1) - tvec(i)))*x(i) + w(i);
end
%%
z = h*x + v; % noisy observations of the state process x(.,.) at sample times in tspan
% A quick segue:

% mean and covariance matrix functions of the process (.,.) obey the
% following differential equations:
xhat0 = 1.5;
P0 = 5;
% [tm, mx] = ode45(@(t,x) (A*x), linspace(tspan(1), tspan(2), 100), xhat0);
% [tp, px] = ode45(@(t,p) A*p + p*A.' + G*Q*G.', linspace(tspan(1), tspan(2), 100), P0);


tl = tiledlayout(2,1,"TileSpacing","compact","Padding","compact");
nexttile(1);
plot(tvec, x, 'DisplayName', 'samples of the process x(.,.)');
box on; hold on;
% plot(tm, mx, 'DisplayName', 'process mean function mx(.)');

% pu_2sigma = mx + 2*sqrt(px);
% pl_2sigma = mx - 2*sqrt(px);
% cshade(tp, pl_2sigma, pu_2sigma);

% hold on;
% plot(tvec, z, 'DisplayName', 'Noisy Measurements of x(.,.) \cdot');
% box on;

legend;

%% kalman filter loop
xhat = zeros(N,1); % corresponding to each measurement there is an estimate
P = zeros(N,1); % uncertainty in the estimate

% x_ are the best estimates of x based on all previous observations, just
% before the new measurement is taken.

tt1vec = [];
xhatvec = [];

tt2vec = [];
phatvec = [];

K = zeros(N,1);
for i=1:N
    if i==1
        % your knowledge of the process just before the first measurement
        % is taken. If you have information about the process in the other
        % previous time < t0, you can propagate it into the time t0(-) using
        % the propagation equations.
        x_ = xhat0; % xhat @ t = t0(-)
        P_ = P0;  % P    @ t = t0(-)
    end

    % 1. process the measurement by updating the estimate and its
    % covariance
    K(i) = P_ * (h / (h*P_*h + R));
    xhat(i) = x_ + K(i)*(z(i) - h*x_);
    P(i) = P_ - K(i)*h*P_;

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
plot(tvec, xhat, 'k-x', 'MarkerSize', 4, 'DisplayName', 'estimates of x(.,.)');
% plot(tt1vec, xhatvec, 'k-');

pu = xhat.' + 2*sqrt(P.');
pl = xhat.' - 2*sqrt(P.');
sp2 = cshade(tvec, pl, pu);
set(sp2, 'DisplayName', '95% confidence region');

dtvec = diff(tvec);
if dtvec(1) ~= dtvec(2)
    title({['Q = ', num2str(Q), ' | R = ', num2str(R)], 'Nonuniform grid'});
else
    title({['Q = ', num2str(Q), ' | R = ', num2str(R)], 'Uniform grid'});
end

nexttile(2);
plot(tvec, P, 'b-o');
plot(tt2vec, phatvec, 'k-');
box on;


% generate multiple realizations of the process:
% realCount = 2e3;
% Y = zeros(realCount, N);
% [V,D] = eig(Qd);
% T = (V*sqrt(D));
% for i=1:realCount
%     w = T*randn(2, N);
%     % verify: cov(w.') - Qd
%     Y(i,:) = lsim(sysd, w, time);
% end

set(gcf, 'Units', 'normalized', 'Position', [0.1219 0.1704 0.7068 0.6759]);
%%
function dxdt = x_odefun(t,x,A)
dxdt = A*x;
end

function dpdt = p_odefun(t,p,A,Q,G)
dpdt = A*p + p*A.' + G*Q*G.';
end