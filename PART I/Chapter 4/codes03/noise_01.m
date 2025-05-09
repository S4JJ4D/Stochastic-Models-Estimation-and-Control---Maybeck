% Simulation of EXAMPLE 4.8 
clear;
a = 2; % a = 1/T
K = 2;

%xdot = -a*x + K/T * w
sys = ss(-a,K*a,1,0);

Ts = .01;
Ad = exp(-a * Ts);
Bd = 1;
Cd = 1;
Dd = 0;

sysd = ss(Ad,Bd,Cd,Dd,Ts);

%%
Q = 1;
Qd = Q * K^2 * 1/2 * a * (1 - exp(-2*a*Ts));

%%
t_max = 30;
time = 0:Ts:t_max;
N = numel(time);

% i = @(tau) expm(sys.A * (Ts - tau)) * sys.B * (sys.B).' * expm((sys.A).' * (Ts - tau));
% propagation of covariance
Px = zeros(1,N);
for i=1:N-1
    Px(i+1) = expm(sys.A * Ts) * Px(i) * expm(sys.A * Ts).' + Qd;
end
%%
% generate multiple realizations of the process:
realCount = 3e3;
Y = zeros(realCount, N);
for i=1:realCount
    w = sqrt(Qd).*randn(1, N);
    Y(i,:) = lsim(sysd, w, time);
end

%% Plot some of the realizations
cla;
idxs = randperm(realCount, 3);
for idx=idxs
plot(time, Y(idx, :));
hold on;
end

%% Compute Ensemble Average Power
EY2 = var(Y,0,1);
figure;
cla;
plot(time, EY2, 'DisplayName', 'Emperical $\mathbf{E}[y^2]$');
hold on;
plot(time, Px, 'DisplayName', 'Theoretical $\mathbf{E}[y^2]$');
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'southeast');
title(num2str(mean(EY2)));
