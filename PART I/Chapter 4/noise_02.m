% simulation of mass-spring-damper system
clear;
[m,b,k] = deal(.01, .4, 10);

sys = ss([0,1;-k/m,-b/m], [0;1], [1/m, 0], 0);

Ts = .01;
Ad = expm(sys.A * Ts);
Bd = eye(2);
Cd = sys.C;
Dd = zeros(1,2);

sysd = ss(Ad,Bd,Cd,Dd,Ts);
%%
% single noise is input to the system
G = [0;1];
Q = eye(1);
f = @(tau) expm(sys.A * (Ts - tau)) * G * Q * G.' * expm((sys.A).' * (Ts - tau));
Qd = integral(f, 0, Ts, 'ArrayValued', true);

%%
t_max = 20;
time = 0:Ts:t_max;
N = numel(time);

% propagation of covariance
Px = zeros(2,2,N);
for i=1:N-1
    Px(:,:,i+1) = expm(sys.A * Ts) * Px(:,:,i) * expm(sys.A * Ts).' + Qd;
end
Py = squeeze(pagemtimes(Cd,pagemtimes(Px, Cd.')));
%%
% generate multiple realizations of the process:
realCount = 2e3;
Y = zeros(realCount, N);
[V,D] = eig(Qd);
T = (V*sqrt(D));
for i=1:realCount
    w = T*randn(2, N);
    % verify: cov(w.') - Qd
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
plot(time, Py, 'DisplayName', 'Theoretical $\mathbf{E}[y^2]$');
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'southeast');
title(num2str(mean(EY2)));
