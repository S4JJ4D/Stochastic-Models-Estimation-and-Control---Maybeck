% simulation of mass-spring-damper system
clear;
[m,b,k] = deal(.01, .4, 10);

sys = ss([0,1;-k/m,-b/m], [0;1], [1/m, 0], 0);

% 1. sample the GWN by approximating it with band-limited GWN
% Specifying w0 dictates the sampling time Ts and variance of the band-limited
% white noise process
Q = 1; % cont. - time WGN
w0 = 50 * bandwidth(sys); % band-limited white noise cut-off frequency
Ts = round(pi/w0,3,'significant');
sigma2 = Q/Ts;
mu = 0;

% 2. discretize the cont-time dynamical system and begin the simulation
sysd = c2d(sys, Ts, 'tustin');

t_max = 10;
time = 0:Ts:t_max;
N = numel(time);

% generate multiple realizations of the process:
realCount = 2e3;
Y = zeros(realCount, N);
for i=1:realCount
    w = sqrt(sigma2).*randn(1, N) + mu;
    % verify: cov(w.') - Qd
    Y(i,:) = lsim(sysd, w, time);
end

%% Plot some of the realizations
% cla;
% idxs = randperm(realCount, 3);
% for idx=idxs
% plot(time, Y(idx, :));
% hold on;
% end

%% Compute Ensemble Average Power
EY2 = var(Y,0,1);
figure;
cla;
plot(time, EY2, 'DisplayName', 'Emperical $\mathbf{E}[y^2]$');
title(num2str(mean(EY2)));
% hold on;
% plot(time, Py, 'DisplayName', 'Theoretical $\mathbf{E}[y^2]$');
% legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'southeast');
