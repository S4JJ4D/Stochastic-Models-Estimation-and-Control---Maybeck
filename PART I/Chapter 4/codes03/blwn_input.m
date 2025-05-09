%% Inputting Band-limited white noise instead of discretizing the system
clear;
close all;

% Studying the statistical properties of the output signal

sf_vec = [.5, .3, .2, .1, 1, 5, 10, 20, 30, 50, 70, 80];
N = numel(sf_vec);
tl = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

for sf = sf_vec
    nexttile;
    [time, EY2] = blwn_test(sf);
    plot(time, EY2, 'DisplayName', 'Emperical $\mathbf{E}[y^2]$');
    title(['sf= ', num2str(sf), ' | E[y^2]= ', num2str(mean(EY2))]);
end

title(tl, 'Target E[y^2] = 0.125');

%%
function [time, EY2] = blwn_test(sf)
[m,b,k] = deal(.01, .4, 10);

sys = ss([0,1;-k/m,-b/m], [0;1], [1/m, 0], 0);

% 1. sample the GWN by approximating it with band-limited GWN
% Specifying w0 dictates the sampling time Ts and variance of the band-limited
% white noise process
Q = 1; % cont. - time WGN
% sf: safety-factor
w0 = sf * bandwidth(sys); % band-limited white noise cut-off frequency
Ts = round(pi/w0,3,'significant');
sigma2 = Q/Ts;
mu = 0;

% 2. discretize the cont-time dynamical system and begin the simulation
sysd = c2d(sys, Ts, 'tustin');

t_max = 20;
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

% Compute Ensemble Average Power
EY2 = var(Y,0,1);

% figure;
% cla;
% plot(time, EY2, 'DisplayName', 'Emperical $\mathbf{E}[y^2]$');
% title(num2str(mean(EY2)));
end