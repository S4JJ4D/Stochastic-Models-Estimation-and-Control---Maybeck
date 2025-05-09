%%
[a,b] = deal(.5,1);
% 1/a = T which is the correlation time of the process
x0 = 2;

A = -a;
B = b;
C = 1;
D = 0;

Ts = .1; % measurement period
R = .5; % measurement noise covariance

T_init = 1; % after T_init, measurement of the process begins

