%%
% xdot = -a x + g w
%    y =  h x + v
% a = 1/T T: correlation time

clear;

Ts = .02; % sampling time of the process
% The equivalent discrete-time system
a = 1;
A = exp(-a*Ts);
B = 0;
G = 1; % noise-input matrix
C = 1; % observation matrix
D = 0;
H = 0; 

x_0 = 2; % initial value for the state at time tspan(1) [Generally unkown to the observer]

xe_0 = 3; % initial estimate of the state just before the first measurement is taken.
P0 = 1; % initial estimate covariance

sys = ss(A,[B G],C,D,Ts,'InputName',{'u' 'w'},'OutputName','y');

Q = 1; % process nosie covariance (cont-time WGN: noise power)
R = .2; % measurement noise covariance (discrete-time WGN)
Qd = Q * G^2 * 1/(2*a) * (1 - exp(-2*a*Ts));

[~,~,~,K,~,~] = kalman(sys,Qd,R);