%% First Order System: Kalman Filter: Using 'kalman' Function
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

sys = ss(A,[B G],C,D,Ts,'InputName',{'u' 'w'},'OutputName','y');

Q = 1; % process nosie covariance (cont-time WGN: noise power)
R = .2; % measurement noise covariance (discrete-time WGN)
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

out = lsim(simModel,[u,w,v],t,[0;0]);

yt = out(:,1);  % true response
ye = out(:,2);  % filtered response
xe = out(:,3);  % estimated state
y = yt + v;     % measured response

figure;
tl = tiledlayout(2,1,"TileSpacing",'compact','Padding','compact');
nexttile, plot(t,yt,'b',t,ye,'k-x'), 
xlabel('Time'), ylabel('Output')
title('Kalman Filter Response')
legend('True','Filtered')
nexttile, plot(t,yt-y,'g',t,yt-ye,'r--'),
xlabel('Number of Samples'), ylabel('Error')
legend('True - measured','True - filtered')

