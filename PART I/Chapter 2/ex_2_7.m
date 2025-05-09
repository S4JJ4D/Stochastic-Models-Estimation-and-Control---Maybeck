%%
clear;
clc;

syms a1 a2 a3 real;
A = [0 -a1;a2*a3 0];
B = [a1 0;0 a3];
C = [0 1];
D = 0;

syms s
G = simplify(C*inv(s*eye(2)-A)*B + D, 'Steps', 10);
pretty(collect(G))


%%
clear;
[a1, a2, a3] = deal(realp('a1', 1), realp('a2', 1) ,realp('a3', 10));
% a2 = realp('a2', 1);
% a3 = realp('a3', 1);
A = [0 -a1;a2*a3 0];
B = [a1 0;0 a3];
C = [0 1];
D = [0 0];

sys = ss(A,B,C,D);

sys.InputName = ...
    {'\omega_d \newline Gyro Drift', 'M_{int} \newline Disturbance Torque'};
sys.OutputName = {'\omega_{im} \newline Platform Angular Velocity'};

step(sys, 10);

