%%
clear;
clc;

syms b m k real;

A1 = [-b/m -1/m;2*k 0];
b1 = [1/m; 0];

T = [0 1;2*k 0];

A2 = simplify(T*A1*inv(T), 'Steps', 10)
b2 = T*b1
c2 = [0 1]; % solution for df/dt
d2 = 0;

syms s 
G = simplify(c2*inv(s*eye(2)-A2)*b2, 'Steps', 10);