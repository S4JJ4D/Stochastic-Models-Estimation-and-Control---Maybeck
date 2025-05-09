%%
clear;
clc;

syms a1 a2 a3 real;
syms s;

A = [0 1 0;0 0 1;-a1 -a2 -a3];

adjoint(s*eye(3) - A)
det(s*eye(3) - A)

b = [0 0 1].';
syms c1 c2 c3 real;
c = [1 0 0];

CM = [b, A*b, A^2 * b]
OM = [c;c*A;c*A^2]



