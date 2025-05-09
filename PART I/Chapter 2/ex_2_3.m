%%
clc;
clear;
close all;
G = tf(1, [1 7 31 25]);
% system characterstic polynomial
a = G.Denominator{:};
% system order
n = numel(a) - 1;

% arbitrary state-space representation
sys = ss(G);

%% Controllable Canonical Form:
% You can simply do it by inspection. Here we adopt a systematic approach
% to perform the transformation

% See SECTION 3.9 of Elbert Hendricks, Ole Jannerup, Paul Haase Sørensen -
% Linear Systems Control_ Deterministic and Stochastic Methods-Springer (2008)
A = sys.A; B = sys.B;

P = zeros(n);
P(:,1) = B;
for i = 2:n
    P(:,i) = A*P(:,i-1) + a(i)*P(:,1);
end
P = P(:,end:-1:1);
Tcc = inv(P);
sys_cc = ss2ss(sys, Tcc);

%% Observable Canonical Form:
% You can simply do it by inspection. Here we adopt a systematic approach
% to perform the transformation

% See SECTION 3.9 of Elbert Hendricks, Ole Jannerup, Paul Haase Sørensen -
% Linear Systems Control_ Deterministic and Stochastic Methods-Springer (2008)

A = sys.A; C = sys.C;
Q = zeros(n);
Q(1,:) = C;
for i = 2:n
    Q(i,:) = Q(i-1,:)*A + a(i)*Q(1,:);
end
Q = Q(end:-1:1,:);
Toc = Q;
sys_oc = ss2ss(sys, Toc);
% do some cleanup
sys_oc.A(abs(sys_oc.A) < 1e-12) = 0;
sys_oc.C(abs(sys_oc.C) < 1e-12) = 0;

%% Observable Companion Form:
% You can simply do it by inspection and some elementary computations by
% hand. Here we adopt a systematic approach to perform the transformation.

OM = obsv(A,C);
T = OM;
sys_ocmp = ss2ss(sys, T);
% do some cleanup
sys_ocmp.A(abs(sys_ocmp.A) < 1e-13) = 0;
sys_ocmp.C(abs(sys_ocmp.C) < 1e-13) = 0;

%% Controllable Companion Form:
% You can simply do it by inspection and some elementary computations by
% hand. Here we adopt a systematic approach to perform the transformation.

CM = ctrb(A,B);
T = inv(CM);
sys_ccmp = ss2ss(sys, T);
% do some cleanup
sys_ccmp.A(abs(sys_ccmp.A) < 1e-13) = 0;
sys_ccmp.C(abs(sys_ccmp.C) < 1e-13) = 0;


%% Jordan/Modal Canonical Form:
% decomposes the system into separate (decoupled) natural modes.
syms s;
G_frac = poly2sym(G.Numerator{:}, s)/poly2sym(G.Denominator{:}, s);
G_part_frac = partfrac(G_frac)

% 1. Method 1
[V,J] = jordan(A);
T = inv(V);
sys_mod = ss2ss(sys,T);
% do some cleanup
sys_mod.A(abs(sys_mod.A) < 1e-13) = 0;
sys_mod.C(abs(sys_mod.C) < 1e-13) = 0;

% 2. Method 2
sys_mod = canon(sys, 'modal');
