%%
clear;
clc;

%% Transfer Function

syms k real;
syms s;

g1(s) = k/(s+6);
G1(s) = simplify(g1/(1+g1));
G2(s) = 1/(s+6);
ff = simplify(G1(s)*G2(s));
% closed-loop transfer function
H = collect(simplify(ff/(1+3*ff)),s)

[n,d] = numden(H)
sol = solve(d,s)
syms k
limit(sol(1),k,inf)
limit(sol(2),k,inf)

%% state-space model;
A = [-6 1;-3*k -6-k];
b = [0;k];
c = [1 0];
d = 0;

%% Root locus
k_range = -10:.01:10;
N = numel(k_range);
sys_poles = zeros(N,1);

close all;
figure;
plot(0,0, 'ko', 'MarkerFaceColor', 'k');
hold on;
xline(0,'k-');
yline(0, 'k-');

to = title('k=');

for k=-4:.1:20
    r = roots([1 (k+12) (9*k+36)]);
    if isreal(r)
        plot(r,[0;0], 'bo');
    else
        plot(r,'bo');
    end
    to.String = ['k=', num2str(k, 4)];
    pause(.02);
end

%%
k = realp('k',1);
g1 = tf(k,[1 6]);
ff = feedback(g1,1)*tf(1,[1 6]);
sys = feedback(ff, 3);


figure;
plot(0,0, 'ko', 'MarkerFaceColor', 'k');
hold on;
xline(0,'k-');
yline(0, 'k-');


for kval=1:1:10
    sys.Blocks.k.Value = kval;
    step(sys);
end

% clear A b;
% A = @(k) [-6 1;-3*k -6-k];
% b = @(k) [0;k];
% c = [1 0];
% d = 0;
% 
% syst= ss(A(0), b(1), c, d);

% for kval=1:1:20
%     syst.A = A(kval);
%     syst.B = b(kval);
%     step(syst, 2);
% end
% 

%% Analytical solution for the roots
[~, den] = numden(H);
sol = solve(den, s);
% if 0 <= k <= 12, we have complex roots, else roots are real








