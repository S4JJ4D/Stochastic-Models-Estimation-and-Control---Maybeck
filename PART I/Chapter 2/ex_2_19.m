%%
close all;
clear;
clc;

[a,b] = deal(1,1);
tspan = [0, 10];
w0 = 1;
N = 50;
ut = linspace(tspan(1), tspan(2), N);
u = zeros(N,3);
x0 = [cos(a*w0*tspan(1)); -sin(a*w0*tspan(1)); w0];

opts_1 = odeset('RelTol',1e-8,'AbsTol',1e-9);
[t,x] = ode45(@(t,x) f(t,x, ut, u, a, b), tspan, x0, opts_1);

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1);
plot(t,x(:,1), 'DisplayName', '$x_1(t)$');
hold on;
ax2 = nexttile(2);
plot(t,x(:,2), 'DisplayName', '$x_2(t)$');
hold on;
ax3 = nexttile(3);
plot(t,x(:,3), 'DisplayName', '$x_3(t)$');
hold on;
ax4 = nexttile(4);


figure;
plot3(x(:,1), x(:,2), x(:,3), 'b-', 'LineWidth', 1.5, ...
    'DisplayName', '$x(t)$');
ax3d = gca;
axis equal;
r0 = vecnorm(x0);

hold on;
[X,Y,Z] = sphere(30);
sp_plt = surf(X*r0,Y*r0,Z*r0, 'FaceAlpha', .8, 'LineStyle', 'none');

%% Perturb the system
e0 = .2*rand(3,1);
x0_perturbed = x0 + e0;
[t1,x1] = ode45(@(t,x) f(t, x, ut, u, a, b), tspan, x0_perturbed, opts_1);
plot(ax1, t1, x1(:,1), 'DisplayName', '$\hat{x}_1(t)$');
plot(ax2, t1, x1(:,2), 'DisplayName', '$\hat{x}_2(t)$');
plot(ax3, t1, x1(:,3), 'DisplayName', '$\hat{x}_3(t)$');
plot3(ax3d, x1(:,1), x1(:,2), x1(:,3), 'k-', 'LineWidth', 1.5, ...
    'DisplayName', '\hat{x}(t)');


[t2,e] = ode45(@(t,e) fe(t,e, ut, u, a, b, w0), t, e0, opts_1);
x2 = e + x;
plot(ax1, t2, x2(:,1), 'DisplayName', '$\tilde{x}_1(t)$');
legend(ax1, 'Interpreter', 'latex', 'FontSize', 14);
plot(ax2, t2, x2(:,2), 'DisplayName', '$\tilde{x}_2(t)$');
legend(ax2, 'Interpreter', 'latex', 'FontSize', 14);
plot(ax3, t2, x2(:,3), 'DisplayName', '$\tilde{x}_3(t)$');
legend(ax3, 'Interpreter', 'latex', 'FontSize', 14);
plot3(ax3d, x2(:,1), x2(:,2), x2(:,3), 'm-', 'LineWidth', 1.5);


plot(ax4, ...
    t2, e(:,1), ...
    t2, e(:,2), ...
    t2, e(:,3));
legend(ax4, {'$e_1(t)$', '$e_2(t)$', '$e_3(t)$'}, 'Interpreter', 'latex', 'FontSize', 14);

%%
% u = [.05 * ut.^2;sqrt(ut);exp(.1*ut)].';
% plot(ut, u(:,1), '-o')
% hold on;
% plot(ut, u(:,2), '-*')
% plot(ut, u(:,3), '-^')
% qt = linspace(tspan(1), tspan(2), 50);
% vq = interp1(ut, u, qt);
% 
% plot(qt, vq(:,1), 'r*')
% hold on;
% plot(ut, u(:,2), '-*', 'MarkerFaceColor', 'auto')
% plot(ut, u(:,3), '-^', 'MarkerFaceColor', 'auto')

%%
function dx = f(t, x, ut, u, a, b)
% size(u) = [N,3]
uq = interp1(ut,u,t); % Interpolate the data set (ut,u) at time t

dx = [...
    a*x(2)*x(3) + b*uq(1);
   -a*x(1)*x(3) + b*uq(2);
    b*uq(3)];
end

function de = fe(t, e, ut, u, a, b, w0)
% size(u) = [N,3]
de = zeros(3,1);
uq = interp1(ut,u,t); % Interpolate the data set (ut,u) at time t

A = [0, a*w0, -a*sin(a*w0*t);
    -a*w0, 0, -a*cos(a*w0*t);
    0 0 0];

B = b*eye(3);
de = A*e + B*uq.';

end

