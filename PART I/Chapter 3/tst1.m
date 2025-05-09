%% Full Column Rank Measurements
% When the atomic measurement is full-column rank, more and more
% measurements of the constant variable (through the same device; constant
% H), yields better and better estimates:
% 1. Full col-rank atomic observations yield full col-rank batch observations:
% H:full col-rank --> [H;H;...;H]: full col-rank
% 2. And thus LS has a unique solution of full col-rank observations (atomic/batch) 
%
% When H is full col-rank, all the the components of x are "captured" by
% the observation vector. In this case H can be square or tall matrix
%
% Measurement Model:
% z = Hx + v
% where H is a full col-rank matrix, x is the column vector of unknown but
% constant parameters to be estimated. z is an atomic vector measurement. 
% v is the Gaussian random noise with specified covariance.
% each measurement is processed once it has been obtained through
% observation: Recursive estimation is used.
% 
% Note that in this case, we seems to be rare in practice, m > n: the size
% of an atomic measurement is larger than the number of parameters.
% Usually, the size of the measurement vector z is smaller than the size of
% x.

clear;

% z \in R^m
% x \in R^n
% v \in R^m
% H \in R^(mxn)

% Measurement Model: model of an atomic measurement
% Note that the measurement model has constant H: we make repeated
% measurements through identical device/model.
H = [1 2;-1 3;0 4]; % full col-rank observation
H = [1 3;4 0]; % row-rank=2=col-count -> full col-rank --> entire state is captured.
% H = [1 0;2 0;3 0]; % row-rank=1 --> 1 feature is estimated.
H = [1 1];


xhat0 = [-5;-5]; % Inital estimate and initial covariance matrix for the parameter

n = size(xhat0, 1); % number of parameters to be estimated
m = size(H,1); % number of measurement variables taken in a single atomic measurement

V = randi(10,n,n);
[V,~] = qr(V); % orthonormalize V
D = diag(randi(5,1,n)); % noise variances in those directions
P0 = V*sqrt(D)*V.'; % almost no initial knowledge about the constant variable;
theta = [8;-3]; % Unknown parameter to be estimated

% Specify noise covariance
% V = [[1;1;1],[1;-1;1],[-1;-1;1]];
% V = [1 -1;1 1];
% V = V/diag(vecnorm(V)); % normalize
V = randi(10,m,m);
[V,~] = qr(V); % orthonormalize V
D = diag(randi(5,1,m)); % noise variances in those directions
Q = V*sqrt(D);
R = V*D/V;
Ri = inv(R);

%%
simOut = sim("estimation_mdl.slx", 'StopTime', '10');
%%
time = simOut.tout;
xhat = squeeze(simOut.logsout.getElement('x-kalman').Values.Data);
z = squeeze(simOut.logsout.getElement('z').Values.Data);
P = simOut.logsout.getElement('P-kalman').Values.Data;
%%
close all;
fg = figure('Units', 'normalized');
fg.Position = [0.0083    0.0704    0.3031    0.8111];
tl = tiledlayout(2,1,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
n = numel(xhat0);
for i=1:n
    plot(time, xhat(i,:));
    hold on;
end
grid on; box on; xlabel('time'); ylabel('x_{hat}');
xl = xlim;
plot(xl, [theta(1), theta(1)], 'k-');
plot(xl, [theta(2), theta(2)], 'k-');

yl = ylim;
vline = plot([0,0],yl, 'k-');
%%
nexttile;
plot(theta(1), theta(2), 'bo', 'marker', 'square', 'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 8, ...
    'DisplayName', 'Constant Parameter');
hold on; grid on; box on; xlabel('x1'); ylabel('x2');
tobj = title({'State-Space', 'Processed Measurements = 0'});


X = cellipse(xhat0, P0, 2); % c=2 contains 0.86466 probability
% X = cellipse(mx(tvec(1)), P(tvec(1)), 1);
cellipse_plt = patch('XData',X(1,:),'YData',X(2,:),...
    'FaceColor', 'r', 'FaceAlpha', .1, ...
    'DisplayName', '2-cellipise');

estimate_plt = plot(xhat0(1), xhat0(2), 'bo', 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'b', 'MarkerSize', 8, ...
    'DisplayName', 'Estimate');

if all(size(H) == [1, 2])
    a = 1;
    endPoints = [-a*H.', a*H.'];
    plot(endPoints(1,:), endPoints(2,:), 'HandleVisibility', 'off');
end
plot(0,0, 'k*', 'MarkerSize', 3, 'HandleVisibility', 'off');
axis equal;
% daspect([1 1 1]);
% X = cellipse(xhat(:,end), P(:,:,end), 2);
% minv = min(X,[],2);
% maxv = max(X,[],2);
% axis([minv(1), maxv(1), minv(2), maxv(2)]);
legend;
%%
pause(1);
N = numel(time);
for i=1:N
    X = cellipse(xhat(:,i), P(:,:,i), 2);
    set(cellipse_plt, 'XData', X(1,:), 'YData', X(2,:));
    set(estimate_plt, 'XData', xhat(1,i), 'YData', xhat(2,i));
    set(vline, 'XData', [time(i) time(i)]);
    tobj.String = {'State-Space', ['Processed Measurements = ', num2str(i-1)]};
    pause(1e-3);
end
%% Batch Processing of Data

% N = numel(time);
% zb = z(:);
% Hb = repmat(H,N,1);
% Rb = R;
% for i=1:N-1
%     Rb = blkdiag(Rb,R);
% end
% % Rb = sparse(Rb)
% % full column rank formulation
% Pni = inv(P0) + (Hb.')/Rb * Hb;
% Pn = inv(Pni)
% xn = (Pn/P0)*xhat0 + ((Pni\Hb.')/Rb)*zb
% 
% xhat(:,end)
% P(:,:,i)
% r
