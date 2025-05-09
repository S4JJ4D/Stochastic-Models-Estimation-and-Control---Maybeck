function [X,k,Z] = cellipse(mu, P, c, varargin)
%C-Ellipse: Covariance Ellipse/Ellipsoid of a Bivariate/Trivariate Normal Distribution
% given the covariance matrix and the mean vector of a bivariate/trivariate
% Gaussian distributoin, the c-ellipse is returned:
% equation of the ellipse/ellipsoid is given by
% (x-mu)' inv(P) (x-mu) - c = 0
% where X is the 2D/3D position vector [x,y]/[x,y,z]
%
% X(1,:) is the x-component of the ellipse
% X(2,:) is the y-component of the ellipse
%
% q1,q2, ..., qn are the unit eigenvectors scaled by their corresponding
% eigenvalues l1,l2,...ln. As such qi = vi * li;
%
% c=1 corresponds the  q-ellipse: ellipsoid passing through q1, q2, ...
% c=2 corresponse the 2q-ellipse: ellipsoid passing through 2*q1, 2*q2, ...
% ...
%
% Table: The probability of x falling into the c-ellipse:
%
%             c=1        c=2        c=3
%           _______    _______    _______
%
%     1D    0.68269    0.95450    0.99730  <--- the well-known 68-95-99 sequence
%     2D    0.39347    0.86466    0.98889
%     3D    0.19875    0.73854    0.97071
%
%
% [X,k,Z] = cellipse(mu,P,c)
% [C,k] = cellipse(mu,P,c,Z)
%
%
%
% example 2D Case:
%
% P = [2.9087, 2.4783;
%     2.4783, 6.0913]
% mu = [1;2];
%
% X = cellipse(mu, P, 1);
% patch('XData',X(1,:),'YData',X(2,:), 'FaceColor', 'r', 'FaceAlpha', .2);
%
% example: 3D Case:
%
% P = [2.9087, 2.4783; 1.12;
%      2.4783, 6.9100; 3.068;
%      1.1200, 3.0680; 4.532];
% mu = [1;2;3];
% [X,k,Z] = cellipse(xhat0, P0, 3); % c=3 contains 0.97071 probability
% ts = trisurf(k,X(:,1),X(:,2),X(:,3),...
%     'Facecolor','red','FaceAlpha',0.5,'LineStyle','-');
% Once Z is known, you can call it as
%     [X,k] = cellipse(xhat(:,i),P(:,:,i),3,Z);
%     set(ts, 'Faces', k, 'Vertices', X);
%

if ~iscolumn(mu)
    mu = mu.';
end

[V,L] = eig(P);
R = (V /(sqrt(L))).';
Ri = V*sqrt(L);

% 1. define y = x - mu
% 2. define z = 1/sqrt(c) * L^(-1/2) * V.' * y ==> y = sqrt(c) * V * L^(1/2) * z
% 3. conclude that (x-mu)' inv(P) (x-mu) - c = 0 can be written as z.' * z = 1
% 4. parameterize z: 2D case: z = [cos(theta);sin(theta)],
%                    3D case: z = [cos(u).*sin(v);sin(u).*sin(v);cos(v)]

if size(P,1) == 2
    theta = 0:.05:2*pi;
    z = [cos(theta);sin(theta)];
    X = (sqrt(c) * Ri * z) + mu;
    nargoutchk(1,1)
else % 3D case
%     fx = @(u,v) sqrt(c) * (Ri(1,1) * cos(u).*sin(v) + Ri(1,2) * sin(u).*sin(v) + Ri(1,3) * cos(v)) + mu(1);
%     fy = @(u,v) sqrt(c) * (Ri(2,1) * cos(u).*sin(v) + Ri(2,2) * sin(u).*sin(v) + Ri(2,3) * cos(v)) + mu(2);
%     fz = @(u,v) sqrt(c) * (Ri(3,1) * cos(u).*sin(v) + Ri(3,2) * sin(u).*sin(v) + Ri(3,3) * cos(v)) + mu(3);
    if nargin == 4
        % It is assumed that Z is supplied in the fourth argument
        Z = varargin{1};
    else
        % One-time computation and generation of Z matrix
        meshN = 20;
        urange = linspace(0,2*pi,meshN);
        urange = urange(1:end-1);
        vrange = linspace(0,pi,meshN);
        Z = zeros(3,numel(vrange),numel(urange));
        i=1;
        for u=urange
            Z(:,:,i) = [cos(u).*sin(vrange);sin(u).*sin(vrange);cos(vrange)];
            i = i+1;
        end
        Z = reshape(Z,3,numel(vrange)*numel(urange));
        Z(:, [meshN*(1:8)+1, meshN*(2:9)]) = []; % removing duplicate points

        if nargout ~=3
            warning("It is advised to capture Z as the 3rd output" + ...
                " and to use it in subsequent calls of this function.");
        end
    end
    X = ((sqrt(c) * Ri * Z) + mu).';
    k = convhull(X);
end


end