function [S2, K, xhat2] = pcsru(S, H, R, xhat, z, options)
%PCSRU The Potter covariance square root update algorithm
%
%   The algorithm is described in
%   {Maybeck, P. S. (1979). Stochastic models, estimation, and control (Vol. 1). Academic press, pp. 373-374}
%
%   S2 = PCSRU(S, H, R) uses the square root factor of the covariance
%   matrix at time t_i^-, just before the measurement takes place, i.e., S1
%   along with the measurement matrix H and noise covariance R and returns
%   the square root factor of the covariance matrix at time t_i^+, just
%   after the measurement.
%
%   [S2,K] = PCSRU(...) also outputs the kalman gain matrix for the case of
%   'Andrews' one-time vector update algorithm. For the potter algorithm K
%   is returned as a null matrix.
%
%   [S2,K,xhat2] = PCSRU(S, H, R, xhat, z) incorporates the measurement
%   vector z and updates the state estimate xhat.
%
%   [...] = PCSRU(..., 'Method', METHOD) uses the specified method for the
%   update of covarinace square roots.
%   METHOD can be:
%           'Potter'   -  uses sequential use of the scalar Potter method.
%           'Andrews'  -  uses the one-time vector update method of Andrews.
%
%   Example:
%
%       S = randi([-10 20], 5, 5);
%       H = randi([-10, 20], 3, 5);
%
%       [Q,~] = qr(randi([-10,20],3,3));  % generate a random orthonormal matrix
%       R = Q*diag([1 2 3])*Q.';          % generate a SPD matrix
%       R = 1/2 * (R + R.');              % making sure that R is symmetric
%
%       S2 = pcsru(S,H,R,'Method','Potter');
%
%       z = 10*rand(3,1);                 % generate some observations
%       xhat = 15*rand(5,1);              % generate initial state estimate
%
%       [S2,~,xhat2] = pcsru(S,H,R,xhat,z,'Method','Andrews');
%

arguments
    S (:,:) {mustBeNonNan, mustBeA(S, {'numeric', 'sym'}), mustBeSquareMatrix(S)}
    H (:,:) {mustBeNonNan, mustBeA(H, {'numeric', 'sym'}), mustBeOfSameNumberOfColumns(H,S,{'H','S'})}
    R (:,:) {mustBeNonNan, mustBeA(R, {'numeric', 'sym'}), mustBeSquareMatrix(R), mustBeOfSameNumberOfRows(H,R,{'H','R'}), mustBeSPD(R,'R')}

    xhat (:,1) {mustBeNonNan, mustBeA(xhat, {'numeric', 'sym'})} = [];
    z    (:,1) {mustBeNonNan, mustBeA(z, {'numeric', 'sym'})} = [];

    options.Method (1,:) char {mustBeMember(options.Method, {'Andrews', 'Potter'})} = 'Potter';
    options.PrintIntermediateResults (1,1) logical = false;
end

[m,n] = size(H);
% m: number of measurements
% n: number of states

% initialize outputs
K = zeros(n,m);

updateStateEstimate = false;
if ~isempty(xhat)
    validateattributes(xhat,{'numeric','sym'},{'numel',n},'','xhat')
    validateattributes(z,{'numeric','sym'},{'numel',m},'','z')
    updateStateEstimate = true;
end


if isa([S,zeros(n,m);H,R], 'sym') || isa(xhat, 'sym') || isa(z, 'sym')
    S = sym(S);
    H = sym(H);
    R = sym(R);
    xhat = sym(xhat);
    z = sym(z);
end

if strcmp(options.Method, 'Potter')
    % Potter method: The most efficient means of performing a vector
    % measurement update is to employ the Potter scalar update repeatedly
    % m times.

    % Check if R is diagonal
    if isa(R, 'sym')
        R_isDiag = all(isAlways(diag(diag(R)) == R), 'All');
    else
        R_isDiag = isdiag(R);
    end

    if ~R_isDiag
        V = chol2(R);
        H = V\H;
        R = eye(m);
    end

    % iterate through measurements
    Si = zeros(n,n,m);
    if isa(S, 'sym')
        Si = sym(zeros(n,n,m+1));
    end
    Si(:,:,1) = S;

    for i=1:m
        [Si(:,:,i+1), k] = cspcsru(Si(:,:,i), H(i,:), R(i,i));
        if updateStateEstimate
            xhat = xhat + k*(z(i) - H(i,:) * xhat);
        end
    end

    S2 = Si(:,:,end);
    xhat2 = xhat;

    if options.PrintIntermediateResults
        Si
    end

else
    % Andrews method: processing an m-dimensional measurement vector in a
    % single update, without requiring diagonalization
    V = chol2(R);
    A = S.'*H.';
    Sigma = chol2(A.'*A + R);
    K = (S*A/(Sigma.'))/Sigma;

    if updateStateEstimate
        xhat2 = xhat + K*(z-H*xhat);
    end
    S2 = S - (S*A/(Sigma.'))/(Sigma+V) * A.';
end

end

% ---------------------------------------------------------------------------------
%% Helper function
function [S2, k] = cspcsru(S, H, R)
%CSPCSRU The Classical scalar form of the Potter covariance square root update algorithm
a = S.' * H.';
b = 1/(a.' * a + R);
gamma = 1/(1+sqrt(b*R));
k = b*S*a;
S2 = S - gamma*k*a.';
end

% ---------------------------------------------------------------------------------
%% Custom validation function
function mustBeSquareMatrix(A)
[m,n] = size(A);
if m ~= n
    erridType = 'mustBeSquareMatrix:NonSquareMatrix';
    msgType = 'Matrix must be square.';
    throwAsCaller(MException(erridType,msgType))
end
end

function mustBeOfSameNumberOfColumns(A, B, inputnames)
if nargin == 2
    inputnames = {'IN(1)','IN(2)'};
end

[~,N] = size(A);
[~,n] = size(B);
if n ~= N
    erridType = 'mustBeOfSameNumberOfColumns:InconsistentDimensions';
    msgType = [inputnames{1}, ' and ', inputnames{2}, ...
        ' must have the same number of columns.'];
    throwAsCaller(MException(erridType,msgType))
end
end

function mustBeOfSameNumberOfRows(A, B, inputnames)
if nargin == 2
    inputnames = {'IN(1)','IN(2)'};
end
[M,~] = size(A);
[m,~] = size(B);
if m ~= M
    erridType = 'mustBeOfSameNumberOfRows:InconsistentDimensions';
    msgType = [inputnames{1}, ' and ', inputnames{2}, ...
        ' must have the same number of rows.'];
    throwAsCaller(MException(erridType,msgType))
end
end

function mustBeSPD(A, inputnames)
% must be symmetric positive definite (SPD)
if nargin == 1
    inputnames = 'IN';
end

if isa(A, 'sym')
    A_isSymmetric = all(isAlways(tril(A) == triu(A)),'All');
else
    A_isSymmetric = issymmetric(A);
end

if ~A_isSymmetric
    erridType = 'mustBeSPD:NonSymmetricMatrix';
    msgType = [inputnames, ' must be a symmetric matrix.'];
    throwAsCaller(MException(erridType,msgType))
else
    sigma = eig(A);
    if any(double(real(sigma))<=0) % if there are any non-positive eigenvalues
        erridType = 'mustBeSPD:NonPDMatrix';
        msgType = [inputnames, ' must be a positive definite matrix.'];
        throwAsCaller(MException(erridType,msgType))
    end
end
end


