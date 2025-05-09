function B = mgs(varargin)
%MGS The modified Gram-Schmidt algorithm
%
%   The algorithm is described in
%   {Maybeck, P. S. (1979). Stochastic models, estimation, and control (Vol. 1). Academic press, pp. 380-381}
%
%   B = MGS(A)
%
%   S = MGS(Phi,S,Gd,Wd)
%
%   Example:
%
%       A = [2 1;2 3;1 1;0 sqrt(2)];
%       B = MGS(A);
%

if nargin == 1
    A = varargin{1};
elseif nargin == 4
    A = [varargin{1}*varargin{2}, varargin{3}*varargin{4}].';
else
    error('Invalid number of input arguments');
end

% In the context of square root filtering this is a (n+s)-by-n matrix
n = size(A,2);

% Initialization:
C = zeros(n,n);

if isa(A, 'sym')
    C = sym(C);
end

for k=1:n
    a = sqrt(A(:,k).' * A(:,k));

    C(k,k) = a;
    for j=k+1:n
        C(k,j) = 1/a * A(:,k).' * A(:,j);
    end
    for j=k+1:n
        A(:,j) = A(:,j) - C(k,j)*(1/a * A(:,k));
    end
end

B = C.';

end

% ---------------------------------------------------------------------------------
%% Helper function
function [S2, K] = cspcsru(S, H, R)
%CSPCSRU The Classical scalar form of the Potter covariance square root update algorithm
a = S.' * H.';
b = 1/(a.' * a + R);
gamma = 1/(1+sqrt(b*R));
K = b*S*a;
S2 = S - gamma*K*a.';
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

if ~issymmetric(A)
    erridType = 'mustBeSPD:NonSymmetricMatrix';
    msgType = [inputnames, ' must be a symmetric matrix.'];
    throwAsCaller(MException(erridType,msgType))
else
    sigma = eig(A);
    if any(real(sigma)<=0) % if there are any non-positive eigenvalues
        erridType = 'mustBeSPD:NonPDMatrix';
        msgType = [inputnames, ' must be a positive definite matrix.'];
        throwAsCaller(MException(erridType,msgType))
    end
end
end


