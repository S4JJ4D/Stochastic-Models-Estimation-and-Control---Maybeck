function B = householder(varargin)
%HOUSEHOLDER The householder transformation algorithm
%
%   The algorithm is described in
%   {Maybeck, P. S. (1979). Stochastic models, estimation, and control (Vol. 1). Academic press, pp. 382-383}
%
%   B = HOUSEHOLDER(A)
%
%   S = HOUSEHOLDER(Phi,S,Gd,Wd)
%
%   Example:
%
%       A = [2 1;2 3;1 1;0 sqrt(2)];
%

if nargin == 1
    A = varargin{1};
elseif nargin == 4
    A = [varargin{1}*varargin{2}, varargin{3}*varargin{4}].';
else
    error('Invalid number of input arguments');
end


[m,n] = size(A); % m=n+s

if isa(A, 'sym')
    uzero = sym(zeros(m,1));
    yzero = sym(zeros(n,1));
else
    uzero = zeros(m,1);
    yzero = zeros(n,1);
end

for k=1:n
    a = sqrt(sum(A(k:end,k).^2)) * sign(A(k,k));
    
    d = 1/(a*(a + A(k,k)));

    u = uzero;
    u(k) = a + A(k,k);
    u(k+1:end) = A(k+1:end,k);

    y = yzero;
    y(k) = 1;
    y(k+1:end) = d*A(:,k+1:end).'*u;

    A = A - u*y.';
end

B = A(1:n,1:n).';

end