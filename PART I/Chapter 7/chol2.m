function L = chol2(A, type)
%CHOL2 Naive implementation of chol2 algorithm
%
%   The algorithm is described in 
%   {Maybeck, P. S. (1979). Stochastic models, estimation, and control (Vol. 1). Academic press, pp. 371-372}
%
%   L = CHOL2(A) returns the Cholesky factor of A as a lower-triangular
%   matrix.
%
%   U = CHOL2(A, 'upper') returns the Cholesky factor of A as an
%   upper-triangular matrix.
%
%   Example:
%      A = [1 2 3;2 8 2;3 2 14];
%      chol2(A)
%      chol2(sym(A),'upper')
%

if nargin == 1
    type = 'lower';
    % the other option is 'upper'
end

n = size(A,1);
if isa(A, 'sym')
    L = sym(zeros(n));
else
    L = zeros(n);
end

if strcmpi(type, 'lower')

    for i=1:n
        for j=1:i
            if i~=j
                sum = 0;
                for k=1:j-1, sum=sum+(L(i,k)*L(j,k)); end
                L(i,j) = 1/L(j,j) * (A(i,j) - sum);
            else
                sum = 0;
                for k=1:i-1, sum=sum+(L(i,k))^2;end
                L(i,j) = sqrt(A(i,i) - sum);
            end
        end
    end

else

    for j=n:-1:1
        for i=j:-1:1
            if i~=j
                sum=0;
                for k=j+1:n, sum=sum+(L(i,k)*L(j,k)); end
                L(i,j) = 1/L(j,j) * (A(i,j) - sum);
            else
                sum=0;
                for k=j+1:n, sum=sum+(L(j,k))^2; end
                L(i,j) = sqrt(A(j,j) - sum);
            end
        end
    end

end

end









