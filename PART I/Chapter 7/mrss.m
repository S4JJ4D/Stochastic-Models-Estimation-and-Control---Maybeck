function S2 = mrss(Phi,S,Gd,Qd)
%MRSS The matrix RSS method
%
%   The algorithm is described in
%   {Maybeck, P. S. (1979). Stochastic models, estimation, and control (Vol. 1). Academic press, pp. 377}
%
%   S2 = MRSS(Phi,S,Gd,Qd) uses the square root factor of the covariance
%   matrix at time t_i^+, just after the measurement, i.e., S1 and returns
%   the square root factor of the covariance matrix at time t_(i+1)^-, just
%   before the next measurement.
%
%   See also CHOL2.

X = Phi*S;
P = X*X.' + (Gd*Qd*Gd.');
S2 = chol2(P);

end

