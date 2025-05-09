function nutvec = nugrid(t1, t2, n)
%NUGRID Nonuniform time grid.
%   NUGRID(t1, t2, N) generates N points between t1 and t2 (end-points are inclusive).
%
%   Example:
%   tvec = nugrid(0, 10, 15); plot(tvec, zeros(1,numel(tvec)), 'r*'); yline(0, 'k-');

% The following approach seems to generate random subintervals that look
% more 'natural' compared to a simple usage of a rand function over the
% entire interval [t1, t2].

tspan = [t1,t2];
tpoints = n; % total number of time points in the time-space, including end
% points: t0=tspan(1), t1, t2, ..., tn = tspan(1)

utvec = linspace(tspan(1),tspan(2), tpoints-1);
t1 = utvec(1:end-1);
t2 = utvec(2:end);
tmid = 1/2 * (t1+t2); % mid-points of the uniform time grid

delta = diff(utvec); 
delta = delta(1); % interval length of the uniform time-grid

sd = delta/3; % set the interval length to 3sigma ~ %99
while true
    nutvec = sort(sd * randn(1,numel(tmid)) + tmid);
    if nutvec(1) > tspan(1) && nutvec(end) < tspan(2)
        % make sure random points do not exceed the given bounds
        break;
    end
end
nutvec = [tspan(1), nutvec, tspan(2)]; % append end-points;


