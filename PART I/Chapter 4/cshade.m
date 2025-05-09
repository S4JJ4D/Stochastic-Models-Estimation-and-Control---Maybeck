function sp = cshade(t, x_l, x_u)
%cshade Shading the area between two curves
%   for the given curves (t,x1) and (t,x2) on the plane, the space between
%   them is shaded.
%
%   Example:
%   f = @(t) cos(t) + sin(t).^2 + 2 + t - log(t+1) + exp(-.9*(t+1));
%   t = (0:.02:10).'; x = f(t);
%   g = @(t) f(t) - (1 + .08*sin(5*t)); h = @(t) f(t) + (1 + .12*cos(6*t));
%   sp = cshade(t, g(t), h(t));


if ~iscolumn(t)
    t = t.';
end
if ~iscolumn(x_l)
    x_l = x_l.';
end
if ~iscolumn(x_u)
    x_u = x_u.';
end

XData = [t;flipud(t)];
YData = [x_u;flipud(x_l)];
% shaded patch
sp = patch('XData', XData', 'YData', YData, ...
    'FaceColor', 'k', 'FaceAlpha', .1, 'LineStyle', 'none');

end