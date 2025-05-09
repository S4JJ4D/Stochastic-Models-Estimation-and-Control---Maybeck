% ------------- MATHEMATICA STYLE PLOTTING ---------------
% Mathematica Style:
% 'FontName', 'Source Code Pro'
% Area Colors: 'FaceColor', '#DAD9EB', 'FaceColor', '#F6E1BE'
% Edge Colors: 'EdgeColor', '#4A457F', 'EdgeColor', '#E5A73C'
% Axes Style:
% set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
%     'XAxisLocation', 'origin', 'YAxisLocation', 'origin', ...
%     'Box', 'on', ...
%     'FontName', 'Source Code Pro');
% set([ax.XAxis, ax.YAxis], 'Color', '#737373', 'TickLabelColor', 'k',...
%    'LineWidth',.8);
% set([ax.XLabel, ax.YLabel], 'Color', 'k')
% --------------------------------------------------------

%%
clear; close all;
t1 = [];
t2 = [];
for i=0:8
    for j=0:9
        if j > i
            t1 = [t1, i];
            t2 = [t2, j];
        end
    end
end
cla;
a_plt = plot(t1, t2, 'bo', 'MarkerFaceColor', 'w');
hold on; grid on; box on;
axis([0, 10, 0, 10]);
ax = gca;
set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Box', 'on', ...
    'FontName', 'Source Code Pro');
set([ax.XAxis, ax.YAxis], 'Color', '#737373', 'TickLabelColor', 'k',...
   'LineWidth',.8);
set([ax.XLabel, ax.YLabel], 'Color', 'k')
xlabel('t_1');
ylabel('t_2');
axis equal;

%% Part a)
[p0,q,r] = deal(1);
tf = 10;

pf = @(t1,t2) ...
     ((p0 + t1*q) * r^2 + (t2-t1)*q*r*(p0 + t1*q + r))/((p0 + t1*q)*r + ((t2-t1)*q + r)*(p0+t1*q+r)) + ...
     (tf-t2)*q;

pf_vec = zeros(1,numel(t1));
for i=1:numel(t1)
    pf_vec(i) = pf(t1(i), t2(i));
    text(t1(i)+.15,t2(i),num2str(pf_vec(i), '%.2f'), 'FontName', 'Source Code Pro');
end
axis([-0.3108    8.8892   -0.3345    9.6655])
title('p(t_f) for each measurement time pair (t_1, t_2)');

%% Part b)
cla;
[p0,q,r] = deal(1);
pthreshold = 6.25;
pf_vec = nan(1,numel(t1));
for i=1:numel(t1)
    p1_prior = p0 + t1(i)*q;
    if p1_prior <= pthreshold
        p1_after = (p1_prior*r)/(p1_prior + r);

        p2_prior = p1_after + (t2(i) - t1(i))*q;
        if p2_prior <= pthreshold
           p2_after = (p2_prior*r)/(p2_prior + r);
        else
            continue;
        end
        pf_vec(i) = p2_after + (tf - t2(i))*q;
        pf_plt_vec(i) = plot(t1(i), t2(i), 'ro', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');
        text(t1(i)+.15,t2(i),num2str(pf_vec(i), '%.2f'), 'FontName', 'Source Code Pro');
    end
end
title({'p(t_f) for each measurement time pair (t_1, t_2) under given constraint', ...
    'p_0 = q = r = 1'});
%% Part c)
cla;
[p0,q,r] = deal(9,1,1);
pthreshold = 6.25;
pf_vec = nan(1,numel(t1));
for i=1:numel(t1)
    p1_prior = p0 + t1(i)*q;
    if p1_prior <= pthreshold
        p1_after = (p1_prior*r)/(p1_prior + r);

        p2_prior = p1_after + (t2(i) - t1(i))*q;
        if p2_prior <= pthreshold
           p2_after = (p2_prior*r)/(p2_prior + r);
        else
            continue;
        end
        pf_vec(i) = p2_after + (tf - t2(i))*q;
        pf_plt_vec(i) = plot(t1(i), t2(i), 'ro', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');
        text(t1(i)+.15,t2(i),num2str(pf_vec(i), '%.2f'), 'FontName', 'Source Code Pro');
    end
end
title({'p(t_f) for each measurement time pair (t_1, t_2) under given constraint', ...
    'p_0 = 9, q = r = 1'});

%% Part d)
cla;
[p0,q,r] = deal(1,9,1);
pthreshold = 6.25;
pf_vec = nan(1,numel(t1));
for i=1:numel(t1)
    p1_prior = p0 + t1(i)*q;
    if p1_prior <= pthreshold
        p1_after = (p1_prior*r)/(p1_prior + r);

        p2_prior = p1_after + (t2(i) - t1(i))*q;
        if p2_prior <= pthreshold
           p2_after = (p2_prior*r)/(p2_prior + r);
        else
            continue;
        end
        pf_vec(i) = p2_after + (tf - t2(i))*q;
        pf_plt_vec(i) = plot(t1(i), t2(i), 'ro', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');
        text(t1(i)+.15,t2(i),num2str(pf_vec(i), '%.2f'), 'FontName', 'Source Code Pro');
    end
end
title({'p(t_f) for each measurement time pair (t_1, t_2) under given constraint', ...
    'q = 9, p_0 = r = 1'});

%% Part e)
cla;
[p0,q,r] = deal(1,1,9);
pthreshold = 6.25;
pf_vec = nan(1,numel(t1));
for i=1:numel(t1)
    p1_prior = p0 + t1(i)*q;
    if p1_prior <= pthreshold
        p1_after = (p1_prior*r)/(p1_prior + r);

        p2_prior = p1_after + (t2(i) - t1(i))*q;
        if p2_prior <= pthreshold
           p2_after = (p2_prior*r)/(p2_prior + r);
        else
            continue;
        end
        pf_vec(i) = p2_after + (tf - t2(i))*q;
        pf_plt_vec(i) = plot(t1(i), t2(i), 'ro', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');
        text(t1(i)+.15,t2(i),num2str(pf_vec(i), '%.2f'), 'FontName', 'Source Code Pro');
    end
end
title({'p(t_f) for each measurement time pair (t_1, t_2) under given constraint', ...
    'r = 9, p_0 = q = 1'});