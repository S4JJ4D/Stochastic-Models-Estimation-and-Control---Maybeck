clear;

a_data = [...
0.039447732	2.008032129;
4.773175542	2.851405622;
7.337278107	3.453815261;
9.861932939	3.935742972;
14.87179487	5.060240964;
19.8816568	6.104417671;
22.4852071	6.586345382;
24.77317554	7.108433735;
30.01972387	8.072289157;
];

b_data = [...
2.04724409400000	5.16000000000000;
6.06299212600000	9;
10.1968503900000	9.92000000000000;
10.1574803100000	15.9600000000000;
5.15748031500000	15.2400000000000;
2.99212598400000	18.1200000000000;
3.93700787400000	26.2000000000000;
11.1023622000000	25.2000000000000;
12.7952755900000	29.1600000000000;
7.28346456700000	23.2800000000000;
];

c_data = [...
15	10;
15	15;
10	15;
20	15;
25	20;
20	20;
15	20;
10	20;
5	20;
10	25;
15	25;
20	25;
15	30;
];

d_data = [...
2.036659878	28.96414343;
8.065173116	24.82071713;
6.150712831	21.87250996;
7.698574338	17.88844622;
13.89002037	20.91633466;
19.30753564	20.67729084;
14.82688391	15.93625498;
23.13645621	13.8247012;
17.92260692	9.960159363;
24.03258656	6.015936255;
28.02443992	6.812749004;
];

figs = {...
    struct('name', 'Fig (A)', 'data', ...
    table(a_data(:,1), a_data(:,2), 'VariableNames', {'x', 'y'}), ...
    'N', size(a_data,1)), ...
    struct('name', 'Fig (B)', 'data', ...
    table(b_data(:,1), b_data(:,2), 'VariableNames', {'x', 'y'}), ...
    'N', size(b_data,1)), ...
    struct('name', 'Fig (C)', 'data', ...
    table(c_data(:,1), c_data(:,2), 'VariableNames', {'x', 'y'}), ...
    'N', size(c_data,1)), ...
    struct('name', 'Fig (D)', 'data', ...
    table(d_data(:,1), d_data(:,2), 'VariableNames', {'x', 'y'}), ...
    'N', size(d_data,1))
}; 
%%

for i = 1:numel(figs)
    % 1. mean estimate
    figs{i}.M = 1/figs{i}.N * sum(figs{i}.data{:,:});

    % 2. variance estimate
    figs{i}.V = 1/(figs{i}.N-1) * sum(figs{i}.data{:,:}.^2) - ...
        figs{i}.N/(figs{i}.N-1) * figs{i}.M.^2;

    % 3. covariance estimate
    figs{i}.C = 1/(figs{i}.N-1) * sum(figs{i}.data.x .* figs{i}.data.y) -...
        figs{i}.N/(figs{i}.N-1) * figs{i}.M(1)*figs{i}.M(2);

    % 4. correlation coefficient estimate
    figs{i}.r = figs{i}.C/sqrt(figs{i}.V(1) * figs{i}.V(2));
end

close all;
fg = figure('Units', 'normalized');
fg.Position = [0.1427 0.1630 0.4906 0.7102];
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for i =1:numel(figs)
    nexttile;
    plot(figs{i}.data.x, figs{i}.data.y, 'ko', 'markerfacecolor', 'k', ...
        'DisplayName', 'Data Points');
    hold on;
    xc = linspace(min(figs{i}.data.x), max(figs{i}.data.x), 20);
    plot(xc, figs{i}.C/figs{i}.V(1) * (xc-figs{i}.M(1)) + figs{i}.M(2), ...
        'LineWidth', 1, 'DisplayName', 'y on x');
    if abs(figs{i}.C) <= 1e-8
        % note that the regression line passes through the mean point
        % (Mx,My):
        plot([figs{i}.M(1), figs{i}.M(1)], ...
            [min(figs{i}.data.y), max(figs{i}.data.y)], ...
            'LineWidth', 1, 'DisplayName', 'x on y');
    else
        plot(xc, figs{i}.V(2)/figs{i}.C * (xc-figs{i}.M(1)) + figs{i}.M(2), ...
            'LineWidth', 1, 'DisplayName', 'x on y');
    end
    axis equal; grid on; box on;
    xlabel('x'); ylabel('y');
    title(string(figs{i}.name) + " $$r_{xy} = " + string(figs{i}.r) + "$$", ...
        'Interpreter', 'latex');
    legend('Interpreter', 'latex', 'FontSize', 10);
end

