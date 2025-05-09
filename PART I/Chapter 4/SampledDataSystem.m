%%
% Useful Functions:
% sldsim, imp2exp, lft, iconnect
%%
clear;
close all;

t_f = 3;
Tc = 1e-3; % cont. time time-step for simulation purposes.
Ts = .1;  % the computer samples data at this particular step-time

% Plant is a continuous-time LTI system (a physical system is continuous in
% nature)
P = tf(1, [1 1 0]);
% A controller is designed in continuous-time:
C = tf(70*[1 2], [1 10]);
% Simulate the system:
tc = 0:Tc:t_f;
cl = feedback(series(C, P), 1);
yc = step(cl, tc);

fig = figure('Units', 'Normalized');
fig.Position = [0.1583 0.1925 0.7297 0.7002];
tl = tiledlayout(2,2,"TileSpacing","compact","Padding","compact");
nexttile(2);
plot(tc, yc, 'DisplayName', ['Plant output in analog setting (driven by the analog controller)', newline, ...
    'not real and only shown as a baseline']);
hold on;
title('Outputs');

% In a computer-controlled system, we treat every signal in the closed-loop
% system as discrete-time signal. By placing the plant inside the control
% loop:

Pd = c2d(P, Ts, 'zoh'); % place the plant into the control loop

% Controller is discretized and placed inside the loop (using forward Euler
% method)
Cd = tf(70*[1 2*Ts-1], [1 10*Ts-1], Ts);

% Form the digital control loop:
cld = feedback(series(Cd, Pd), 1);
% Since all signals in a digital control system are sampled at constant
% rate, data is only available at times separated by Ts
td = 0:Ts:t_f;
yd = step(cld, td');
plot(td, yd, '*', 'DisplayName', ...
    ['$y[n]$: Sampled plant output in the digital control system', ...
    newline, ...
    ': plant output seen by the computer: output picked by the sensor']);
hold on;

% To access u[n], we cannot use 'cld' object, because it represents the
% resulting closed-loop control system where internal signals, like u[n]
% and e[n] are inaccessible. To extract u[n], we use 'imp2exp' approach:
A = [1 0 -Cd 0;-Pd 1 0 0;0 -1 -1 1];
yidx = [1;2;3];  % [u;y;e] 
uidx = 4;        % [r] 
Sys = imp2exp(A,yidx,uidx); 
Sys.InputName = {'r'}; 
Sys.OutputName = {'u';'y';'e'};
ud = step(Sys(1), td');

% search for the word 'resampling' in the following comments which
% discusses a way to retrieve sampled-time signals within the lft control
% structure.


% Set up the continuous interconnection and calculate the sampled data response with sdlsim.
 %
 % Following control system is the target for 'sdlsim' function. sdlsim is used 
 % to emulate the sampled-data (hybrid) control system where the plant is a
 % continuous-time system and the controller is a digital system. Plant is
 % preceeded by a ZOH (D/A) and followed by A/D before contacting the
 % controller.
 %
 % In fact, 'sdlsim' is able to "simulate" the following sampled data (hybrid)
 % control system where <Ts> is the specified sample time for the digital
 % control system and <Tc> is chosen to be much smaller than <Ts> to emulate
 % the continuous-time nature of the plant and other analog components.
 %
 %
 %                               ┌────────────────────────────────────┐
 %                               │                                    │
 %                               │                                    │
 %                               │          <Ts>                      │           <Tc>
 %         <Ts>                  │      ┌───────────┐        <Ts>     │      ┌────────────┐
 % r(t)   ┌─────┐ r[n]   ┌───┐   │  e[n]│           │ u[n]  ┌─────┐   │  u(t)│            │ y(t)
 %  ─────►│ A/D ├───────►│SUM├───┼──────┤    C(z)   ├──────►│ ZOH ├───┼─────►│    G(s)    ├───┐
 %        └─────┘      + └─▲─┘   │      │           │       └─────┘   │      │            │   │
 %                        -│     │      └───────────┘                 │      └────────────┘   │
 %                         │     │                                    │                       │
 %                         │     │                                    │                       │
 %                         │     │                                    │                       │
 %                         │     │                                    │                       │
 %                         │     │      <Ts>                          │                       │
 %                         │     │ y[n]┌─────┐                        │                       │
 %                         └─────┼─────┤ A/D ├────────────────────────┼───────────────────────┘
 %                               │     └─────┘                        │
 %                               │                                    │
 %                               └────────────────────────────────────┘
 %                                              COMPUTER
 %
 %
 %
 %
 %
 % To adapt to sdlsim input arguments, the following equivalent control
 % system is used: 
 % w is the exposed input and v is the exposed output.
 % 
 % The overall structure of the control system is shown here:
 %
 %
 %                                                                        
 %                         ┌─────────────────┐                        
 %                u1       │     ██████╗     │       y1               
 %      ──────────────────▶│     ██╔══██╗    ├───────────────────────▶
 %                         │     ██████╔╝    │       y2               
 %                         │     ██╔═══╝     ├───────────────────┐    
 %                u2       │     ██║         │       y3          │    
 %           ┌────────────▶│     ╚═╝         ├────────────────┐  │    
 %           │             └─────────────────┘                │  │    
 %           │                                                │  │    
 %           │                                                │  │    
 %           │                                                │  │    
 %           │                                                │  │    
 %           │                                                │  │    
 %           │                                                │  │    
 %           │                                                │  │    
 %           │                                                │  │    
 %           │                                                │  │    
 %           │                                                │  │    
 %           │             ┌─────────────────┐                │  │    
 %           │             │                 │                │  │    
 %           │             │    ██╗  ██╗     ├──┬┬┬┬─────┐    │  │    
 %           │  ┌─────┬──┬┬┤    ██║ ██╔╝     │◀─┼┼┤│ A/D │◀───┘  │    
 %           └──┤ ZOH │◀─┼┼┤    █████╔╝      ├──┴┴┴┴─────┘       │    
 %              └─────┴──┴┴┤    ██╔═██╗      ├──┬┬┬┬─────┐       │    
 %                         │    ██║  ██╗     │◀─┼┼┤│ A/D │◀──────┘    
 %                         │    ╚═╝  ╚═╝     ├──┴┴┴┴─────┘            
 %                         └─────────────────┘                        
 %                                                                    
 %                                                                        
 %
 % The internal plumbing in systems P and K is as follows:
 %
 %
 %                            M
 %                     ┌─────────────┐
 %                     │             │
 %         w=r=u1      │        ┌─┐  │    y1=y(t)=v
 %      ───────────────┼──┐  ┌─►│P├──┼────────────────────►
 %                     │  │  │  └─┘  │
 %                     │  │  │       │
 %                     │  │  │       │
 %                     │  │  │       │    y2
 %                     │  └──┼───────┼─────────────────────┐
 %                     │     │       │                     │
 %                     │     │       │                     │
 %             u2      │     │  ┌─┐  │    y3               │
 %      ┌──────────────┼─────┴─►│P├──┼────────────────┐    │
 %      │              │        └─┘  │                │    │
 %      │              │             │                │    │
 %      │              └─────────────┘                │    │
 %      │                                             │    │
 %      │                                             │    │
 %      │                                             │    │
 %      │                                             │    │
 %      │                                             │    │
 %      │                    K                        │    │
 %      │           ┌───────────────────┐             │    │
 %      │           │                   │             │    │
 %      │           │           ┌──┐    ├────┬─────┐  │    │
 %      │           │     +-----┤-C│<---┤◄├┼┼│ A/D │◄─┘    │
 %      │           │     |     └──┘    ├────┴─────┘       │
 %      │           │     |             │                  │
 %      │  ┌───┬────┤   ┌─▼─┐           │                  │
 %      └──┤ZOH│◄├┼┼┤<--┤SUM│           │                  │
 %         └───┴────┤   └─▲─┘           │                  │
 %                  │     |     ┌──┐    ├────┬─────┐       │
 %                  │     +-----┤ C│<---┤◄├┼┼│ A/D │◄──────┘
 %                  │           └──┘    ├────┴─────┘
 %                  │                   │
 %                  └───────────────────┘
 % 
 % 
 % All the signals in this diagram, except for the ones shown by the symbol
 % ◄├┼┼ (signals sampled with Ts) are included in the output of 'sdlsim' function.
 % To access discrete-time signals, we can resample u2, y2, and y3, at
 % sampling frequence of the controller, Ts.
 % %%% Resampling of a sdlsim signal named 'x':
 % %%% idx = 1:Ts/Tc:numel(x);
 % %%% xd = x(idx);
 % 
 % Notes on 'sdlsim' function:
 % 1. It is required that nusys > nu and nysys > ny, where
 % nusys: number of input signals to P
 % nysys: number of output signals from P
 % nu: number of output signals from K
 % ny: number of input signals to K
 % 
 %
 %
 %           ╭────────────────╮          
 %   nusys   │                │ nysys    
 % ═════════▶│       P        ╠════════▶ 
 %           │                │          
 %           ╰────────────────╯          
 %
 %                                      
 %           ╭────────────────╮          
 %   nu      │                │ ny       
 %  ◀════════╣       K        │◀═════════
 %           │                │          
 %           ╰────────────────╯          
 %
 %
 %
 % 2. 'sdlsim' automatically connects internal signals, exposing 'nusys-nu'
 % input signals and 'nysys-ny' output signals. (There's also an
 % undocumented way of calling lft which uses names of the signals for
 % connections. Execute the command 'edit lft' in command window to see
 % more information. I tried to use it, but it seems to have issues.)
 % This is an example where
 % nusys = 6;
 % nysys = 5;
 % nu = 2;
 % ny= 3;
 % In sdlsim terminology,
 % Exposed inputs are     w = [u1, u2, u3, u4];
 % Exposed outputs are   vt = [y1 y2];
 % Inputs to K are       yt = [y3 y4 y5];      % notice the order of entries
 % Outputs from K are    ut = [u5 u6];         % notice the order of entries
 %
 %                          ┌────────────┐
 %                w1 = u1   │            │
 % ────────────────────────►│            │
 %                w2 = u2   │            │y1 = v1
 % ────────────────────────►│            ├─────────────────────►
 %                w3 = u3   │            │y2 = v2
 % ────────────────────────►│      P     ├─────────────────────►
 %                w4 = u4   │            │y3
 % ────────────────────────►│            ├─────────┐
 %                     u5   │            │y4       │
 %                ┌────────►│            ├──────┐  │
 %                │    u6   │            │y5    │  │
 %                │   ┌────►│            ├───┐  │  │
 %                │   │     └────────────┘   │  │  │
 %                │   │                      │  │  │
 %                │   │                      │  │  │
 %                │   │                      │  │  │
 %                │   │     ┌────────────┐   │  │  │
 %                │   └─────┤            │◄──┘  │  │
 %                │         │            │      │  │
 %                └─────────┤      K     │◄─────┘  │
 %                          │            │         │
 %                          │            │◄────────┘
 %                          │            │
 %                          └────────────┘
 %
 % ------------------------------------------------------------------------------------
 % ------------------------------------------------------------------------------------
 %
 % To design K, we note that the following control system,
 %
 %
 %                 ┌────────┐       ┌─────────┐
 % r    ┌───┐      │        │       │         │
 %  ───►│SUM├─────►│  C(s)  ├──────►│   G(s)  ├───┐
 %    + └─▲─┘      │        │       │         │   │
 %      - │        └────────┘       └─────────┘   │
 %       y│                                       │
 %        │                                       │
 %        │                                       │
 %        │                                       │
 %        └───────────────────────────────────────┘
 %
 %
 % is in fact, equivalent to the following system where K is a
 % 2-input, 1-output transfer matrix:
 %
 % 
 % 
 %            r          ┌────────┐       ┌─────────┐
 %            ──────────►│        │       │         │
 %                       │  K(s)  ├──────►│   G(s)  ├───┐
 %                ┌─────►│        │       │         │   │
 %              y │      └────────┘       └─────────┘   │
 %                │                                     │
 %                │                                     │
 %                │                                     │
 %                │                                     │
 %                └─────────────────────────────────────┘
 %
 %
 % K can be expressed as [C(s) -C(s)]:
 %
 %
 %                  ┌───────────┐
 %        r         │           │
 % ────────────────►│   C(s)    ├─────────────┐
 %                  │           │             │
 %                  └───────────┘             │
 %                                          + │            ┌──────────┐
 %                                         ┌──▼──┐         │          │
 %                                         │ SUM ├────────►│   G(s)   ├───┐
 %                                         └──▲──┘         │          │   │
 %                  ┌───────────┐           + │            └──────────┘   │
 %        y         │           │             │                           │
 %        ┌────────►│   -C(s)   ├─────────────┘                           │
 %        │         │           │                                         │
 %        │         └───────────┘                                         │
 %        │                                                               │
 %        │                                                               │
 %        │                                                               │
 %        │                                                               │
 %        │                                                               │
 %        └───────────────────────────────────────────────────────────────┘
 %
 %


t = tc'; 
w = ones(size(t)); 

K = [Cd -Cd]; % extended controller to match sdlsim structure: A 2-input, 1-output system
M = [0 P;1 0;0 P]; % extended plant to match sdlsim structure: A 2-input, 3-output system

% Using Tc to specify integration time of the continuous-time system
[vt,yt,ut,t] = sdlsim(M,K,w,t,[],[],[],Tc); 
% actual plant output is unobservable in real life. It is only inspectable
% in simulation settings.
plot(vt{1},vt{2}, 'LineWidth', 1.1, 'DisplayName', '\bf $y(t)$: Actual plant output in the digital control system');
xlabel('Time');
ylabel('Amplitude');
grid on;
legend('Location','southeast', 'Interpreter', 'latex');

nexttile(1);
plot(td, ud, 'r*', 'DisplayName', '$u[n]$: Output of the digital controller');
hold on;
title('Inputs');
plot(ut{1},ut{2}, 'b-', 'LineWidth', 1, 'DisplayName', ['\bf $u(t)$: Output of the digital controller after ZOH:', newline, ...
    '\bf right before injecting to the plant']);
xlabel('Time');
ylabel('Amplitude');
grid on;
legend('Interpreter','latex');

nexttile(3, [1, 2]);
text(0,0,...
[ ...
'                                  <Ts>                         <Tc>',               newline, ...
'         <Ts>                 ┌───────────┐     <Ts>      ┌────────────┐',          newline, ...
' r(t)   ┌─────┐     ┌───┐e[n] │           │u[n]┌─────┐u(t)│            │y(t)',      newline, ...
'  ─────►│ A/D ├────►│SUM├────►│    C(z)   ├───►│ ZOH ├───►│    G(s)    ├───┐',      newline, ...
'        └─────┘   + └─▲─┘     │           │    └─────┘    │            │   │',      newline, ...
'                     -│       └───────────┘               └────────────┘   │',      newline, ...
'                      │     <Ts>                                           │',      newline, ...
'                      │    ┌─────┐                                         │',      newline, ...
'                      └────┤ A/D ├─────────────────────────────────────────┘',      newline, ...
'                       y[n]└─────┘',                                                newline, ...
 ], 'FontName', 'Cascadia Code');

axis([-0.1641, 0.8359, -0.4689, 0.5311]);
box on;
ax = gca;
set(ax, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

title(tl, ...
    ['Bold signals are cont. and are not observable by digital computer in their cont. form', ...
    newline, sprintf('Ts = %.3f', Ts)]);