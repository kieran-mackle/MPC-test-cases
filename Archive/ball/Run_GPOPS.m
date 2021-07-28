% ===================== Simple 1D Ball GPOPS  ===================== %
% ================================================================= %
clearvars; clc;
tic;
% ----------------------------------------------------------------- %
%                    Construct bounds and guess                     %
%------------------------------------------------------------------ %
% Phase 1: Position change
i = 1;
bounds.phase(i).initialtime.lower  = 0;
bounds.phase(i).initialtime.upper  = 0;
bounds.phase(i).finaltime.lower    = 0;
bounds.phase(i).finaltime.upper    = 50;

bounds.phase(i).initialstate.lower = [0, 0, 0];         % Initial state
bounds.phase(i).initialstate.upper = [0, 0, 0];

bounds.phase(i).state.lower        = [-10, -10, -1.5];  % State constraints
bounds.phase(i).state.upper        = [10, 10, 1.5];

bounds.phase(i).finalstate.lower   = [5,0,0];           % Final state
bounds.phase(i).finalstate.upper   = [5,0,0];

bounds.phase(i).control.lower      = [-1];             % Control rate constraints
bounds.phase(i).control.upper      = [1];

guess.phase(i).control             = [1.5; 0];

guess.phase(i).time                = [0; 30];
guess.phase(i).state(:,1)          = [0; 5];
guess.phase(i).state(:,2)          = [0; 0];
guess.phase(i).state(:,3)          = [0.8; 0.8];

bounds.phase(i).integral.lower     = zeros(1,1);
bounds.phase(i).integral.upper     = 10e4*ones(1,1);
guess.phase(i).integral            = zeros(1,1);

% Phase 1: Position hold
i = 2;
bounds.phase(i).initialtime.lower  = 0;
bounds.phase(i).initialtime.upper  = 50;
bounds.phase(i).finaltime.lower    = 0;
bounds.phase(i).finaltime.upper    = 100;

bounds.phase(i).initialstate.lower = [5, 0, 0];
bounds.phase(i).initialstate.upper = [5, 0, 0];

bounds.phase(i).state.lower        = [-10, -10, -1.5];
bounds.phase(i).state.upper        = [10, 10, 1.5];

bounds.phase(i).finalstate.lower   = [5,0,0];
bounds.phase(i).finalstate.upper   = [5,0,0];

bounds.phase(i).control.lower      = bounds.phase(1).control.lower;
bounds.phase(i).control.upper      = bounds.phase(1).control.upper;

guess.phase(i).control             = [0; 0];

guess.phase(i).time                = [13; 33];
guess.phase(i).state(:,1)          = [5; 5];
guess.phase(i).state(:,2)          = [0; 0];
guess.phase(i).state(:,3)          = [0; 0];

bounds.phase(i).integral.lower     = zeros(1,1);
bounds.phase(i).integral.upper     = 10e3*ones(1,1);
guess.phase(i).integral            = zeros(1,1);


% Phase 1-2 continuity constraints
bounds.eventgroup(1).lower         = [zeros(1,4)];
bounds.eventgroup(1).upper         = [zeros(1,4)];

bounds.eventgroup(2).lower         = [20];
bounds.eventgroup(2).upper         = [20];


% ----------------------------------------------------------------- %
%                    Construct GPOPS-II input                       %
%------------------------------------------------------------------ %
setup.name                  = '1D-ball';
setup.functions.continuous  = @multiphase_ball_dynamics;
setup.functions.endpoint    = @GPOPS_objective;
setup.nlp.solver            = 'ipopt';
setup.bounds                = bounds;
setup.guess                 = guess;
% setup.auxdata               = auxdata;
setup.displaylevel          = 2;


% ----------------------------------------------------------------- %
%                       Call GPOPS-II solver                        %
%------------------------------------------------------------------ %
GPOPS_out = gpops2(setup);
toc;
% ----------------------------------------------------------------- %
%                            Plot Output                            %
%------------------------------------------------------------------ %
% figure(4);
% clf;
% figure(5);
% clf;
% figure(6);
% clf;

for phase = 1:2
solution        = GPOPS_out.result.solution.phase(phase);

output.time     = solution.time;
output.control  = solution.control;
output.state    = solution.state;
output.Z        = solution.state(:,1);

figure(1);
set(gcf,'color','w');
hold on;
grid on;
title('GPOPS input');
ylabel('Input force (N)');
xlabel('Time (s)');
% stairs(output.time, output.control, 'k-')
stairs(output.time, output.state(:,3), 'b-');

figure(2);
set(gcf,'color','w');
sgtitle('Physical state trajectory');
subplot(2,1,1);
hold on;
grid on;
plot(output.time, output.state(:,1), 'b-');
ylabel('x-position (m)');
subplot(2,1,2);
hold on;
grid on;
plot(output.time, output.state(:,2), 'b-');
ylabel('velocity (m)');
xlabel('Time (s)');


figure(3);
set(gcf,'color','w');
hold on;
title('Output');
grid on;
plot(output.time, output.Z, 'b-');
ylabel('x-position (m)');
xlabel('Time (s)');

end