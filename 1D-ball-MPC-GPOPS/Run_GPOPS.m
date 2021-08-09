% ===================== Simple 1D Ball GPOPS  ===================== %
% ================================================================= %
% This script runs GPOPS to find the optimal control sequence for a 
% minimum time manouevre to a reference position of 5m.
%
% State variables are: position, velocity, force rate
% Objective: minimise final time.

clearvars;


auxdata.mass                    = 10;

% ----------------------------------------------------------------- %
%                    Construct bounds and guess                     %
%------------------------------------------------------------------ %
bounds.phase.initialtime.lower  = 0;
bounds.phase.initialtime.upper  = 0;
bounds.phase.finaltime.lower    = 0;
bounds.phase.finaltime.upper    = 30;

bounds.phase.initialstate.lower = [0, 0, 0];         % Initial state
bounds.phase.initialstate.upper = [0, 0, 0];

bounds.phase.state.lower        = [-10, -10, -10];  % State constraints
bounds.phase.state.upper        = [10, 10, 10];

bounds.phase.finalstate.lower   = [5,0,0];           % Final state
bounds.phase.finalstate.upper   = [5,0,0];

bounds.phase.control.lower      = [-10];             % Control rate constraints
bounds.phase.control.upper      = [10];

guess.phase.control             = [1.5; 0];
guess.phase.time                = [0; 10];
guess.phase.state(:,1)          = [0; 5];
guess.phase.state(:,2)          = [0; 0];
guess.phase.state(:,3)          = [0.8; 0.8];


% ----------------------------------------------------------------- %
%               Provide mesh refinement and method                  %
%------------------------------------------------------------------ %
% mesh.method                 = 'hp-PattersonRao';
% mesh.tolerance              = 1e-5;
% mesh.maxiterations          = 8;
% mesh.colpointsmin           = 5;
% mesh.colpointsmax           = 10;
% M                           = 15;
% mesh.phase.fraction         = (1/M)*ones(1,M);
% mesh.phase.colpoints        = 2*ones(1,M);

% ----------------------------------------------------------------- %
%                    Construct GPOPS-II input                       %
%------------------------------------------------------------------ %
setup.name                  = '1D-ball';
setup.functions.continuous  = @ball_dynamics;
setup.functions.endpoint    = @min_time;
setup.nlp.solver            = 'ipopt';
setup.bounds                = bounds;
setup.guess                 = guess;
setup.auxdata               = auxdata;
% setup.mesh                  = mesh;
setup.displaylevel          = 2;

% ----------------------------------------------------------------- %
%                       Call GPOPS-II solver                        %
%------------------------------------------------------------------ %
GPOPS_out = gpops2(setup);

% ----------------------------------------------------------------- %
%                            Plot Output                            %
%------------------------------------------------------------------ %
figure(1);
clf;
figure(2);
clf;
figure(3);
clf;

for phase = 1:1
solution        = GPOPS_out.result.solution.phase(phase);

gpops_output.time     = solution.time;
gpops_output.control  = solution.control;
gpops_output.state    = solution.state;
gpops_output.Z        = solution.state(:,1);

figure(1);
set(gcf,'color','w');
sgtitle('GPOPS input');
subplot(2,1,1);
hold on;
grid on;
ylabel('Input force (N)');
stairs(gpops_output.time, gpops_output.state(:,3), 'b-');
subplot(2,1,2);
hold on;
grid on;
stairs(gpops_output.time, gpops_output.control, 'b-');
ylabel('Force rate (N/s)');
xlabel('Time (s)');

figure(2);
set(gcf,'color','w');
sgtitle('Physical state trajectory');
subplot(2,1,1);
hold on;
grid on;
plot(gpops_output.time, gpops_output.state(:,1), 'b-');
ylabel('x-position (m)');
% legend('MPC', 'GPOPS')

subplot(2,1,2);
hold on;
grid on;
plot(gpops_output.time, gpops_output.state(:,2), 'b-');
ylabel('velocity (m)');
xlabel('Time (s)');
% legend('MPC', 'GPOPS')

figure(3);
set(gcf,'color','w');
hold on;
title('Output');
grid on;
plot(gpops_output.time, gpops_output.Z, 'b-');
ylabel('x-position (m)');
xlabel('Time (s)');
% legend('MPC', 'GPOPS')

end

