% ===================== Simple 1D Ball GPOPS  ===================== %
% ================================================================= %
% Single phase GPOPS example with 1D ball
% Differences to MPC include:
%   - dynamics function: contains an extra state (Fdot) to allow
%     constraint of rates in GPOPS

mpc_output = output;
clearvars -except mpc_output; clc;

auxdata.mass                    = 10;

% ----------------------------------------------------------------- %
%                    Construct bounds and guess                     %
%------------------------------------------------------------------ %
bounds.phase.initialtime.lower  = 0;
bounds.phase.initialtime.upper  = 0;
bounds.phase.finaltime.lower    = 30;
bounds.phase.finaltime.upper    = 30;

bounds.phase.initialstate.lower = [0, 0, 0];         % Initial state
bounds.phase.initialstate.upper = [0, 0, 0];

bounds.phase.state.lower        = [-10, -10, -10];  % State constraints
bounds.phase.state.upper        = [10, 10, 10];

bounds.phase.finalstate.lower   = [-10,-10,-10];           % Final state
bounds.phase.finalstate.upper   = [10,10,10];

bounds.phase.control.lower      = [-100];             % Control rate constraints
bounds.phase.control.upper      = [100];

bounds.phase.integral.lower     = zeros(1,1);
bounds.phase.integral.upper     = 10e5*ones(1,1);

guess.phase.control             = [1.5; 0];
guess.phase.time                = [0; 30];
guess.phase.state(:,1)          = [0; 5];
guess.phase.state(:,2)          = [0; 0];
guess.phase.state(:,3)          = [0.8; 0.8];

guess.phase.integral            = zeros(1,1);

% ----------------------------------------------------------------- %
%               Provide mesh refinement and method                  %
%------------------------------------------------------------------ %
% mesh.method                 = 'hp-PattersonRao';
mesh.tolerance              = 1e-5;
mesh.maxiterations          = 8;
mesh.colpointsmin           = 5;
mesh.colpointsmax           = 10;
M                           = 15;
mesh.phase.fraction         = (1/M)*ones(1,M);
mesh.phase.colpoints        = 2*ones(1,M);

% ----------------------------------------------------------------- %
%                    Construct GPOPS-II input                       %
%------------------------------------------------------------------ %
setup.name                  = '1D-ball';
setup.functions.continuous  = @ball_dynamics;
setup.functions.endpoint    = @GPOPS_objective;
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
% figure(4);
% clf;
% figure(5);
% clf;
% figure(6);
% clf;

for phase = 1:1
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
legend('MPC', 'GPOPS')

figure(2);
set(gcf,'color','w');
sgtitle('Physical state trajectory');
subplot(2,1,1);
hold on;
grid on;
plot(output.time, output.state(:,1), 'b-');
ylabel('x-position (m)');
legend('MPC', 'GPOPS')

subplot(2,1,2);
hold on;
grid on;
plot(output.time, output.state(:,2), 'b-');
ylabel('velocity (m)');
xlabel('Time (s)');
legend('MPC', 'GPOPS')

figure(3);
set(gcf,'color','w');
hold on;
title('Output');
grid on;
plot(output.time, output.Z, 'b-');
ylabel('x-position (m)');
xlabel('Time (s)');
legend('MPC', 'GPOPS')

end