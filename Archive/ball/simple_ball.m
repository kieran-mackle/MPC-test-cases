% ======================================================================= %
%                           Constrained Ball MPC                          %
% ======================================================================= %

clearvars;
tic;
deg = pi/180;
rad = 180/pi;

% ----------------------------------------------------------------------- %
% Define MPC Parameters
% ----------------------------------------------------------------------- %
params.timestep     = 0.5;
params.horizon      = 10;
params.sim_time     = 50;

% ----------------------------------------------------------------------- %
% Define function handles
% ----------------------------------------------------------------------- %
functions.dynamics  = @ball_dynamics;
functions.output    = @ball_output;

% ----------------------------------------------------------------------- %
% Define initial state
% ----------------------------------------------------------------------- %
x0  = [0, 0];
u0  = [0.0]';

initial.state       = x0;
initial.control     = u0;

target_ref          = [5]';

auxdata.mass        = 10; 

% ----------------------------------------------------------------------- %
% Define cost and constraint matrices
% ----------------------------------------------------------------------- %
cost.output         = eye(1);
cost.control        = eye(1);

penalty_method      = 'quadratic';
penalty_weight      = 1e0; % For slack variables
hard_only           = true;

constraints.hard.rate   = [-1, 1];
constraints.hard.input  = [-1.5, 1.5];
constraints.hard.output = [0, 8];

constraints.soft.rate   = [-1, 1];
constraints.soft.input  = [-1.5, 1.5];
constraints.soft.output = [0, 7];

% ----------------------------------------------------------------------- %
% Define model mismatch auxdata
% ----------------------------------------------------------------------- %
model_mismatch          = false;
mismatch_auxdata.mass   = 3;

% ----------------------------------------------------------------------- %
% Solve MPC QP Problem
% ----------------------------------------------------------------------- %
output = solve_slack_mpc(initial, params, cost, constraints, functions, ...
                         auxdata, target_ref, penalty_method,           ...
                         penalty_weight, hard_only,                     ...
                         model_mismatch, mismatch_auxdata);

toc;
% ----------------------------------------------------------------------- %
% Extract results
% ----------------------------------------------------------------------- %
figure(1);
clf;
figure(2);
clf;
figure(3);
clf;

figure(1);
set(gcf,'color','w');
hold on;
grid on;
title('MPC input');
ylabel('Input force (N)');
xlabel('Time (s)');
stairs(output.time, output.control, 'b-');
% patch('vertices', [0, constraints.soft.input(1); 
%                    0, constraints.hard.input(1); 
%                    params.sim_time, constraints.hard.input(1); 
%                    params.sim_time, constraints.soft.input(1)], ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.input(1), ...
%      constraints.hard.input(1)], 'r-', 'LineWidth', 2);
% patch('vertices', [0, constraints.soft.input(2); 
%                    0, constraints.hard.input(2); 
%                    params.sim_time, constraints.hard.input(2); 
%                    params.sim_time, constraints.soft.input(2)], ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.input(2), ...
%      constraints.hard.input(2)], 'r-', 'LineWidth', 2);

figure(2);
set(gcf,'color','w');
sgtitle('Physical state trajectory');
subplot(2,1,1);
hold on;
grid on;
plot(output.time, output.state(1:2:end), 'b-');
ylabel('x-position (m)');
subplot(2,1,2);
hold on;
grid on;
plot(output.time, output.state(2:2:end), 'b-');
ylabel('velocity (m)');
xlabel('Time (s)');


figure(3);
set(gcf,'color','w');
hold on;
title('Output');
grid on;
plot(output.time, output.Z, 'b-');
% plot([0,params.sim_time], [target_ref, target_ref], 'b--');
% patch('vertices', [0, constraints.soft.output(1); 
%                    0, constraints.hard.output(1); 
%                    params.sim_time, constraints.hard.output(1); 
%                    params.sim_time, constraints.soft.output(1)], ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.output(1), ...
%      constraints.hard.output(1)], 'r-', 'LineWidth', 2);
% patch('vertices', [0, constraints.soft.output(2); 
%                    0, constraints.hard.output(2); 
%                    params.sim_time, constraints.hard.output(2); 
%                    params.sim_time, constraints.soft.output(2)], ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.output(2), ...
%      constraints.hard.output(2)], 'r-', 'LineWidth', 2);
ylabel('x-position (m)');
xlabel('Time (s)');

