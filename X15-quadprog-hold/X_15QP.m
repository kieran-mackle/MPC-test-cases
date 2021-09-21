% ======================================================================= %
%                              X-15 MPC Control                           %
% ======================================================================= %
addpath '/home/kieran/Documents/MATLAB/MPC'

clearvars;
deg = pi/180;
rad = 180/pi;

% ----------------------------------------------------------------------- %
% Define MPC Parameters
% ----------------------------------------------------------------------- %
params.timestep     = 0.5;
params.horizon      = 75;
params.sim_time     = 20;
convex_solver       = 'quadprog';

% ----------------------------------------------------------------------- %
% Define initial state
% ----------------------------------------------------------------------- %
load('./../../Optimisation/Results/Config1/6DOF/20kmAltHold/20kmAltHold.mat')

out = output.result.solution.phase;
x0  = out.state(1,:);
u0  = out.control(1,:);

initial.state       = x0;
initial.control     = u0;

% Define control model
control_model.auxdata       = auxdata;
control_model.dynamics      = @Ellip6DOF;
control_model.output        = @X15_outputs;
reference_function          = @X15_reference;

% Define plant model
plant_model                 = control_model;


% ----------------------------------------------------------------------- %
% Define cost and constraint matrices
% ----------------------------------------------------------------------- %
cost_weightings.output         = 1e3* [1, 0, 0, 0, 0;   % Altitude
                                       0, 0, 0, 0, 0;   % FDA
                                       0, 0, 0, 0, 0;   % THR
                                       0, 0, 0, 0, 0;   % Ma
                                       0, 0, 0, 0, 1];  % FPA
cost_weightings.control        = eye(length(initial.control));

% Quadprog constraint handling options
constraints.type = 'soft';      % None, soft, hard or mixed
penalty_method = 'quadratic';   % Quadratic or linear
penalty_weight      = 1e3;


constraints.hard.rate    = [-10*deg, 10*deg;
                            -0.2,   0.2];
constraints.hard.input   = [-10*deg, 10*deg;
                            -0.2,   0.2];
constraints.hard.output  = [   0,      0     ;
                            -40*deg, 40*deg  ;
                               0,      1     ;
                               0,      0     ;
                               0,      0     ];

constraints.soft.rate    = [-10*deg, 10*deg;
                            -0.2,   0.2];
constraints.soft.input   = [-10*deg, 10*deg;
                            -0.2,   0.2];
constraints.soft.output  = [   0,      0     ;
                            -40*deg, 40*deg  ;
                               0,      1     ;
                               0,      0     ;
                               0,      0     ];
                           
% NaN for hard constraints
constraints.weights.hard_rate   = [0, 0;
                                   0, 0];
constraints.weights.hard_input  = [0, 0;
                                   0, 0];
constraints.weights.hard_output = [0, 0;
                                   0, 0;
                                   0, 0;
                                   0, 0;
                                   0, 0];

% ----------------------------------------------------------------------- %
% Construct MPC Input
% ----------------------------------------------------------------------- %
mpc_input.control_model     = control_model;
mpc_input.cost              = cost_weightings;
mpc_input.constraints       = constraints;
mpc_input.penalty_method    = penalty_method;
mpc_input.penalty_weight    = penalty_weight;
mpc_input.params            = params;
mpc_input.initial           = initial;
mpc_input.solver            = convex_solver;

% ----------------------------------------------------------------------- %
% Construct Forward Simulation Input
% ----------------------------------------------------------------------- %
dyn_input.phase.state       = initial.state;
sim_input.dynamics_input    = dyn_input;
sim_input.plant_model       = plant_model;
sim_input.Ts                = params.timestep;

% ----------------------------------------------------------------------- %
% Construct Input
% ----------------------------------------------------------------------- %
input.mpc_input             = mpc_input;
input.sim_input             = sim_input;
input.reference_function    = reference_function;

% ----------------------------------------------------------------------- %
% Solve MPC QP Problem
% ----------------------------------------------------------------------- %
output                      = mpc_control(input);

% ----------------------------------------------------------------------- %
% Plot Results
% ----------------------------------------------------------------------- %
% Control Inputs
% ------------------
figure(1); clf; set(gcf,'color','w');
sgtitle('Optimal MPC input sequence');
subplot(1,2,1);
hold on;
grid on;
stairs(output.time, output.Z(2, :)*rad);
% row = 2; scaler = rad;
% patch('vertices', [0, constraints.soft.output(row, 1); 
%                    0, constraints.hard.output(row, 1); 
%                    params.sim_time/scaler, constraints.hard.output(row, 1); 
%                    params.sim_time/scaler, constraints.soft.output(row, 1)]*scaler, ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.output(row, 1), ...
%      constraints.hard.output(row, 1)]*scaler, 'r-', 'LineWidth', 2);
% patch('vertices', [0, constraints.soft.output(row, 2); 
%                    0, constraints.hard.output(row, 2); 
%                    params.sim_time/scaler, constraints.hard.output(row, 2); 
%                    params.sim_time/scaler, constraints.soft.output(row, 2)]*scaler, ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.output(row, 2), ...
%      constraints.hard.output(row, 2)]*scaler, 'r-', 'LineWidth', 2);
title('Flap deflection angle');
xlabel('Time (s)');
ylabel('\delta_f (deg)');

subplot(1,2,2);
hold on;
grid on;
stairs(output.time, output.Z(3, :)*100);
% row = 3; scaler = 100;
% patch('vertices', [0, constraints.soft.output(row, 1); 
%                    0, constraints.hard.output(row, 1); 
%                    params.sim_time/scaler, constraints.hard.output(row, 1); 
%                    params.sim_time/scaler, constraints.soft.output(row, 1)]*scaler, ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.output(row, 1), ...
%      constraints.hard.output(row, 1)]*scaler, 'r-', 'LineWidth', 2);
% patch('vertices', [0, constraints.soft.output(row, 2); 
%                    0, constraints.hard.output(row, 2); 
%                    params.sim_time/scaler, constraints.hard.output(row, 2); 
%                    params.sim_time/scaler, constraints.soft.output(row, 2)]*scaler, ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.output(row, 2), ...
%      constraints.hard.output(row, 2)]*scaler, 'r-', 'LineWidth', 2);
title('Thrust throttle setting');
xlabel('Time (s)');
ylabel('\delta_T (%)');


% System Outputs
% ------------------
figure(2);
clf; set(gcf,'color','w');
sgtitle('Simulated output trajectory');
subplot(1,3,1);
hold on;
grid on;
plot(output.time, output.Z(1,:)*1e-3, 'k-');
% plot([0,params.sim_time], [target_ref(1), target_ref(1)]/1e3, 'b--');
title('Altitude');
xlabel('Time (s)');
ylabel('h (km)');

subplot(1,3,2);
hold on;
grid on;
plot(output.time, output.Z(4,:), 'k-');
% plot([0,params.sim_time], [target_ref(4), target_ref(4)], 'b--');
title('Mach number');
xlabel('Time (s)');
ylabel('Ma');

subplot(1,3,3);
hold on;
grid on;
plot(output.time, output.Z(5,:)*rad, 'k-');
% plot([0,params.sim_time], [target_ref(5), target_ref(5)]*rad, 'b--');
title('Flight path angle');
xlabel('Time (s)');
ylabel('\gamma (deg)');





