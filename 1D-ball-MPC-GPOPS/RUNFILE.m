% ======================================================================= %
%                           Constrained Ball MPC                          %
% ======================================================================= %
% This is the main runfile for this directory.
% The case attempts to integrate the GPOPS minimum time manouevre with the
% MPC controller, by providing a velocity reference to track. 



% ----------------------------------------------------------------------- %
%                                RUN GPOPS
% ----------------------------------------------------------------------- %
% Run GPOPS to get optimal reference
% -------------------------------------
Run_GPOPS;

% Extract GPOPS output reference
% -------------------------------------
t = gpops_output.time;
x = gpops_output.state(:,1);
p = polyfit(t,x,5);
dp = polyder(p);

% Calculate climb phase velocity (pen on paper)
% ----------------------------------------------
% figure();
% plot(t, polyval(p, t));

phase1_t = [0, 1.3];
phase1_x = [0, 0];

phase2_t = [phase1_t(2), 4.3];
phase2_x = [phase1_x(2), 5];

phase3_t = [phase2_t(2), t(end)];
phase3_x = [phase2_x(2), 5];


% Visualise phases of flight
% -------------------------------------
figure(3);
hold on;
grid on;
title("Velocity reference");
plot(phase1_t, phase1_x, 'r');
plot(phase2_t, phase2_x, 'r');
plot(phase3_t, phase3_x, 'r');


% From the plot above, there are three phases to the trajectory, two of
% which are constant state. In the climb phase, the velocity is:
v_phase2 = (phase2_x - phase1_x)/(phase2_t - phase1_t);
% 1.4032e+00


% Run MPC Controller
% -------------------------------------
% ball_mpc;




