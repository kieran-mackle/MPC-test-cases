% ======================================================================= %
%                           Constrained Ball MPC                          %
% ======================================================================= %
% This is the main runfile for this directory.
% The case attempts to integrate the GPOPS minimum time manouevre with the
% MPC controller, by providing a velocity reference to track. 


% Run GPOPS to get optimal reference
% -------------------------------------
Run_GPOPS;


% Run MPC Controller
% -------------------------------------
ball_mpc;




