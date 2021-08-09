function J = eval_obj(output, reference, optU)
% Manually evaluates MPC objective function
dU                      = [optU(1); optU(2:end) - optU(1:end-1)];
tracking_error          = (output - reference)' * (output - reference);
controller_movement     = dU' * dU;
J                       = tracking_error + controller_movement; 