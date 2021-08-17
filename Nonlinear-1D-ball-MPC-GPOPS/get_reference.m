function output = get_reference(mpc_input)
% Provides controller with reference and cost matrices

% Unpack input
% t           = mpc_input.t;
% set_points  = mpc_input.set_points;
% 
% % Update reference and cost matrices
% if t < 5.5
%     s_ref       = 0;
%     v_ref       = polyval(set_points, t);
%     r           = [s_ref, v_ref];
%     
%     mpc_input.cost.output   = 1e5*[0,0; 
%                                    0,1]; % Only track velocity
% else
%     r = [5,0];
%     mpc_input.cost.output   = 1e5*[1,0; 
%                                    0,1]; 
% end

% Nominal reference
r = [5,0];

% Current state
X = ball_output(mpc_input.initial.state);

% Offset error
e = r - X;

% Modify reference
% r_m = r + 0.1*e;

r_m = r + 5*e * [0,1;0,0];

mpc_input.reference   = r_m;

output      = mpc_input;