function output = get_reference(mpc_input)
% Provides controller with reference and cost matrices

% Unpack input
t           = mpc_input.t;
% set_points  = mpc_input.set_points;

% Update reference and cost matrices
if t < 0.6
    r           = [0, 0];
    mpc_input.cost.output   = 1e0*[1,0; 0,1];
    
elseif (t >= 0.6) && (t < 4.3)
    r           = [0, 1.4032];
    mpc_input.cost.output   = 1e4*[0,0; 
                                   0,1]; % Only track velocity
else 
    r = [5, 0];
    mpc_input.cost.output   = 1e5*[1,0; 
                                   0,1]; 
end

% Modify MPC reference
mpc_input.reference   = r;

output      = mpc_input;