function output = get_reference(mpc_input)

% Unpack input
t           = mpc_input.t;
set_points  = mpc_input.set_points;

% Update reference and cost matrices
if t < 5.5
    s_ref       = 0;
    v_ref       = polyval(set_points, t);
    r           = [s_ref, v_ref];
    
    mpc_input.cost.output   = 1e5*[0,0; 
                               0,1]; % Only track velocity
else
    r = [5,0];
    mpc_input.cost.output   = 1e5*[1,0; 
                               0,1]; 
end

mpc_input.reference   = r;

output      = mpc_input;