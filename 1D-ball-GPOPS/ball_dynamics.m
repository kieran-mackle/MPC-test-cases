function output = ball_dynamics(input)

X       = input.phase.state;
Fdot    = input.phase.control;
m       = input.auxdata.mass; 

s       = X(:,1);
v       = X(:,2);

F       = X(:,3);
a       = F./m;

s_dot   = v;
v_dot   = a;

ref_errorQ2         = (5 - s) .* (5 - s);
controller_movement = Fdot .* Fdot;

output.dynamics     = [s_dot, v_dot, Fdot];
output.integrand    = [ ref_errorQ2 + controller_movement ];


end