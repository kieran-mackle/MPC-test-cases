function output = ball_dynamics(input)

X       = input.phase.state;
Fdot    = input.phase.control;
m       = input.auxdata.mass; 
k       = input.auxdata.k;

s       = X(:,1);
v       = X(:,2);

F       = X(:,3);
a       = (F - k*s)./m;

s_dot   = v;
v_dot   = a;

output.dynamics     = [s_dot, v_dot, Fdot];


end