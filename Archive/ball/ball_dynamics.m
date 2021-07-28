function output = ball_dynamics(input)

X       = input.phase.state;
U       = input.phase.control;
m       = input.auxdata.mass; 

s       = X(:,1);
v       = X(:,2);

F       = U(:,1);
a       = F./m;

s_dot   = v;
v_dot   = a;

output.dynamics = [s_dot, v_dot];

end