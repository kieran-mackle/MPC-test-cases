function output = msd_dynamics(input)

X       = input.phase.state;
U       = input.phase.control;
m       = input.auxdata.mass;
k       = input.auxdata.k;
c       = input.auxdata.c;

s       = X(:,1);
v       = X(:,2);

% ADD SPRING AND DAMPER

F       = U(:,1) - k*s - c*v;
a       = F./m;

s_dot   = v;
v_dot   = a;

output.dynamics = [s_dot, v_dot];

end