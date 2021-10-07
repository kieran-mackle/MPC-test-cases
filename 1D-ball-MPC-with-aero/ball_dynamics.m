function output = ball_dynamics(input)

X       = input.phase.state;
U       = input.phase.control;
m       = input.auxdata.mass; 

% State variables
s       = X(:,1);
v       = X(:,2);

% Transfer funtion
Ma_eq   = 6 + 0.01*v.^3;

% Forces and acceleration
[~,Cd,~] = GetAero(input.auxdata, 0, Ma_eq, 0);
D       = 0.5 * 1.225 * v.^2 * 1 * Cd;
F       = U(:,1) - D;
a       = F./m;

% Output
s_dot   = v;
v_dot   = a;

output.dynamics = [s_dot, v_dot];

end