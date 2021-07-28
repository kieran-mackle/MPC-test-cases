function output = GPOPS_objective(input)
% t1_0 = input.phase(1).initialtime;
t1_f = input.phase(1).finaltime;
% X1_0 = input.phase(1).initialstate;
X1_f = input.phase(1).finalstate;

t2_0 = input.phase(2).initialtime;
t2_f = input.phase(2).finaltime;
X2_0 = input.phase(2).initialstate;
% X2_f = input.phase(2).finalstate;



J = input.phase(1).integral + input.phase(2).integral;

output.objective = J;

% Linkage constraints:
output.eventgroup(1).event = [X2_0(1:3) - X1_f(1:3), t2_0 - t1_f];
output.eventgroup(2).event = [t2_f - t2_0];

end