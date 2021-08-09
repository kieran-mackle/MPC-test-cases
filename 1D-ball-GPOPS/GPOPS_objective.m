function output = GPOPS_objective(input)
% tf      = input.phase.finaltime;

J = input.phase(1).integral;

output.objective = J;

end