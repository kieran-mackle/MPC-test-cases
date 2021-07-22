function J = numerical_jac(varargin)

f           = varargin{1};  % Function handle
x           = varargin{2};  % Function input
f0          = f(x);         % Initial function value
pert        = 1e-9;         % Perturbation

% Initialisation
x_len       = length(x);
f_len       = length(f0);
J           = zeros(f_len, x_len);

for i = 1:x_len
        x_pert              = x;
        x_pert(i)           = x(i) + pert;

        f1                  = f(x_pert);

        partial             = (f1 - f0) / pert;
        J(:,i)              = partial;
end

end