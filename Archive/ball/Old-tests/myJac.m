function J = myJac(varargin)
%MYJAC A callable function for the numerical Jacobian evaluation of a
%GPOPS-II dynamics function input script.
%   This function numerically evaluates the Jacobian of a function f, about
%   the point x.
%  
%   J = myJac(f, x, ID)
%       Inputs:
%           For an anlytical function input:
%           - f: A function of x
%           - x: State vector
%           If a GPOPS-II dynamics function is inputted:
%           - f: dynamics function handle (eg @Ellip6DOF)
%           - x: GPOPS-II input struct
%           - ID: distinguish to find Jacobian for state or control 
%
%       Output: 
%           - J: Numerical Jacobian of f(x) about x
%
%       If there are three inputs, the code will assume that the supplied
%       function, f, is a GPOPS-II dynamics function, and so will
%       automatically use struct notation to index the dynamics output.
%
%   Example usage:
%       Analytical function input:
%           f = @(x)[2*x(1)^3+x(2)-x(3), 5*x(1)+x(2)^2, 2*x(2)+3*x(3)^3];
%           x = [1,-2,2];
%           J = myJac(f,x)
% 
%           J = 
%               6.0000    1.0000   -1.0000
%               5.0000   -4.0000         0
%                    0    2.0000   36.0000
%   
%       GPOPS-II dynamics function input:
%           A = myJac(@Ellip6DOF,input,'state')
%           B = myJac(@Ellip6DOF,input,'control')

f           = varargin{1};
x           = varargin{2};
f0          = f(x);
pert        = 1e-9;


if length(varargin) == 2
    flag  = 0;
else
    flag = 1;
    ID          = varargin{3};
    input       = x;
    f0          = f0.dynamics;
    x           = input.phase.(ID);
end

x_len               = length(x);
f_len               = length(f0);
J                   = zeros(f_len, x_len);

for i = 1:x_len
        x_pert              = x;
        x_pert(i)           = x(i) + pert;
        if flag == 1
            input.phase.(ID)    = x_pert;
            f_pert              = f(input);
            f1                  = f_pert.dynamics;
        else
            f1                  = f(x_pert);
        end
        partial             = (f1 - f0) / pert;
        J(:,i)              = partial;
end

end

