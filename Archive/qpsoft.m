function [X,lambda,how,epsilon]=qpsoft(H,f,A,B,vlb,vub,X,weights,normtype,verbosity)
%QPSOFT Quadratic programming with soft constraints. 
%   X=QPSOFT(H,f,A,b) solves the quadratic programming problem:
%
%            min 0.5*x'Hx + f'x  subject to:  Ax <= b 
%             x    
%
%       X=QPSOFT(H,f,A,b,VLB,VUB) defines a set of lower and upper
%       bounds on the design variables, X, so that the solution  
%       is always in the range VLB <= X <= VUB.
%
%       X=QPSOFT(H,f,A,b,VLB,VUB,X0) sets the initial starting point to X0.
%
%       [x,LAMBDA]=QPSOFT(H,f,A,b) returns the set of Lagrangian multipliers,
%       LAMBDA, at the solution.
%
%       [X,LAMBDA,HOW] = QPSOFT(H,f,A,b) also returns a string HOW that 
%       indicates error conditions at the final iteration.
%
%       [X,LAMBDA,HOW,EPSILON] =
%         QPSOFT(H,f,A,b,VLB,VUB,X0,WEIGHTS,NORMTYPE) allows
%         soft constraints to be added. The NORMTYPE used to penalise the
%         constraint violation can be chosen to be '1-norm', '2-norm',
%         'mixed-norm' (using both 1-norm and 2-norm) or 'inf-norm'. 
%         WEIGHTS is a one- or two-column matrix depending on NORMTYPE. 
%         Each row of WEIGHTS represents penalty weights associated with
%         the corresponding constraint row in Ax <= b. For the
%         'mixed-norm' the first and second columns are for the 1-norm
%         and 2-norm calculations, respectively. If NORMTYPE is
%         'inf-norm', then WEIGHTS must be a scalar.
%         The constraint violations are returned in EPSILON. 
%         The new problem that is solved for is either 
%
%              min 0.5*x'Hx + f'x  + 0.5*e'Pe + Qe 
%              x,e
%                             or
%
%              min 0.5*x'Hx + f'x  + R||e|| 
%              x,e                         inf
%  
%                    subject to:  Ax <= b + e 
%                                  e >= 0
%
%         where P is the weight vector corresponding to the 1-norm
%         violation, Q is the weight matrix corresponding to the 2-norm
%         violation and R is the scalar that weights the inf-norm violation.
%
%       X=QPSOFT(H,f,A,b,VLB,VUB,X0,WEIGHTS,NORMTYPE,DISPLAY) controls the
%       level of warning messages displayed.  All messages can be
%       turned off with DISPLAY = -1 and constraint information with
%       DISPLAY = 0.
%
%       QPSOFT produces warning messages when the solution is either unbounded
%       or infeasible. 

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.33 $  $Date: 1998/08/31 22:29:21 $
%   Andy Grace 7-9-90. Mary Ann Branch 9-30-96.

%   Original QP.M code changed in order to handle soft constraints.
%   Now also works with QUADPROG.M - change source code to enable/disable
%   Copyright (c) 1999 by Eric C. Kerrigan
%   $Revision: 1.03$  $Date: 20/12/1999 10:55:00
%   $Revision: 1.04$  $Date: 23/12/1999 13:37:00, JMM

% Handle missing arguments

% MOD ECK 03/02/99
if nargin < 10, verbosity = 1; % MOD ECK 10/3/99
  if nargin < 9, normtype = [];
    if nargin < 8, weights = [];
end, end, end
% MOD
 
if nargin < 7, X = []; 
  if nargin < 6, vub = []; 
	if nargin < 5, vlb = [];
end, end, end
[ncstr,nvars]=size(A);
nvars = max([length(f),length(H),nvars]); % In case A is empty

if isempty(verbosity), verbosity = 1; end % MOD ECK 10/3/99
if isempty(X), X=zeros(nvars,1); end
if isempty(A), A=zeros(0,nvars); end
if isempty(B), B=zeros(0,1); end

% Expect vectors
f=f(:);
B=B(:);
X=X(:);

if  norm(H,'inf')==0 | isempty(H)
  H=[]; 
  % Really a lp problem
  % caller = 'lp'; % MOD ECK 02/03/99
else
  % caller = 'qp'; % MOD ECK 02/03/99
  % Make sure it is symmetric
  if norm(H-H') > eps
	if verbosity > -1
	  disp('WARNING: Your Hessian is not symmetric.  Resetting H=(H+H'')/2')
	end
    H = (H+H')*0.5;
  end
end

% This bit is where the soft constraints functionality has been added

%%% JMM MOD, 23.12.99:
if strcmp(normtype,'inf-norm') & weights==inf,
  weights=[];
end
%%% END OF JMM MODS, 23.12.99.

if ~isempty(weights)
  
  % Check to see if all all elements of weights are non-negative
  if any(any(weights < 0))
	error('All elements of WEIGHTS must be non-negative')
  end
    
  % Check to see if weights has correct number of rows
  [nrows,ncols] = size(weights);
  if strcmp(normtype,'1-norm') | strcmp(normtype,'2-norm') | strcmp(normtype,'mixed-norm')
	if nrows ~= ncstr
	  error(['WEIGHTS must have ', num2str(ncstr), ' rows.']) 
	end
  end
  
  % Check to see if weights has correct number of columns and append if
  % necessary
  if strcmp(normtype,'1-norm')
	if ncols > 1
      error('WEIGHTS must be a column vector for NORMTYPE 1-norm.')
	end
	weights = [weights, zeros(nrows,1)];  
  elseif strcmp(normtype,'2-norm')
	if ncols > 1
	  error('WEIGHTS must be a column vector for NORMTYPE 2-norm.')
	end
	weights = [zeros(nrows,1), weights];
  elseif strcmp(normtype,'mixed-norm')
	if ncols > 2
      error('WEIGHTS can only have a maximum of two columns for NORMTYPE mixed-norm.')
	end
	if ncols== 1
	  weights = [weights, weights];
	  if verbosity > -1
		disp('WARNING: The same WEIGHTS are being used for both 1-norm and 2-norm.')
      end
    end
  elseif strcmp(normtype,'inf-norm')
	if (nrows > 1) | (ncols > 1)
	  error('WEIGHTS must be a scalar for NORMTYPE inf-norm')
	end
	weights = [weights*ones(ncstr,1), zeros(ncstr,1)];
  else
	error('NORMTYPE must be 1-norm, 2-norm, mixed-norm or inf-norm');
  end
  
  weights=weights'; % in order to work with find and any
  soft = find(any(weights~=0) & ~any(weights==inf));
  numsoft = length(soft);
  hard = find(any(weights==inf));
  numhard = length(hard);
  remove = find(~any(weights~=0));
  numremove = length(remove);
  keep = find(any(weights~=0));
  numkeep = length(keep);
  
  if strcmp(normtype,'inf-norm')
	minusones = -ones(ncstr,1);
	numsoft=1;
  else
	minusones = zeros(ncstr,numsoft);
    for i = 1:numsoft
	  minusones(soft(i),i) = -1;
    end
  end

%  if ~strcmp(normtype,'inf-norm')
	
	% output number of hard and soft constraints
	if numhard ~= 0 & verbosity > 0
	  disp('The following constraint(s) are hard:')
	  disp(hard)
    end
    if numsoft ~= 0 & verbosity > 0
      disp('The following constraint(s) are soft:')
	  disp(soft)
    end
	if numremove ~= 0  & verbosity > 0
	  disp('The following constraint(s) will be removed:')
	  disp(remove)
	end
%  end
  
  % remove constraints
  A = A(keep,:);
  B = B(keep,:);
  
  if (numsoft ~= 0)
	
	minusones = minusones(keep,:); % MOD Moved down here to fix Matlab 4
                                   % error if minusones is empty ECK 10/3/99
	
	% Append the A matrix and B vector
    A = [           A,             minusones;
	     zeros(numsoft,nvars),  -eye(numsoft)];

    B = [         B;
	     zeros(numsoft,1)];
	
    % Append the H matrix and f vector
    
    weights = weights(:,soft);
    if ~isempty(weights), % JMM MOD,22.12.99
      if strcmp(normtype,'inf-norm')
	  weights = weights(:,1);
      end
	
      [nrows,ncols] = size(H);

      H = [            H,          zeros(nrows,numsoft); 
	     zeros(numsoft,ncols),   diag(weights(2,:))];
      f = [      f;
	     weights(1,:)'];
    end % of ~isempty(weights), % JMM MOD,22.12.99

	% set upper and lower bounds on variables. Note that epsilon is
	% bounded above by 1e6
    if ~isempty(vlb) & ~isempty(vub)	 
      vlb = vlb(:);
      vub = vub(:);
      vlb = [vlb; zeros(numsoft,1)];
      vub = [vub; 1e6*ones(numsoft,1)];
    end
  
    if ~isempty(X)
      X = [         X; 
	       zeros(numsoft,1)];	 
    end
  end
else
  weights=[];
end

% run new QP problem

qpsolver = 'quadprog'; % Change this to enable QUADPROG.M
%qpsolver = 'qpsoft';

if strcmp(qpsolver,'qp')

  [X,lambda,how]=qp(H,f,A,B,vlb,vub,X,[],verbosity); 

elseif strcmp(qpsolver,'quadprog'), % Now runs using QUADPROG.M, instead of QP.M

  if verbosity <= 0
    options=optimset('Display','off');
  else
	options=optimset('Display','iter');
  end

  [X,fval,exitflag,output,lambda]=quadprog(H,f,A,B,[],[],vlb,vub,X,options);

  if exitflag > 0
    how ='ok';
  else
    how ='infeasible';
  end

  lambda = [lambda.ineqlin; lambda.lower; lambda.upper];

else

  error(['QP solver ',qpsolver,' not supported.']) % JMM MOD,21.12.99

end % if

% return vector of optimized variables and vector of constraint violations.
% MOD ECK 10/3/99

epsilon = [];
if ~isempty(weights)
  epsilon = zeros(ncstr,1); 
  for i = 1:numsoft
    epsilon(soft(i)) = X(nvars+i);
  end
  temp=lambda;
  lambda = zeros(ncstr,1);
  for i = 1:numkeep
	lambda(keep(i)) = temp(i);
  end
end
X = X(1:nvars);
