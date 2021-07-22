function [y,u,ym,x,xm]=scmpcnl2(pmod,imod,ywt,uwt,blocks,p,tvec, ...
                        setpts,ulim,ylim,Kest,ydist,mdist, ...
                        umdist,udist,x0,u0,hint,xm0,suwt,sywt,normtype)
%SCMPCNL2 simulates MPC controller with nonlinear Simulink plant and soft constraints.
%
% [y,u,ym,x,xm]=scmpcnl2(pmod,imod,ywt,uwt,M,P,tvec, ...
%                        r,ulim,ylim,Kest,z,d,w,wu, ...
%                        x0,u0,hint,xm0,suwt,sywt,normtype)
%                     
% This function has been modified from the MPC Toolbox function SCMPCNL for use 
% with the book 'Predictive Control with Constraints' by J.M.Maciejowski,
% (Prentice Hall, 2001) - see Appendix C of the book.
%
% SCMPCNL2 designs an MPC-type controller for constrained problems.
% Simulates the closed-loop system with hard and/or soft constraints
% (inputs and outputs) using quadratic programming.
% This version uses a SIMULINK plant model that can include
% nonlinear dynamics, etc.
%
% The initial values of the estimator can now be set on the
% command line. XM0 is used as the estimator's initial states.  The
% states of the estimator are also passed as outputs to variable XM, in
% the augmented form (see MPCAUGSS.M) (ECK 27/4/99).
% Soft constraint functionality added. Requires QPSOFT.M (ECK 03/03/99)
% Some Simulink plants might be time-dependent and the ability to
% specify the start time of the simulation has been included (ECK
% 12/1/99).
% Function changed to use the QP.M routine in the Optimization Toolbox, instead
% of the DANTZGMP.M routine in the MPC Toolbox (ECK 11/12/98).
% This version uses `sim(...)' to simulate a Simulink MDL model.
% Received from Larry Ricker 20 may 1998.
% Also this version breaks cleanly from simulation if QP infeasibility 
% encountered. (JMM 20.5.98)
%
% ---------usage------------------
% Inputs:
%  pmod:     Plant model in SIMULINK format
%  imod:     Internal model in MOD format -- used as basis for 
%            controller design.
%  ywt:      Penalty weighting for setpoint tracking.
%  uwt:      Penalty weighting for changes in manipulated variables.
%  M:        Number of moves OR sequence of blocking factors.
%  P:        Length of prediction horizon.
%  tvec:     One of [TEND] or [T0 TEND] where T0 and TEND are
%            the start and finishing times of the simulation. Default T0=0.
%  r:        Setpoint sequence, one column for each output.
%  ulim:     [Ulow Uhigh delU].  NOTE:  delU must be finite.
%  ylim:     [Ylow Yhigh].
%  Kest:     Estimator gain matrix.
%  z:        Measurement noise.
%  d:        Measured disturbance (for feedforward control).
%  w:        General unmeasured disturbance.
%  wu:       Unmeasured disturbance added to manipulated variables.
%  x0:       Initial states of plant (optional).
%  u0:       Initial KNOWN plant inputs u(-1) (optional).  Must include 
%            manipulated variables and measured disturbances, if any.
%  hint:     Fixed step size for interation of nonlinear ODEs (optional).
%  xm0:      Initial states of the estimator in augmented state space
%            form (see MPCAUGSS.M) (optional). Default is XM0=0.
%  suwt:     Penalty weights for constraint violations in manipulated
%            variables (optional). Takes the same format as ulim for 1-norm and
%            2-norm. Use the form [SUWT1 SUWT2] for mixed-norm. Scalar
%            for inf-norm. Empty argument or 'inf' indicates that the
%            constraint(s) are to be kept hard. A 0 indicates removal of
%            the constraint(s).  
%  sywt:     Penalty weights for constraint violations in output
%            variables (optional). Takes the same format as ylim for 1-norm and
%            2-norm. Use the form [SYWT1 SYWT2] for mixed-norm. Empty
%            argument or 'inf' indicates that the constraint(s) are to be kept
%            hard. A 0 indicates removal of the constraint(s).
%  normtype: Norm type to be used in penalising constraint
%            violations (optional). One of '1-norm', '2-norm',
%            'mixed-norm' or 'inf-norm'.
% Outputs:   y (system response), u (manipulated variables),
%            ym (model response), x (state trajectory of plant),
%            xm (state trajectory of internal model in augmented form).
%
% See also PLOTALL, PLOTEACH, SCMPC, SMPCCL, SMPCCON, SMPCEST, SMPCSIM.
%
% DISCLAIMER: This function has been modified from the original function 
%             SCMPCNL of the MPC Toolbox, with permission. The MathWorks, Inc.
%             cannot provide any support for this modified function.

% Copyright (c) 1994-98 by The MathWorks, Inc.
% $Revision: 1.6 $% Copyright (c) 1999 by Eric C. Kerrigan
% $ Revision 2.03 $ $ Date: 01/12/1999 16:50:03 $
% Changed by J.M.Maciejowski, 06.06.2001.

% +++ Check input arguments for errors and set default values. +++

if nargin == 0
   disp('USAGE:  [y,u,ym,x,xm]=scmpcnl2(pmod,imod,ywt,uwt,M,P,tvec,r,ulim,ylim,...')
   disp('                         Kest,z,d,w,wu,x0,u0,hint,xm0,suwt,sywt,normtype)')
   return
elseif nargin < 9
   error('Too few input arguments')
end

if isempty(pmod) | isempty (imod)
   error('Plant and internal model must be supplied.')
end

% Get internal model specifications.

[phii,gami,ci,di,minfoi]=mod2ss(imod);
tsamp=minfoi(1);
ni=minfoi(2);
nui=minfoi(3);
nvi=minfoi(4);
mi=nui+nvi;
nwi=minfoi(5);
if nwi > 0             % Strip off unmeasured disturbance inputs if they exist.
   gami=gami(:,1:mi);
   di=di(:,1:mi);
end
nymi=minfoi(6);
nyui=minfoi(7);
nyi=nymi+nyui;

% Check for errors and inconsistencies in the models.

if any(any(di(:,1:nui)))
   error(['IMOD:  first nu=',int2str(nui),' columns of D must be zero.'])
end

[inp,outp,xp0]=chknl(pmod,minfoi);
nwp=inp-mi;

if isempty(p)
   p=1;
elseif p < 1
   error('Specified prediction horizon is less than 1')
end

if isempty(ywt)
   ywt=ones(1,nyi);
   nywt=1;
else
   [nywt,ncol]=size(ywt);
   if ncol ~= nyi | nywt <= 0
      error('YWT is wrong size')
   end
   if any(any(ywt < 0))
      error('One or more elements of YWT are negative')
   end
end

if isempty(uwt),
   uwt=zeros(1,nui);
   nuwt=1;
else
   [nuwt,ncol]=size(uwt);
   if ncol ~= nui | nuwt <= 0
      error('UWT is wrong size')
   end
   if any(any(uwt < 0))
      error('UWT is negative')
   end
end

if isempty(setpts)
   nset=1;
   setpts=zeros(1,nyi);
else
   [nset,ncol]=size(setpts);
   if ncol ~= nyi
      error('Setpoint input matrix has incorrect dimensions')
   end
end

if isempty(blocks)
   blocks=ones(1,p);
   nb=p;
else
   [nrow,nb]=size(blocks);
   if nrow ~= 1 | nb < 1 | nb > p
      error('M vector is wrong size')
   end
   if any(blocks < 1)
      error('M contains an element that is < 1')
   end

   if nb == 1

%  This section interprets "blocks" as a number of moves, each
%  of one sampling period duration.

      if blocks > p
         disp('WARNING: M > P.  Truncated.')
         nb=p;
      elseif blocks <= 0
         disp('WARNING: M <= 0.  Set = 1.')
         nb=1;
      else
         nb=blocks;
      end
      blocks=[ones(1,nb-1) p-nb+1];

   else

% This section interprets "blocks" as a vector of blocking factors.

      sumblocks=sum(blocks);
      if sumblocks > p
               disp('WARNING:  sum(M) > P.')
               disp('          Moves will be truncated at P.')
               nb=find(cumsum(blocks) > p);
               nb=nb(1);
               blocks=blocks(1,1:nb);
      elseif sumblocks < p
         nb=nb+1;
         blocks(nb)=p-sumblocks;
         disp('WARNING:  sum(M) < P.  Will extend to P.')
      end
   end
end

% MOD ECK 13/1/99 - Sets initial time t0 
% MOD ECK 10/3/99 - Changed variable name from TSTART to T0
if length(tvec) == 2
  t0 = tvec(1);
  tend = tvec(2);
  if t0 > tend
	   error('Start time T0 > Finish time TEND');
  end
else
  if length(tvec) == 1
    t0=0;
    tend=tvec;
  else
	error('TVEC has incorrect size: has to be of form [TEND] or [T0 TEND]');
  end
end
clear tvec % MOD ECK 10/3/99
% End of MOD 13/1/99

% Check the constraint specifications.  First set up some indices to pick out
% certain columns of the ulim and ylim matrices.

iumin=[1:nui];     % Points to columns of ulim containing umin.
iumax=iumin+nui;   % Points to columns of ulim containing umax.
idumax=iumax+nui;  % Points to columns of ulim containing delta u max.
iymin=[1:nyi];     % Points to columns of ylim containing ymin.
iymax=iymin+nyi;   % Points to columns of ylim containing ymax.

% Now check the values supplied by the user for consistency.

[nulim,ncol]=size(ulim);
if ncol ~= 3*nui | nulim <= 0
   error('ULIM matrix is empty or wrong size.')
elseif any(any(ulim(:,idumax) < 0))
   error('A constraint on DELTA U was < 0')
elseif any(any(ulim(:,iumax)-ulim(:,iumin) < 0))
   error('A lower bound on U was greater than its upper bound')
end

% When using the DANTZGMP routine for the QP problem, we must have all
% bounds on delta u finite.  A bound that is finite but large can cause
% numerical problems.  The following loop checks for this.

ichk=0;
for i=idumax
   ifound=find(ulim(:,i) > 1e6);
   if ~ isempty(ifound)
      ichk=1;
      ulim(ifound,i)=1e6*ones(length(ifound),1);
   end
end
if ichk      
   disp('One or more constraints on delta_u were > 1e6.')
   disp('Reduced to 1e6 to prevent numerical problems in QP.')
end

if nargin > 9
   if isempty(ylim)
      ylim=[-inf*ones(1,nyi) inf*ones(1,nyi)];
   else
      [nylim,ncol]=size(ylim);
      if ncol ~= 2*nyi | nylim <= 0
         error('YLIM matrix is wrong size')
      elseif any(any(ylim(:,iymax)-ylim(:,iymin) < 0))
         error('A lower bound on y was greater than its upper bound')
      end
   end
else
   ylim=[-inf*ones(1,nyi) inf*ones(1,nyi)];
end

if nargin > 10
   if isempty(Kest)
      Kest=[zeros(ni,nymi)
                eye(nymi)
           zeros(nyui,nymi)];
   else
      [nrow,ncol]=size(Kest);
      if nrow ~= ni+nyi | ncol ~= nymi
         error('Estimator gain matrix is wrong size')
      end
   end
else
   Kest=[zeros(ni,nymi)
            eye(nymi)
         zeros(nyui,nymi)];
end

% +++ If there are unmeasured outputs, add columns of zeros
%     to the estimator gain matrix.

if nyui > 0
   Kest=[Kest zeros(ni+nyi,nyui)];
end

if nargin > 11
   if isempty(ydist)
      nyd=1;
      ydist=zeros(1,nymi);
   else
      [nyd,ncol]=size(ydist);
      if ncol ~= nymi
         error(['Z (measurement noise) matrix must have ',num2str(nymi),...
                ' columns.'])
      end
   end
else
   nyd=1;
   ydist=zeros(1,nymi);
end

if nargin > 12
   if ~isempty(mdist)
      [nmd,ncol]=size(mdist);
      if nvi == 0
         error('V input was supplied, but model has no measured disturbances.')
      elseif ncol ~= nvi
         error(['Measured disturbance (V) matrix must have ',num2str(nvi),...
                ' columns.'])
      end
   else
      nmd=1;
      mdist=zeros(1,nvi);
   end
else
   nmd=1;
   mdist=zeros(1,nvi);
end

if nargin > 13
   if ~isempty(umdist)
      [numd,ncol]=size(umdist);
      if nwp == 0 & ncol > 0
         error('W input was specified, but PMOD does not include model for it!')
      elseif ncol ~= nwp
         error(['Unmeasured disturbance (W) matrix must have ',num2str(nwp),...
                ' columns.'])
      end
   else
      numd=1;
      umdist=zeros(1,nwp);
   end
else
   numd=1;
   umdist=zeros(1,nwp);
end

if nargin > 14
   if isempty(udist)
      nud=1;
      udist=zeros(1,nui);
   else
      [nud,ncol]=size(udist);
      if ncol ~= nui
         error(['WU disturbance matrix must have ',num2str(nui),...
                ' columns.'])
      end
   end
else
   nud=1;
   udist=zeros(1,nui);
end

if nargin > 15
   if ~isempty(x0)
      [rows,cols]=size(x0);
      if cols ~= 1 | rows ~= length(xp0)
         error(['x0 must be a vector, ',num2str(length(xp0)),' by 1.'])
      end
      xp0=x0;
   end
end

if nargin > 16
   if ~ isempty(u0)
      [rows,cols]=size(u0);
      if cols ~= 1 | rows ~= nui+nvi
         error(['u0 must be a vector, ',num2str(nui+nvi),' by 1.'])
      end
      manvold=u0(1:nui,1);
      vold=u0(nui+1:nui+nvi,1);
   else
      manvold=zeros(nui,1);
      vold=zeros(nvi,1);
      disp('++ WARNING:  initial values of MVs and MDs taken as zero.')
   end
else
   manvold=zeros(nui,1);
   vold=zeros(nvi,1);
   disp('++ WARNING:  initial values of MVs and MDs taken as zero.')
end
deltav=[];

y0=[];

if nargin > 17
   if ~ isempty(hint)
      [rows,cols]=size(hint);
      if rows ~= 1 | cols ~= 1
         error('HINT must be a scalar')
      end
   else
      hint=tsamp/100;       % default choice of integration step size
   end
else
   hint=tsamp/100;       % default choice of integration step size
end

% Start of MOD ECK 27/4/99
if nargin > 18
  if ~isempty(xm0)
      [rows,cols]=size(xm0);
      if cols > 1 | rows ~= ni+nyi
         error(['XM0 must be in augmented form,',num2str(ni+nyi),' by 1.'])
	  end
	else
	  xm0=zeros(ni+nyi,1);
	  %disp('WARNING: Internal model initial value XM0 = 0')
	end
else
   xm0=zeros(ni+nyi,1);
   %disp('WARNING: Internal model initial value XM0 = 0')
 end
% End of MOD ECK 27/4/99

% Soft constraints ECK 02/03/99
if nargin == 20 | nargin == 21
  error('Not enough input arguments - NORMTYPE must be specified as well.')
end

if nargin > 21
  if isempty(normtype) & (~isempty(suwt) | ~isempty(sywt))
	error('NORMTYPE must be specified.')
  elseif ~isempty(normtype) & isempty(suwt) & isempty(sywt) 
	error('SUWT and/or SYWT must be specified.')
  end
  if ~isempty(normtype) & (~isempty(suwt) | ~isempty(sywt))
	
	[nrowsuwt,ncolsuwt]=size(suwt);
	[nrowsywt,ncolsywt]=size(sywt);
    
	% Determine if correct NORMTYPE has been specified
	if strcmp(normtype,'inf-norm')
	  if ncolsuwt > 1 | nrowsuwt > 1 | ~isempty(sywt)
	    error('For inf-norm, SUWT must be a scalar and SYWT must be empty.')
	  end
	  if suwt < 0
		error('SUWT must be positive')
	  end
	  softconstraints=2; % indicates that inf-norm is being used
	  sweights=suwt;
	elseif ~(strcmp(normtype,'1-norm') | strcmp(normtype,'2-norm') | strcmp(normtype,'mixed-norm'))
	   error('NORMTYPE must be 1-norm, 2-norm, mixed-norm or inf-norm.')
	else
	   softconstraints=1; % indicates 1-norm, 2-norm or mixed-norm
	end
    
	if softconstraints==1 % if NORMTYPE is 1-norm, 2-norm or mixed-norm
	  % Check to see if elements of SUWT and SYWT are non-negative
	  if any(any(suwt < 0))
        error('One or more elements of SUWT are negative')
      end
	  if any(any(sywt < 0))
        error('One or more elements of SYWT are negative')
      end
	
	  % Check to see if SUWT and SYWT has the correct size. If one of them
	  % is empty, then keep corresponding constraints hard.
	  if isempty(suwt)
		suwt = inf*ones(1,3*nui); 
		[nrowsuwt,ncolsuwt]=size(suwt);
	  end
	  if isempty(sywt)
		sywt = inf*ones(1,2*nyi);
		[nrowsywt,ncolsywt]=size(sywt);
	  end
	  
	  if ~(ncolsuwt == 3*nui | ncolsuwt == 2*3*nui)
        error('SUWT is the wrong size.')
      end
	  if ~(ncolsywt == 2*nyi | ncolsywt == 2*2*nyi)
        error('SYWT is the wrong size.')
      end
	  
	  if ncolsuwt == 2*3*nui & (strcmp(normtype,'1-norm') | strcmp(normtype,'2-norm'))
	    error(['SUWT is only allowed to have ', num2str(3*nui),' columns for 1-norm and 2-norm.'])
	  end 
      if ncolsywt == 2*2*nyi & (strcmp(normtype,'1-norm') | strcmp(normtype,'2-norm'))
	    error(['SYWT is only allowed to have ', num2str(2*nyi),' columns for 1-norm and 2-norm.'])
	  end
	  
	  if strcmp(normtype,'mixed-norm') & ncolsuwt == 3*nui
	    disp('WARNING: Same SUWT is being used for both 1-norm and 2-norm.')
	    suwt = [suwt suwt];
		[nrowsuwt,ncolsuwt]=size(suwt);
	  end
	  if strcmp(normtype,'mixed-norm') & ncolsywt == 2*nyi
	    disp('WARNING: Same SYWT is being used for both 1-norm and 2-norm.')
	    sywt = [sywt sywt];
		[nrowsywt,ncolsywt]=size(sywt);
	  end
	end
	
  else
	normtype =[];
    sweights =[];
	softconstraints=0;
  end
else
  normtype =[];
  sweights =[];
  softconstraints=0;
end
% End of soft constraints ECK 02/03/99

if nargin > 22 % MOD ECK 02/03/99
   error('Too many input arguments')
end

% ++++ Beginning of controller design calculations. ++++

% The following index vectors are used to pick out certain columns
% or rows in the state-space matrices.

iu=[1:nui];          % columns of gami, di related to delta u. % MOD ECK 13/1/99
iv=[nui+1:nui+nvi];  % points to columns for meas. dist. in gamma.
iym=[1:nymi];       % index of the measured outputs.

% +++ Augment the internal model state with the outputs.

[PHI,GAM,C,D,N]=mpcaugss(phii,gami,ci,di);

% +++ Calculate the basic projection matrices +++

pny=nyi*p;          % Total # of rows in the final projection matrices.
mnu=nb*nui;             % Total number of columns in final Su matrix.

Cphi=C*PHI;
Sx=[      Cphi
   zeros(pny-nyi,N)];
Su=[    C*GAM(:,iu)
    zeros(pny-nyi,nui)];
if nvi > 0
   Sv0=[ C*GAM(:,iv)
        zeros(pny-nyi,nvi)    ];
else
   Sv0=[]; 
end

r1=nyi+1;
r2=2*nyi;
for i=2:p
   if nvi > 0
      Sv0(r1:r2,:)=Cphi*GAM(:,iv);
   end
   Su(r1:r2,:)=Cphi*GAM(:,iu);
   Cphi=Cphi*PHI;
   Sx(r1:r2,:)=Cphi;
   r1=r1+nyi;
   r2=r2+nyi;
end

Sdel=eye(nui); % Sdel is to be a block-lower-triangular matrix in which each
               % block is an identity matrix.  Used in constraint definition.
eyep=eye(nyi); % eyep is a matrix containing P identity matrices (dimension nyi)
               % stacked one on top of the other.
for i=2:p
   eyep=[eyep;eye(nyi)];
end
for i=2:nb
   Sdel=[Sdel;eye(nui)];
end

% If number of moves > 1, fill the remaining columns of Su and Sdel, 
% doing "blocking" at the same time.

if nb > 1
   k = nui;
   blocks=cumsum(blocks);
   for i = 2:nb
      row0=blocks(i-1)*nyi;
      row1=(i-1)*nui;
      Su(row0+1:pny,k+1:k+nui)=Su(1:pny-row0,1:nui);
      Sdel(row1+1:mnu,k+1:k+nui)=Sdel(1:mnu-row1,1:nui);
      k=k+nui;
   end
end

% Set up weighting matrix on outputs.  Q is a column vector
% containing the diagonal elements of the weighting matrix, SQUARED.

irow=0;
for i=1:p
   Q(irow+1:irow+nyi,1)=ywt(min(i,nywt),:)';
   irow=irow+nyi;
end
Q=Q.*Q;

% Set up weighting matrix on manipulated variables.  R
% is a column vector containing the diagonal elements, SQUARED.

uwt=uwt+10*sqrt(eps);  %for numerical stability
irow=0;
for i=1:nb
   R(irow+1:irow+nui,1)=uwt(min(i,nuwt),:)';
   irow=irow+nui;
end
R=R.*R;

% Usually, some of the general inequality constraints are not used.
% This section sets up index vectors for each type of constraint to
% pick out the ones that are actually needed for the problem.  This
% helps to minimize the size of the QP.

% First set up column vectors containing the bounds for each type of
% constraint over the entire prediction horizon.  For the inputs, the
% resulting vectors must be length mnu.  For outputs, length is pny.

umin=ulim(:,iumin)';
umin=umin(:);         % Stetches the matrix out into one long column
umax=ulim(:,iumax)';
umax=umax(:);
dumax=ulim(:,idumax)';
dumax=dumax(:);
ymin=ylim(:,iymin)';
ymin=ymin(:);
ymax=ylim(:,iymax)';
ymax=ymax(:);
% Soft constraints ECK 11/2/99
if softconstraints == 1
  temp=suwt(:,iumin)';
  uminswt=temp(:);
  temp=suwt(:,iumax)';
  umaxswt=temp(:);
  temp=suwt(:,idumax)';
  dumaxswt=temp(:);
  temp=sywt(:,iymin)';
  yminswt=temp(:);
  temp=sywt(:,iymax)';
  ymaxswt=temp(:);
  if ncolsuwt == 2*3*nui
    temp=suwt(:,3*nui+iumin)';
    uminswt=[uminswt, temp(:)];
    temp=suwt(:,3*nui+iumax)';
    umaxswt=[umaxswt, temp(:)];
    temp=suwt(:,3*nui+idumax)';
    dumaxswt=[dumaxswt, temp(:)];
  end
  if ncolsywt == 2*2*nyi
    temp=sywt(:,2*nyi+iymin)';
    yminswt=[yminswt, temp(:)];
    temp=sywt(:,2*nyi+iymax)';
    ymaxswt=[ymaxswt, temp(:)];
  end
  
  % Work out which constraints in ULIM are relaxed/removed
  if ncolsuwt == 3*nui
	suwt = [suwt suwt];
  end
  if nulim > nrowsuwt
    suwt = [suwt; ones(nulim-nrowsuwt,1)*suwt(end,:)];
  end
  if nrowsuwt > nulim
    ulim = [ulim; ones(nrowsuwt-nulim,1)*ulim(end,:)];
  end
  weights1 = suwt(:,1:3*nui)';
  weights2 = suwt(:,3*nui+(1:3*nui))';
  weights = [weights1(:) weights2(:)]'; % in order to work with find and any
  ulim = ulim'; 
  ulim = [ulim(:) ulim(:)]';
  
  soft = find(any(weights~=0) & ~any(weights==inf) & any(ulim~=inf) & any(ulim~=-inf));
  numsoft = length(soft); 
  hard = find(any(weights==inf) & any(ulim~=inf)  & any(ulim~=-inf));
  numhard = length(hard);
  remove = find(~any(weights~=0) & any(ulim~=inf) & any(ulim~=-inf) );
  numremove = length(remove);
  
  % output number of hard and soft constraints in ULIM
  if numhard ~= 0
	disp('The following constraint(s) in ULIM are hard:')
	disp(hard)
  end
  if numsoft ~= 0
    disp('The following constraint(s) in ULIM are soft:')
	disp(soft)
  end
  if numremove ~= 0
	disp('The following constraint(s) in ULIM will be removed:')
	disp(remove)
  end

  % Work out which constraints in YLIM are relaxed/removed
  if ncolsywt == 2*nyi
	sywt = [sywt sywt];
  end
  if nylim > nrowsywt
    sywt = [sywt; ones(nylim-nrowsywt,1)*sywt(end,:)];
  end
  if nrowsywt > nylim
    ylim = [ylim; ones(nrowsywt-nylim,1)*ylim(end,:)];
  end
  weights1 = sywt(:,1:2*nyi)';
  weights2 = sywt(:,2*nyi+(1:2*nyi))';
  weights = [weights1(:) weights2(:)]'; % in order to work with find and any
  ylim = ylim'; 
  ylim = [ylim(:) ylim(:)]';
  
  soft = find(any(weights~=0) & ~any(weights==inf) & any(ylim~=inf) & any(ylim~=-inf));
  numsoft = length(soft); 
  hard = find(any(weights==inf) & any(ylim~=inf)  & any(ylim~=-inf));
  numhard = length(hard);
  remove = find(~any(weights~=0) & any(ylim~=inf) & any(ylim~=-inf) );
  numremove = length(remove);
  
  % output number of hard and soft constraints in YLIM
  if numhard ~= 0
	disp('The following constraint(s) in YLIM are hard:')
	disp(hard)
  end
  if numsoft ~= 0
    disp('The following constraint(s) in YLIM are soft:')
	disp(soft)
  end
  if numremove ~= 0
	disp('The following constraint(s) in YLIM will be removed:')
	disp(remove)
  end
  
  clear temp weights hard soft remove numhard numsoft numremove weights1 weights2 suwt sywt
end
% ECK 11/2/99
clear ulim ylim       % Releases memory no longer needed.

lenu=length(umin);
if lenu > mnu         % Has user specified more bounds than necessary?
   disp('WARNING:  too many rows in ULIM matrix.')
   disp('          Extra rows deleted.')
   umin=umin(1:mnu);
   umax=umax(1:mnu);
   dumax=dumax(1:mnu);
elseif lenu < mnu     % If fewer rows than needed, must copy last one.
   r2=[lenu-nui+1:lenu];
   for i=1:round((mnu-lenu)/nui)
      umin=[umin;umin(r2,:)];
      umax=[umax;umax(r2,:)];
      dumax=[dumax;dumax(r2,:)];
   end
end

leny=length(ymin);
if leny > pny         % Has user specified more bounds than necessary?
   disp('WARNING:  too many rows in YLIM matrix.')
   disp('          Extra rows deleted.')
   ymin=ymin(1:pny);
   ymax=ymax(1:pny);
elseif leny < pny     % If fewer rows than needed, must copy last one.
   r2=[leny-nyi+1:leny];
   for i=1:round((pny-leny)/nyi)
      ymin=[ymin;ymin(r2,:)];
      ymax=[ymax;ymax(r2,:)];
   end
end

if softconstraints == 1 % Soft constraints ECK 12/2/99
  lenuswt=length(uminswt);
  if lenuswt > mnu         % Has user specified more weights than necessary?
     disp('WARNING:  too many rows in U part of SWT matrix.')
     disp('          Extra rows deleted.')
     uminswt=uminswt(1:mnu,:);
     umaxswt=umaxswt(1:mnu,:);
     dumaxswt=dumaxswt(1:mnu,:);
  elseif lenuswt < mnu     % If fewer rows than needed, must copy last one.
     r2=[lenuswt-nui+1:lenuswt];
     for i=1:round((mnu-lenuswt)/nui)
       uminswt=[uminswt;uminswt(r2,:)];
       umaxswt=[umaxswt;umaxswt(r2,:)];
       dumaxswt=[dumaxswt;dumaxswt(r2,:)];
     end
  end

  lenyswt=length(yminswt);
  if leny > pny         % Has user specified more weights than necessary?
     disp('WARNING:  too many rows in Y part of SWT matrix.')
     disp('          Extra rows deleted.')
     yminswt=yminswt(1:mnu,:);
     ymaxswt=ymaxswt(1:mnu,:);
  elseif lenyswt < pny     % If fewer rows than needed, must copy last one.
     r2=[lenyswt-nyi+1:lenyswt];
     for i=1:round((pny-lenyswt)/nyi)
       yminswt=[yminswt;yminswt(r2,:)];
       ymaxswt=[ymaxswt;ymaxswt(r2,:)];
     end
   end
end
% ECK 12/2/99

% The bounds on delta u must always be included in the problem.  The
% other bounds should only be included as constraints if they're finite.
% Generate vectors that contain a list of the finite constraints.

iumin=find(umin ~= -inf);
iumax=find(umax ~=  inf);
iymin=find(ymin ~= -inf);
iymax=find(ymax ~=  inf);

% Delete the infinite values.  At the same time, form the coefficient
% matrix for the inequality constraints.  Do this by picking out only
% the equations actually needed according to the lists established above.
% Finally, calculate the constant part of the RHS of the inequality
% constraints for these equations.

A=eye(mnu);        % These are the equations that are always present.
rhscon=2*dumax;    % They are the bounds on delta u.  A is the coefficient
                   % matrix and rhscon is the constant part of the RHS.
				   
if softconstraints == 1 % Soft constraints ECK 11/2/99
  sweights = dumaxswt; 
end

if ~ isempty(iumin)    % Add equations for lower bound on u
   umin=umin(iumin);
   A=[A;-Sdel(iumin,:)];
   rhscon=[rhscon;-Sdel(iumin,:)*dumax-umin];
   if softconstraints == 1 % Soft constraints ECK 11/2/99
     sweights=[sweights;uminswt(iumin,:)]; 
   end
else
   umin=[];
end
if ~ isempty(iumax)    % Add equations for upper bound on u
   umax=umax(iumax);
   A=[A;Sdel(iumax,:)];
   rhscon=[rhscon;Sdel(iumax,:)*dumax+umax];
   if softconstraints == 1 % Soft constraints ECK 11/2/99
     sweights=[sweights;umaxswt(iumax,:)];
   end
else
   umax=[];
end
if ~ isempty(iymin)    % Add equations for lower bound on y
   ymin=ymin(iymin);
   A=[A;-Su(iymin,:)];
   rhscon=[rhscon;-Su(iymin,:)*dumax-ymin];
   if softconstraints == 1 % Soft constraints ECK 11/2/99
     sweights=[sweights;yminswt(iymin,:)];
   end
else
   ymin=[];
end
if ~ isempty(iymax)    % Add equations for upper bound on y
   ymax=ymax(iymax);
   A=[A;Su(iymax,:)];
   rhscon=[rhscon;Su(iymax,:)*dumax+ymax];
   if softconstraints == 1 % Soft constraints ECK 11/2/99
     sweights=[sweights;ymaxswt(iymax,:)];
   end
else
   ymax=[];
end

[nc,dumdum]=size(A);   % Save total number of inequality constraints.

% +++ Define the matrices needed for the QP +++


SuTQ=Su'*diag(Q);
B=SuTQ*Su+diag(R);
clear Su
a=B'*dumax;   % This is a constant term that adds to the initial basis
              % in each QP
% The following 3 lines must be included for the dantzgmp routine ECK 11/12/98
% B=inv(B);    
% TAB=[-B   B*A' ;A*B   -A*B*A'];
% clear A B

%  ++++ SIMULATION SECTION ++++

% Initialization of states, etc.

%xi=zeros(ni+nyi,1);   % States of the augmented internal model.
xi=xm0;                % MOD ECK 27/4/99
%if (xm0 == 0) & (xmp == 0)  % MOD ECK 15/1/99
%  xi = zeros(ni+nyi,1);   
%else
%  xi=[xm0-xmp 
%      ci*xm0];              % MOD ECK 11/12/98
%end
%xmsum=xmp;                  % MOD ECK 12/1/99 States will be updated
                            % later on
							
up=[];

if nvi > 0
   v=mdist(1,:)';
   up=[up;v];             % Is this necessary? up is updated a few lines
                          % further down ECK 11/12/98
else
   v=[];
end
if nwp > 0
   w=umdist(1,:)';
else
   w=[];
end

% t0=0;                          % Time at start of simulation.
                                 % Simulink model might be time-dependent ECK 11/12/98
up=[manvold+udist(1,:)';v;w];
if exist([pmod '.mdl']) == 4 % MOD ECK 03/03/99
   % Get initial noise-free plant outputs using Simulink
   plntissl=1;
	%tvec=[0 0]';   % MOD ECK 11/12/98
	tvec=[t0 t0]';  % MOD
	simopts=simset('OutputVariables','ty',...
		'MaxStep',tsamp,...
		'OutputPoints','specified',...
		'InitialState',xp0);
	[T,X,Y]=sim(pmod,tvec,simopts,[tvec [up';up']]);
	yp=Y(end,:)';  %re-order outputs
elseif exist(pmod) % MOD ECK 03/03/99
   plntissl=0;
	yp=feval(pmod,t0,xp0,up,3);      % Initial outputs of the plant.  Uses
                                 	% SIMULINK system call with FLAG=3.
else
   error('PMOD must be an M-file or MDL-file.') % MOD ECK 03/03/99
end

IKC=eye(ni+nyi)-Kest*C;
  
% Simulation

Hwait=waitbar(0,'Simulation in progress...');
nstep=max(fix((tend-t0)/tsamp)+1,2);  % Number of sampling periods in
                                      % simulation. MOD ECK 12/1/99
for i=1:nstep
   waitbar(i/nstep);

   ypnew=yp;

   ypnew(1:nymi,1)=yp(1:nymi,1)+ydist(min(i,nyd),:)';   % add measurement noise.
   setpt=setpts(min(i,nset),:)';         % current setpoints

% Calculate starting basis vector for the QP

   xi=IKC*xi+Kest*ypnew;  % measurement update for state estimator.
   y0=Sx*xi;
   if nvi > 0
      v=mdist(min(i,nmd),:)';            % current measured disturbances.
      deltav=v-vold;
      vold=v;
      y0=y0 + Sv0*deltav;
   end

   rhsa=a+SuTQ*(eyep*setpt-y0);

% Update the RHS of the inequality constraints

   rhsc=zeros(mnu,1);           
   del=Sdel(:,1:nui)*manvold;      % vector of previous value of manip. vars.

   if ~ isempty(iumin)    % Equations for lower bound on u
      rhsc=[rhsc;del(iumin,:)];
   end
   if ~ isempty(iumax)    % Equations for upper bound on u
      rhsc=[rhsc;-del(iumax,:)];
   end
   if ~ isempty(iymin)    % Equations for lower bound on y
      rhsc=[rhsc;y0(iymin,:)];
   end
   if ~ isempty(iymax)    % Equations for upper bound on y
      rhsc=[rhsc;-y0(iymax,:)];
   end

   rhsc=rhsc+rhscon;      % Add on the constant part computed earlier.

% Set up and solve the QP;

% The following block uses DANTZGMP commented out by ECK 11/12/98
%   basisi=[    -TAB(1:mnu,1:mnu)*rhsa
%           rhsc-TAB(mnu+1:mnu+nc,1:mnu)*rhsa];
%   ibi=-[1:mnu+nc]';
%   ili=-ibi;
%   [basis,ib,il,iter]=dantzgmp(TAB,basisi,ibi,ili);
%   if iter < 0
%      disp('Infeasible QP.  Check constraints.');          % MODS by JMM 20.5.98
%      disp(['Step = ',int2str(i),',  Time =',num2str(t0)]) % MODS                
%      disp('Simulation terminating prematurely!')          % MODS       
%      break                                                % MODS
%    end
%   deltau=[];
%   for j=1:nui
%      if il(j) <= 0
%         deltau(j,1)=-dumax(j,1);
%      else
%        deltau(j,1)=basis(il(j))-dumax(j,1);
%      end
%   end

% Start of QP using QPSOFT.M routine ECK 11/2/99

   H=(B+B')/2; % This is to ensure that the Hessian is symmetric 
               % and that QPSOFT.M does not return the corresponding warning message
               % because of truncation error. Be careful that the correct
               % non-inverted B matrix be used here, not as used by DANTZGMP.
   
   [deltau,lambda,how,epsilon]=qpsoft(H,-rhsa,A,rhsc,zeros(mnu,1),2*dumax,[],sweights,normtype,0);
   if ~strcmp(how,'ok')
      disp(['Problem with QP. The message is: ', how])
      disp('Abandoning simulation')              
      break                                      
   end
   deltau=deltau-dumax;
   deltau=deltau(1:nui);

% End of QP using QPSOFT.M routine ECK 11/2/99

   manvold=deltau+manvold;
   
% Input and state updates

   ui=[deltau;deltav];  
   ud=udist(min(i,nud),:)';
   if nwp > 0
      w=umdist(min(i,numd),:)';
   end
   up=[manvold+ud;v;w];
   
   y(i,:)=yp';
   if nargout > 1
      u(i,:)=manvold';
   end
   if nargout > 2
      ym(i,:)=xi(ni+1:ni+nyi,:)';
   end
   if nargout > 3
      x(i,:)=xp0';
	end
   if nargout > 4                % MOD ECK 27/4/99 
	 xm(i,:)=xi';             
   end                           % End of MOD

   if i == nstep
      break       % Loop exits from here.  
   end            % Prevents extra integration of model equations.
  
   xi=PHI*xi+GAM*ui;          % State update for state estimator.
   
   if ~plntissl
      % Plant model is not in Simulink format.
		% Use EULERFS function to integrate the plant equations for one
		% sampling period.  Assumes ALL inputs are constant over this time
		% period.  If you have SIMULINK, substitute RK45 or equivalent
		% for the EULERFS function for faster execution and better
		% reliability.

  		tvec=[t0 t0+tsamp];         % Starting and ending time.
  		xp0=eulerfs(pmod,tvec,xp0,up,hint);
  		t0=t0+tsamp;
      yp=feval(pmod,t0,xp0,up,3); % New plant output (noise free);
   else
		% Plant is a Simulink model in MDL format. ECK 11/12/98	
		% Use SIMULINK to integrate the plant model for one sampling
		% period.  All simulation options (such as integration method)
		% come from the Simulink model, except overrides in SIMOPTS.

   	tvec=[t0 t0+tsamp]';        % Starting and ending time.
  	 	simopts=simset('OutputVariables','ty',...
			'MaxStep',tsamp,...
			'OutputPoints','specified',...
			'InitialState',xp0,'FinalStateName','xp0');
   	[T,X,Y]=sim(pmod,tvec,simopts,[tvec [up';up']]);
   	xp0=xp0(:);
   	t0=t0+tsamp;
   	yp=Y(end,:)'; % New plant output (noise free);
	end
end
close(Hwait)

% +++ End of SCMPCNL +++
