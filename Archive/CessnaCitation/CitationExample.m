% ----------------------------------------------------------------------- %
% Cessna Citation Aircraft Example from Predictive Control (Maciejowskie) %
% ----------------------------------------------------------------------- %


% Define open-loop aircraft dynamics model
A = [-1.2822    0.0     0.98    0.0;
     0.0        0.0     1.0     0.0;
     -5.4293    0.0   -1.8366   0.0;
     -128.2    128.2    0.0     0.0];

B = [-0.3;
      0.0;
     -17.0;
      0.0];
  
C = [   0     1     0   0;
        0     0     0   1;
     -128.2  128.2  0   0];

D = [0.0;
     0.0;
     0.0];

% Define Plant
plant   = ss(A, B, C, D);
x0      = zeros(4,1);

% Manipulated variables
MV = struct('Min',{20,-15,-30,20},'Max',{15,30,20});

% Output variables
OV = struct('Min',{[-0.5;-Inf],[-100;-Inf]},...
            'Max',{[0.5;Inf],[100;Inf]},...
            'ScaleFactor',{1,200});

% Controller tuning weights
Weights = struct('MV',[0 0],'MVRate',[0.1 0.1],'OV',[10 10]);

% MPC controller
Ts      = 0.5;            % Sample time
p       = 10;             % Prediction horizon
m       = 3;              % Control horizon
mpcobj  = mpc(plant, Ts, p, m, Weights, MV, OV);






