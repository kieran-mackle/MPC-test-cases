% Aircraft with unstable poles
% ----------------------------


% Define open-loop aircraft dynamics model
A = [-0.0151 -60.5651 0 -32.174;
     -0.0001 -1.3411 0.9929 0;
      0.00018 43.2541 -0.86939 0;
      0      0       1      0];
  
B = [-2.516 -13.136;
     -0.1689 -0.2514;
     -17.251 -1.5766;
      0        0];
  
C = [0 1 0 0;
     0 0 0 1];
 
D = [0 0;
     0 0];
 
% Define Plant
plant = ss(A,B,C,D);
x0 = zeros(4,1);

% The open-loop response is unstable, as shown by the poles:
% pole(plant)


% Manipulated variables
MV = struct('Min',{-25,-25},'Max',{25,25});


% Output variables
OV = struct('Min',{[-0.5;-Inf],[-100;-Inf]},...
            'Max',{[0.5;Inf],[100;Inf]},...
            'ScaleFactor',{1,200});

% Controller tuning weights
Weights = struct('MV',[0 0],'MVRate',[0.1 0.1],'OV',[10 10]);


% MPC controller
Ts = 0.05;          % Sample time
p = 10;             % Prediction horizon
m = 2;              % Control horizon
mpcobj = mpc(plant,Ts,p,m,Weights,MV,OV);


% Simulate in Simulink
% mdl = 'mpc_aircraft';
% open_system(mdl)
% sim(mdl)
sim('mpc_aircraft_example')











