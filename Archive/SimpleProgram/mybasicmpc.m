% ===================================================================== %
%                              BASIC MPC                                %
% ===================================================================== %

% Time-constant
Tref = 6;

% Sampling interval
if Tref == 0
    Ts = 1;
else
    Ts = Tref/10;
end


% Define plant ---------------------------------- %
% Change below for new plant
nump = 1;
denp = [1,-1.4, 0.45];
plant = tf(nump, denp, Ts);
% Change above for new plant model. 

plant = tf(plant);
nump = get(plant, 'num'); nump = nump{:};
denp = get(plant, 'den'); denp = denp{:};

nnump = length(nump) - 1;
ndenp = length(denp) - 1;

if nump(1) ~= 0
    error("Plant must be strictly proper");
end

if any(abs(roots(denp)) > 1)
    disp("Warning: Unstable plant")
end


% Define model
% Change below for new model
model = plant;
% Change above for new model
model = tf(model);
numm = get(model, 'num'); numm = numm{:};
denm = get(plant, 'den'); denm = denm{:};
nnumm = length(numm) - 1;
ndenm = length(denm) - 1;

if numm(1) ~= 0
    error("Model must be strictly proper");
end

if any(abs(roots(denm)) > 1)
    disp("Warning: Unstable model")
end

nump = [zeros(1, ndenp - nnump -1), nump];
numm = [zeros(1, ndenm - nnumm -1), numm];

% Define prediction horizon
if Tref == 0
    P = 5;
else
    P = round(0.8*Tref/Ts);
end

% Define control horizon
M = 1;

% Compute model step response values
stepresp = step(model, [0:Ts:max(P)*Ts]);
theta = zeros(length(P), M);
for j = 1:length(P)
    theta(j,:) = [stepresp(P(j):-1:max(P(j) - M + 1, 1))', zeros(1, M-P(j))];
end

S = stepresp(P);

% Compute reference error factor at coincidence points
if Tref == 0
    errfac = zeros(length(P), 1);
else
    errfac = exp(-P*Ts/Tref);
end


% Simulation parameters
if Tref == 0
    tend = 100*Ts;
else
    tend = 10*Tref;
end

nsteps = floor(tend/Ts);
tvec = (0:nsteps-1)' * Ts;

% Define set-point vector
setpoint = ones(nsteps+max(P), 1);

% Define vectors to hold input and output signals
uu = zeros(nsteps,1); 
yp = zeros(nsteps,1); 
ym = zeros(nsteps,1);

% Initial conditions
umpast = zeros(ndenm,1); 
uppast = zeros(ndenp,1);
ympast = zeros(ndenm,1); 
yppast = zeros(ndenp,1);



% Simulation
for k = 1:nsteps
   % Define reference trajectory at coincidence points:
  errornow = setpoint(k)-yp(k);
  reftraj = setpoint(k+P) - errornow*errfac;

  % Free response of model over prediction horizon:
  yfpast = ympast;
  ufpast = umpast;
  for kk=1:max(P) % Prediction horizon
    ymfree(kk) = numm(2:nnumm+1)*ufpast - denm(2:ndenm+1)*yfpast;
    yfpast=[ymfree(kk); yfpast(1:length(yfpast)-1)];
	ufpast=[ufpast(1); ufpast(1:length(ufpast)-1)];
  end
  
  % Compute input signal uu(k):
  if k>1
    dutraj = theta\(reftraj-ymfree(P)');
    uu(k) = dutraj(1) + uu(k-1);  
  else
	dutraj = theta\(reftraj-ymfree(P)');
	uu(k) = dutraj(1) + umpast(1);
  end
  
  % Simulate plant:
  % Update past plant inputs
  uppast = [uu(k); uppast(1:length(uppast)-1)];
  yp(k+1) = -denp(2:ndenp+1)*yppast + nump(2:nnump+1)*uppast; % Simulation
  % Update past plant outputs
  yppast = [yp(k+1);yppast(1:length(yppast)-1)];


  % Simulate model:
  % Update past model inputs
  umpast = [uu(k); umpast(1:length(umpast)-1)];
  ym(k+1) = -denm(2:ndenm+1)*ympast + numm(2:nnumm+1)*umpast;  % Simulation
  % Update past model outputs
  ympast = [ym(k+1);ympast(1:length(ympast)-1)];

end

disp('***** Results from script file BASICMPC:');
disp(['Tref = ', num2str(Tref), ', Ts = ', num2str(Ts), ...
      ' P = ', int2str(P), ' (steps), M = ', int2str(M)]);
diffpm = get(plant-model, 'num');

if diffpm{:} == 0
    disp('Model = Plant');
else
    disp('Plant-model mismatch');
end


figure(1);
subplot(211);
plot(tvec, yp(1:nsteps), '-', tvec, setpoint(1:nsteps), '--');
grid;
title('Plant output (solid) and set-point (dashed)')
xlabel('Time')

subplot(212);
stairs(tvec, uu, '-');
grid;
title('Input');
xlabel('Time');







