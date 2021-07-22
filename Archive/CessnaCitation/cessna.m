% ----------------------------------------------------------------------- %
% Cessna Citation Aircraft Example from Predictive Control (Maciejowskie) %
% ----------------------------------------------------------------------- %
clearvars;
deg = pi/180;
rad = 180/pi;

% Time-constant
Tref = 5;

% Sampling interval
Ts = 0.5;

% Define prediction horizon
Hp = 20;

% Define control horizon
Hu = 10;


% Define open-loop dynamics model
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

n = size(A, 1);     % Number of states
m = size(B, 2);     % Number of control inputs
p = size(C, 1);     % Number of outputs

r = [0, 40, 0]';    % Set-point
Q = eye(p);         % State weightings
R = eye(m);         % Control weightings

xk = [ 0*deg, 0, 0, 0 ]'; % Initial state


% Define plant
plant = ss(A, B, C, D);
plant = c2d(plant, Ts);
% 
% Ad = plant.A;
% Bd = plant.B;
% Cd = plant.C;
% Dd = plant.D;

Ad = plant.A;
Bd = plant.B;
Cd = plant.C;
Dd = plant.D;

% Create full simulation vectors
Cc = repmat({Cd}, 1, Hp+1);
CC = blkdiag(Cc{:});
Dc = repmat({Dd}, 1, Hp+1);
DD = blkdiag(Dc{:});
rc = repmat({r'}, 1, Hp+1);
RR = cell2mat(rc)';
Qc = repmat({Q}, 1, Hp+1);
QQ = blkdiag(Qc{:});
big_Rc = repmat({R}, 1, Hp);
big_R = blkdiag(big_Rc{:});


for i = 1:Hp+1
    delta(i,i:i+1) = [1,-1];
end

% What are the states?
% How to impose state constraints?
% How to do forward sim

% Run cvx 
cvx_begin
    variable x(n, Hp+1)
    variable u(m*(Hp+1), 1)
    
    XX = reshape(x, [n*(Hp+1),1]);
    UU = reshape(u, [m*(Hp+1),1]);
    YY = CC*XX + DD*UU;
    
    du = u(2:end) - u(1:end-1);
    dx = x(2:end) - x(1:end-1);
    
    minimize( (YY-RR)' * QQ * (YY-RR) + du'*big_R*du)
    subject to
        x(:, 1) == xk;
        -15*deg <= u(:) <= 15*deg % elevator angle
        -30*deg*Ts <= du(:) <= 30*deg*Ts % Elevator slew rate
        -30*deg*Ts <= u(2) - u(1) <= 30*deg*Ts % Elevator slew rate
        
    for i = 1:Hp
        subject to
            x(:, i+1) == Ad * x(:, i) + Bd * u(m*(i-1)+1:m*i);
            -20*deg<= YY(p*(i-1)+1,1) <= 20*deg  % pitch angle
            
    end
    
cvx_end


% Visualise results
t = 0:Ts:Hp*Ts;
opt_x = reshape(x, [n*(Hp+1),1]);
opt_y = CC*opt_x;
opt_y = reshape(opt_y, [p,Hp+1]);

elevator = u;
pitch    = opt_y(1,:);
altitude = opt_y(2,:);
alt_rate = opt_y(3,:);

figure(1);
clf;
set(gcf,'color','w');

subplot(2,2,1);
hold on;
grid on;
title("Pitch angle");
plot(t, pitch*rad);
xlabel("Time (s)");
ylabel("Pitch (deg)");

subplot(2,2,2);
hold on;
grid on;
title("Altitude");
plot(t, altitude);
xlabel("Time (s)");
ylabel("Altitude (m)");

subplot(2,2,3);
hold on;
grid on;
title("Altitude rate");
plot(t, alt_rate);
xlabel("Time (s)");
ylabel("Altitude rate (m/s)");

subplot(2,2,4);
hold on;
grid on;
title("Elevator angle");
stairs(t, elevator*rad);
xlabel("Time (s)");
ylabel("Elevator angle (deg)");


% Forward simulate using control output
fs = forward_sim(Ad,Bd,Cd,Dd,xk,u);
fs_pitch = fs(1,:);
fs_alt = fs(2,:);
fs_alt_rate = fs(3,:);

figure(2);
clf;
set(gcf,'color','w');

subplot(3,1,1);
hold on;
grid on;
title("Pitch angle");
plot(t, fs_pitch*rad);
xlabel("Time (s)");
ylabel("Pitch (deg)");

subplot(3,1,2);
hold on;
grid on;
title("Altitude");
plot(t, fs_alt);
xlabel("Time (s)");
ylabel("Altitude (m)");

subplot(3,1,3);
hold on;
grid on;
title("Altitude rate");
plot(t, fs_alt_rate);
xlabel("Time (s)");
ylabel("Altitude rate (m/s)");






