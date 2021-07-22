% ----------------------------------------------------------------------- %
% Cessna Citation Aircraft Example from Predictive Control (Maciejowskie) %
% ----------------------------------------------------------------------- %
% Inclusion of time varying reference

clearvars;


deg = pi/180;
rad = 180/pi;

% Time-constant
Tref = 5; 

% Sampling interval
Ts = 0.5;

% Define prediction horizon
Hp = 10;

% Define control horizon
Hu = 3;

% Total time horizon (?)
total_time = 20;


% Define open-loop dynamics model (continuous)
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

Q = eye(p);         % State weightings
R = eye(m);         % Control weightings


% Construct reference
% r = [0, 40, 0]';    % Set-point for period 1
% r_2 = [0, 20, 0]';    % Set-point for period 2


% Define plant and convert to discrete-time domain
plant = ss(A, B, C, D);
plant = c2d(plant, Ts);

Ad = plant.A;
Bd = plant.B;
Cd = plant.C;
Dd = plant.D;

% Create full simulation vectors
Cc = repmat({Cd}, 1, Hp+1);
CC = blkdiag(Cc{:});
Dc = repmat({Dd}, 1, Hp+1);
DD = blkdiag(Dc{:});
Qc = repmat({Q}, 1, Hp+1);
QQ = blkdiag(Qc{:});
big_Rc = repmat({R}, 1, Hp);
big_R = blkdiag(big_Rc{:});


for i = 1:Hp+1
    delta(i,i:i+1) = [1,-1];
end


% Simulate Trajectory
N = total_time/Ts;
t = zeros(1,N);
opt_x = zeros(n, N);
opt_y = zeros(p, N);
opt_u = zeros(m, N);
set_points = zeros(p, N);

x_0 = zeros(n,1); % Initial state
u_0 = zeros(m,1);

% f = waitbar(0);
x_k = x_0;
u_k = u_0;
r = [0, 40, 0]';

for k = 1:N
    t(k) = Ts*(k);
    
    % Update waitbar
%     waitbar(k/(N+1), f, sprintf('Solving MPC Problem (%d%%)', floor(k/(N+1)*100)));
    
    % Calculate reference vector
    
    if t(k) < 2000
        r = [0, 40, 0]';
    else
        r = [0, 20, 0]';
    end
    
    set_points(:,k) = r;
    rc = repmat({r'}, 1, Hp+1);
    RR = cell2mat(rc)';
    
    % Add results to vectors
    opt_u(:,k) = u_k;
    opt_x(:,k) = x_k;
    opt_y(:,k) = Cd*x_k + Dd*u_k;
    
    % Run cvx
    cvx_begin quiet
        variable x(n, Hp+1)
        variable u(m*(Hp+1), 1)

        XX = reshape(x, [n*(Hp+1),1]);
        UU = u;
        YY = CC*XX + DD*UU;

        du = u(m+1:end) - u(1:end-m);
        
        minimize( (YY-RR)' * QQ * (YY-RR) + du'*big_R*du)
        subject to
            x(:, 1) == x_k;
            -15*deg <= u(:) <= 15*deg % elevator angle
            -30*deg*Ts <= du(:) <= 30*deg*Ts % Elevator slew rate
            -30*deg*Ts <= u(m+1:2*m) - u(1:m) <= 30*deg*Ts % Initial elevator slew rate

            
        for i = 1:Hu
            subject to
                x(:, i+1) == Ad * x(:, i) + Bd * u(m*(i-1)+1:m*i);
                -20*deg<= YY(p*(i-1)+1,1) <= 20*deg  % pitch angle
        end
        
        for i = Hu+1:Hp
            subject to
                u(m*(i-1)+1:m*i) == u(m*(i-2)+1:m*(i-1));
                x(:, i+1) == Ad * x(:, i) + Bd * u(m*(i-1)+1:m*i);
                -20*deg<= YY(p*(i-1)+1,1) <= 20*deg  % pitch angle
        end
    cvx_end
    
    % Extract latest results
    x_k1 = x(:,1);
    u_k1 = u(1:m);
    
    % Update state
    x_k = Ad*x_k1 + Bd*u_k1;
    u_k = u_k1;
    
end
% close(f)

elevator = opt_u;
pitch    = opt_y(1,:);
altitude = opt_y(2,:);
alt_rate = opt_y(3,:);


figure(1);
clf;
set(gcf,'color','w');
sgtitle("MPC Control Output");

subplot(2,2,1);
hold on;
grid on;
title("Pitch angle");
plot(t, pitch*rad);
stairs(t, set_points(1,:), 'k--');
xlabel("Time (s)");
ylabel("Pitch (deg)");

subplot(2,2,2);
hold on;
grid on;
title("Altitude");
plot(t, altitude);
stairs(t, set_points(2,:), 'k--');
xlabel("Time (s)");
ylabel("Altitude (m)");

subplot(2,2,3);
hold on;
grid on;
title("Altitude rate");
plot(t, alt_rate);
stairs(t, set_points(3,:), 'k--');
xlabel("Time (s)");
ylabel("Altitude rate (m/s)");

subplot(2,2,4);
hold on;
grid on;
title("Elevator angle");
stairs(t, elevator*rad);
xlabel("Time (s)");
ylabel("Elevator angle (deg)");



% % Forward simulate using control output
% fs = forward_sim(Ad,Bd,Cd,Dd,x_0,elevator);
% fs_pitch = fs(1,:);
% fs_alt = fs(2,:);
% fs_alt_rate = fs(3,:);
% 
% figure(2);
% clf;
% set(gcf,'color','w');
% sgtitle("Forward Simulation Output");
% 
% subplot(2,2,1);
% hold on;
% grid on;
% title("Pitch angle");
% plot(t, fs_pitch*rad);
% xlabel("Time (s)");
% ylabel("Pitch (deg)");
% 
% subplot(2,2,2);
% hold on;
% grid on;
% title("Altitude");
% plot(t, fs_alt);
% stairs(t, set_points(2,:), 'k--');
% xlabel("Time (s)");
% ylabel("Altitude (m)");
% 
% subplot(2,2,3);
% hold on;
% grid on;
% title("Altitude rate");
% plot(t, fs_alt_rate);
% xlabel("Time (s)");
% ylabel("Altitude rate (m/s)");
% 
% subplot(2,2,4);
% hold on;
% grid on;
% title("Elevator angle");
% stairs(t, elevator*rad);
% xlabel("Time (s)");
% ylabel("Elevator angle (deg)");

