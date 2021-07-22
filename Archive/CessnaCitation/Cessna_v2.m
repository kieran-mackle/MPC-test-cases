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
Hp = 10;

% Define control horizon
Hu = 3;

% Total time horizon (?)
N = 50;

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

% Define initial state
pitch = 0*deg;
altitude = 0;
alt_rate = 0;
elevator = -10*deg;


% Define plant
plant = ss(A, B, C, D);
plant = c2d(plant, Ts);


% Create full simulation vectors
Cc = repmat({C}, 1, Hp+1);
CC = blkdiag(Cc{:});
Dc = repmat({D}, 1, Hp+1);
DD = blkdiag(Dc{:});
rc = repmat({r'}, 1, Hp+1);
RR = cell2mat(rc)';
Qc = repmat({Q}, 1, Hp+1);
QQ = blkdiag(Qc{:});


% Simulate Trajectory

opt_x = zeros(n, N);
opt_y = zeros(p, N);
opt_u = zeros(m, N);

xk = [ pitch, altitude, alt_rate, elevator ]';

for k = 1:N
    
    % Run cvx 
    cvx_begin
        variable x(n, Hp+1)
        variable u(m, Hp+1)

        XX = reshape(x, [n*(Hp+1),1]);
        UU = reshape(u, [m*(Hp+1),1]);
        YY = CC*XX + DD*UU;

        minimize( (YY-RR)' * QQ * (YY-RR) )
        subject to
            x(:, 1) == xk;

        for i = 1:Hu-1
            subject to
                x(:, i+1) == A * x(:, i) + B * u(:, i);
        end

        for i = Hu:Hp
            subject to
                u(:, i) == u(:, i-1);
                x(:, i+1) == A * x(:, i) + B * u(:, i);
        end

    cvx_end
    
    
    % Extract latest results
    x_out = reshape(x, [n*(Hp+1),1]);
    u_out = reshape(u, [m*(Hp+1),1]);
    y_out = CC*x_out + DD*u_out;
    y_out = reshape(y_out, [p,Hp+1]);

    elevator = u(:,1);
    pitch    = y_out(1, 1);
    altitude = y_out(2, 1);
    alt_rate = y_out(3, 1);
    
    % Update state
    xk = [ pitch, altitude, alt_rate, elevator ]';
    
    % Construct output
    opt_y(:,k) = [ pitch, altitude, alt_rate ]';
    opt_u(:,k) = elevator;
    
end

% Visualise results
% t = 0:Ts:Hp*Ts;
% opt_x = reshape(x, [n*(Hp+1),1]);
% opt_y = CC*opt_x;
% opt_y = reshape(opt_y, [p,Hp+1]);

% elevator = opt_u(:,:);
% pitch    = opt_y(1,:);
% altitude = opt_y(2,:);
% alt_rate = opt_y(3,:);
% 
% 
% figure(1);
% clf;
% set(gcf,'color','w');
% 
% subplot(2,2,1);
% hold on;
% grid on;
% title("Pitch angle");
% plot(t, pitch*rad);
% xlabel("Time (s)");
% ylabel("Pitch (deg)");
% 
% subplot(2,2,2);
% hold on;
% grid on;
% title("Altitude");
% plot(t, altitude);
% xlabel("Time (s)");
% ylabel("Altitude (m)");
% 
% subplot(2,2,3);
% hold on;
% grid on;
% title("Altitude rate");
% plot(t, alt_rate);
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



