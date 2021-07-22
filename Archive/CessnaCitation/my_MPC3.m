% ----------------------------------------------------------------------- %
% Cessna Citation Aircraft Example from Predictive Control (Maciejowskie) %
% ----------------------------------------------------------------------- %
clearvars;

% Then wrap into function
% Then use on X-15

deg = pi/180;
rad = 180/pi;

% Sampling interval
Ts = 0.5;

% Define prediction horizon
Hp = 10;

% Define control horizon
Hu = 3;

% Total time horizon (?)
N = 20;


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

[opt_u, opt_y, t] = solve_MPC(A,B,C,D,Hp,Hu,Ts, 1e-4);

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


% % Forward simulate using control output
% fs = forward_sim(Ad,Bd,Cd,Dd,x_0,elevator);
% fs_pitch = fs(1,:);
% fs_alt = fs(2,:);
% fs_alt_rate = fs(3,:);
% 
% figure(2);
% % clf;
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