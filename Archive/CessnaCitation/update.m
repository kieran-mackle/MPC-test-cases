% ----------------------------------------------------------------------- %
% Cessna Citation Aircraft Example from Predictive Control (Maciejowskie) %
% ----------------------------------------------------------------------- %



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


% model is of the form x(k+1) = Ax(k) + Bu(k)


% Cost function
% V(k) = sum abs(z_cap(k+i) - r(k+i))^2 + sum abs(u_cap(k+i))^2 



% MPC Update
y_k_1   = A*x_k + B*u_k + f_0;
y_k     = C*x_k + D*u_k + d;



























