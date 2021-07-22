%% Assignment 2: Run each in section and then run simulink model!

%Back up parameters
Ms = 311.905644;
Ks = 456.7803532;
Bs = 22.95995478;
Mus = 59.91490365;
Kus = 50574.47491;


%%  Discretise SS Model
A = [0 0 1 0; 0 0 0 1; -(Ks+Kus)/Mus Ks/Mus -Bs/Mus Bs/Mus; Ks/Ms -Ks/Ms Bs/Ms -Bs/Ms];
B = [0 0; 0 0; -1/Mus Kus/Mus; 1/Ms 0];
C = [1 0 0 0; -1 1 0 0; 0 0 1 0; 0 0 0 1];
D = [0 -1; 0 0; 0 0; 0 0];
T = 0.05;

[sysd, G] = c2d(ss(A,B,C,D), T, 'zoh');

Aq = sysd.A;
Bq = sysd.B;
Cq = sysd.C;
Dq = sysd.D;

%% Now augment the SS Model
Ce = [0 0 0 -1];
A_aug = [Aq zeros(4,1); Ce 1];
B_aug = [Bq; zeros(1, 2)];
B1 = B_aug(:,1);
B2 = B_aug(:,2);
C_aug = [C zeros(4, 1); zeros(1,4) 1];
D_aug = [D; zeros(1, 2)];
D2 = D_aug(:,2);

% Cost matrices 
delta = 1;
gamma = 5;
Q = [1 0 0 0 0; 0 delta 0 0 0; 0 0 1 0 0; 0 0 0 gamma 0; 0 0 0 0 gamma];
P = Q;
r = 1e-8;

% Constraint matrices
Mi = [0 1 0 0 0; 0 -1 0 0 0; 0 0 0 1 0; 0 0 0 -1 0];
Gi = [0.2; 0.2; 0.4; 0.4];
Ei = [1; -1];
Hi = [2500; 2500];


%% High horizon time

% Time step
N = 20;
R = r*eye(N);
Omega = zeros(5*N, 5*N);
Lambda = zeros(5*N, 5*N);
Phi = zeros(5*N, 5);
Gamma_D = zeros(5*N,N);
Gamma_F = zeros(5*N, N);
Gamma_Z = zeros(5*N, N);

M = zeros(4*N, 5*N);
G = zeros(4*N, 1);
H = zeros(2*N, 1);
E = zeros(2*N, N);
for i = 1:N
    Omega(5*(i-1)+1:5*i, 5*(i-1)+1: 5*i) = Q; 
    Lambda(5*(i-1)+1:5*i, 5*(i-1)+1: 5*i) = C_aug; 
    Phi(5*(i-1)+1:5*i, :) = A_aug^i;
    Gamma_D(5*(i-1)+1:5*i, i) = D2;
    M(4*(i-1)+1:4*i, 5*(i-1)+1:5*i) = Mi;
    G(4*(i-1)+1:4*i) = Gi;
    H(2*(i-1)+1:2*i) = Hi;
    E(2*(i-1)+1:2*i, i) = Ei;
    for j = 1:i
        Gamma_F(5*(i-1)+1:5*i, j) = A_aug^(i - j)*B1;
        Gamma_Z(5*(i-1)+1:5*i, j) = A_aug^(i - j)*B2;
    end
    Gamma_F(5*(i-1)+1: 5*i, i+1: N) = zeros(5, (N-i));
    Gamma_Z(5*(i-1)+1: 5*i, i+1: N) = zeros(5, (N-i));
end

