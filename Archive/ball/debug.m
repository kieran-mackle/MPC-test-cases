% Extract delta-u terms
du_k        = dUbar(1);
du_kp1      = dUbar(2);
du_kp2      = dUbar(3);

% Predicted state at k+1
xbar_predicted(1:2);
ubar_k      = ubar_km1 + du_k;
xbar_kp1    = Ad*xbar_k + Bd*ubar_k + f0d;

% Predicted state at k+2
xbar_predicted(3:4);
ubar_kp1    = ubar_km1 + du_k + du_kp1;
xbar_kp2    = Ad*xbar_kp1 + Bd*ubar_kp1 + f0d;

% Predicted state at k+3
xbar_predicted(5:6);
ubar_kp2    = ubar_km1 + du_k + du_kp1 + du_kp2;
xbar_kp3    = Ad*xbar_kp2 + Bd*ubar_kp2 + f0d;


mpc_pred = xbar_predicted; %[xbar_predicted(1:2), xbar_predicted(3:4),xbar_predicted(5:6)]';
manual = [xbar_kp1; xbar_kp2; xbar_kp3];

[mpc_pred, manual]

% h=2;
% [mpc_pred(n*(h-1)+1:n*(h-1)+n,:), manual(n*(h-1)+1:n*(h-1)+n,:), ...
%     mpc_pred(n*(h-1)+1:n*(h-1)+n,:) + (h-1)*Ad*f0d]