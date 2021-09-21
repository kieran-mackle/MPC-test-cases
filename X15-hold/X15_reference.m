function output = X15_reference(mpc_input)

% fpa = mpc_input.initial.state(11);
% thr = mpc_input.initial.state(12);
% 
% if fpa > 40*pi/180
%     % Flight path angle too high
%     mpc_input.constraints.hard.input(1,2) = 0;
% elseif fpa < -40*pi/180
%     mpc_input.constraints.hard.input(1,1) = 0;
% else
%     mpc_input.constraints.hard.input(1,:) = [-10, 10]*pi/180;
% end
% 
% if thr > 1
%     % Flight path angle too high
%     mpc_input.constraints.hard.input(2,2) = 0;
% elseif thr < 0.1
%     mpc_input.constraints.hard.input(2,1) = 0;
% else
%     mpc_input.constraints.hard.input(1,:) = [-0.2,   0.2];
% end


mpc_input.reference = [20e3, 0, 0, 6, 0]';

output = mpc_input;

