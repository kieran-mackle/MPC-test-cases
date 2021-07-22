% test script to make big matrices
clearvars


% Define prediction horizon
P = 10;



C = [   0     1     0   0;
        0     0     0   1;
     -128.2  128.2  0   0];
xk = [ 0 * pi / 180, 0, 0, -10 * pi / 180 ]'; % Initial state


p = size(C, 1);     % Number of outputs

r = [0, 40, 0]';    % Reference state
Q = eye(p);         % State weightings


Cc = repmat({C}, 1, P);
CC = blkdiag(Cc{:});
xkc = repmat({xk'}, 1, P);
XX = cell2mat(xkc)';






% CC  = zeros(P*size(C,2), P*size(C,1));
% for i = 1:P
%     CC(i:i+2, i:i+3) = C
% end




% for i = 1:p
%     CC(:,:,i) = C;
%     XX(:,:,i) = xk;
%     YY(:,:,i) = CC(:,:,i) * XX(:,:,i);
%     RR(:,:,i) = r;
%     QQ(:,:,i) = Q;
%     
%     VV(:,:,i) = (YY(:,:,i)-RR(:,:,1))' * QQ(:,:,i) * (YY(:,:,i)-RR(:,:,1));
% end





