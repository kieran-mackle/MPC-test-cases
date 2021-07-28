function output = multiphase_ball_dynamics(input)

for i = 1:2
    X       = input.phase(i).state;
%     U       = input.phase(i).control;
    Fdot    = input.phase(i).control;
    m       = 10; 

    s       = X(:,1);
    v       = X(:,2);

%     F       = U(:,1);
    F       = X(:,3);
    a       = F./m;

    s_dot   = v;
    v_dot   = a;

    ref_errorQ2 = 50*(5 - s) .* (5 - s);
%     controller_movement = [U(1); U(2:end) - U(1:end-1)] .* ...
%                           [U(1); U(2:end) - U(1:end-1)];
    controller_movement = 0*Fdot .* Fdot;
    
%     output(i).dynamics = [s_dot, v_dot];
    output(i).dynamics = [s_dot, v_dot, Fdot];
    output(i).integrand = [ ref_errorQ2 + controller_movement];

end

end