function output = X15_outputs(input)
    latd    = input(1);
    dist    = input(3);
    uD      = input(4);
    vD      = input(5);
    wD      = input(6);
    fda     = input(11);
    thr     = input(12);

    % Calculate altitude
    % ------------------
    a       = 6378137.0;
    f       = 3.33528106e-3;
    Re      = a.*(1-(f/2).*(1-cos(2.*latd)) + ...
              (5.*f.^2./16).*(1-cos(4.*latd)));
    h       = -(dist + Re);

    % Calculate Mach number
    % ---------------------
    gamma           = 1.4;
    R               = 2.8705e+02;
    vBED            = [uD,vD,wD];
    V               = sqrt(sum(vBED.^2,2));
    [Temp,~,~]      = GetAtmo(h);
    sos             = sqrt(gamma.*R.*Temp);
    Ma              = V./sos;
    
    % Calculate flight path angle
    % ---------------------------
    delta           = GetDelta(h, latd);
    T_DG            = TM_DG(delta);
    vBEG            = T_DG' * vBED';
    [V,~,fpa]       = car2pol(vBEG);
    
    output = [h, fda, thr, Ma, fpa];
end