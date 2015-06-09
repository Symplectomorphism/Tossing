function [t,x] = tossing()

t0 = 0; tf = 0.5;
x0 = [0;0;0.03; zeros(3,1)]; % States for object translation
x0 = [x0; reshape(eye(3), 9, 1); zeros(9,1)]; % States for object orientation
x0 = [x0; 0.1; 0]; % States for the radial deformation
x0 = [x0; 0.0; 0]; % States for the vertical deformation
x0 = [x0; zeros(6,1)]; % States for hand translation
x0 = [x0; reshape(eye(3), 9, 1); zeros(9,1)]; % States for hand orientation
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
par.mo = 0.5;           % [kg] 
par.mh = 1.0;           % [kg]
% par.r0 = 0.3048;   % [m]
par.nu = 0.35;       
par.E = 5e+06;       % [Pa] = [kg/m/s/s]
par.rho = 1.2e+03;   % [kg/m/m/m]
par.gamma = 9.81;        % [m/s/s]
par.I = diag([1/4*par.mo*x0(25)^2, 1/4*par.mo*x0(25)^2, 1/2*par.mo*x0(25)^2]);    % [kg m^2]

[t,y] = ode45(@EoM, [t0, tf], x0, options, par);

% Interpolate for fixed time steps
tstep = 0.01;
tspan = 0:tstep:tf;
x = interp1(t, y, tspan, 'spline');
t = tspan';
lent = length(t);

theta = zeros(lent,1);
for i = 1:lent
    temp = reshape(x(i,7:15), 3, 3)';
    theta(i) = atan2(temp(2,1), temp(1,1));
end

figure(1), clf
subplot(2,1,1)
plot(t, x(:,25))
xlabel('t')
ylabel('r [m]')
title('Dough radius')
subplot(2,1,2)
plot(t, x(:,27))
xlabel('t')
ylabel('h [m]')
title('Dough height')

figure(2), clf
subplot(2,2,4)
plot(t, theta)
xlabel('t')
ylabel('\theta', 'Interpreter', 'TeX')
title('Yaw rotation')
subplot(2,2,1)
plot(t, x(:,23))
xlabel('t')
ylabel('\omega_x', 'Interpreter', 'TeX')
subplot(2,2,2)
plot(t, x(:,18))
xlabel('t')
ylabel('\omega_y', 'Interpreter', 'TeX')
subplot(2,2,3)
plot(t, x(:,19))
xlabel('t')
ylabel('\omega_z', 'Interpreter', 'TeX')

    function dx = EoM(t,x,par)
        %% Parameters and states
        
        mo = par.mo;
        mh = par.mh;
        gamma = par.gamma;
        nu = par.nu;
        E = par.E;
        rho = par.rho;
        
        % Orientation and angular velocity of the object
        Rpo = reshape(x(7:15), 3, 3)';
        wpohat = reshape(x(16:24), 3, 3)';
        Rpod = Rpo*wpohat;
        gpo = [Rpo, x(1:3); zeros(1,3), 1];
        
        % Inertia matrix with radius and height as variables
        Io = mo/18*diag( [3*x(25)^2 + x(27)^2, 3*x(25)^2 + x(27)^2, 1/6*x(25)^2] );
        Iop = mo/18*diag( [6*x(25)*x(26) + 2*x(27)*x(28), 6*x(25)*x(26) + ...
                        2*x(27)*x(28), 1/3*x(25)*x(26)] );
        % Inertia matrix with radius as the only variable
%         Ip = m*x(25)/2*x(26)*[1,0,0;0,1,0;0,0,2];

        Ioinv = inv(Io);
        
        % Orientation and the angular velocity of the hand
        Rph = reshape(x(35:43), 3, 3)';
        wphhat = reshape(x(44:52), 3, 3)';
        Rphd = Rph*wphhat;
        gph = [Rph, x(29:31); zeros(1,3), 1];
%         gph_Inverse = [Rph', -Rph'*x(29:31); zeros(1,3), 1];
        gph_Inverse = homInv(gph);
        
        % Inertia matrix of the hand
        Ih = diag([1.5, 1, 0.5])*1e-02;
        
        
        % Kinematics
        Roc = eye(3); poc = zeros(3,1);
        goc = [Roc, poc; zeros(1,3), 1];
%         goc_Inverse = [Roc', -(Roc')*poc; zeros(1,3), 1];
        goc_Inverse = homInv(goc);
%         Adgoc = [Roc, skewmat(poc)*Roc; zeros(3,3), Roc];
%         Adgoc_Inverse = [Roc', -(Roc')*skewmat(poc); zeros(3,3), Roc'];
        Adgoc = Ad(goc);
        Adgoc_Inverse = Ad(goc_Inverse);
        
        
        ghc = gph_Inverse*gpo*goc;
        Rhc = ghc(1:3, 1:3);
        phc = ghc(1:3,4);
%         ghc_Inverse = [Rhc', -(Rhc')*phc; zeros(1,3), 1];
        ghc_Inverse = homInv(ghc);
%         Adghc = [Rhc, skewmat(phc)*Rhc; zeros(3,3), Rhc];
%         Adghc_Inverse = [Rhc', -(Rhc')*skewmat(phc); zeros(3,3), Rhc'];
        Adghc = Ad(ghc);
        Adghc_Inverse = Ad(ghc_Inverse);

        % Mass matrix
        M = blkdiag(mo*eye(3), Io, mh*eye(3), Ih);
        
        % Natural forces
        phio = zeros(6,1);
        phio(1:3) = cross(skewvec(wpohat), mo*x(4:6)) + mo*gamma*(Rpo')*[0;0;1];
        phio(4:6) = -Io*skewvec( wpohat*Ioinv*Iop + Iop*Ioinv*wpohat ) + ...
                    trace(Ioinv*Iop)*Io*skewvec(wpohat) + ...
                    cross( skewvec(wpohat), Io*skewvec(wpohat) );
                
        phih = zeros(6,1);
        phih(1:3) = cross(skewvec(wphhat), mh*x(32:34)) + mo*gamma*(Rph')*[0;0;1];
        phih(4:6) = cross( skewvec(wphhat), Ih*skewvec(wphhat) );
        
        phi = [phio; phih];
        
        % Constraint directions
        Bc = [eye(3), zeros(3,1); zeros(3,3), [0; 0; 1]];
        A = (Bc')*[Adgoc_Inverse, -Adghc_Inverse];
        
        % Augmented system
        Ma = [M, -A'; -A, zeros(4,4)];
        
        
        Vphb = [x(32:34); skewvec(wphhat)];
        Vpob = [x(4:6); skewvec(wpohat)];
        Vocb = zeros(6,1);
        Vhcb = -Adghc_Inverse*Vphb + Adgoc_Inverse*Vpob + Vocb;
        
        Vhcbhat = [skewmat(Vhcb(4:6)), Vhcb(1:3); zeros(1,3), 0];
        Vphbhat = [skewmat(Vphb(4:6)), Vphb(1:3); zeros(1,3), 0];
        Vocbhat = zeros(4,4);
        Vpobhat = [skewmat(Vpob(4:6)), Vpob(1:3); zeros(1,3), 0];
        
        temp1 = bracket(Vhcbhat, ghc_Inverse*Vphbhat*ghc);
        temp2 = bracket(Vocbhat, goc_Inverse*Vpobhat*goc);
        
        Fi = (Bc')*( [temp1(1:3,4); skewvec(temp1(1:3,1:3))] - ...
                     [temp2(1:3,4); skewvec(temp2(1:3,1:3))] );
        
        %% Forces and torques
%         % Grasp matrix
%         G = (Adgoc_Inverse')*Bc;
       

        if t < 0.5
            f = (mo+mh)*[0; 0; 2*gamma - 4*gamma*t];
            tau = [0; 0; 0.2];
        else
            f = zeros(3,1);
            tau = zeros(3,1);
        end
        
        Fh = [f; tau];
        F = [zeros(6,1); Fh];
        
        S = Ma \ [F - phi; Fi];
%         S(13:16)

        if S(3) < 1e-06
            keyboard
        end
        
        %% Equations of motion
        dx = zeros(52, 1);
        dx(1:3) = Rpo*x(4:6);
        dx(4:6) = S(1:3);
        dx(7:15) = reshape(Rpod', 9, 1);
        dx(16:24) = reshape( skewmat(S(4:6))', 9, 1);
        dx(25) = x(26);
        dx(26) = -100*x(26) + 0.1*x(25)*wpohat(2,1)^2;
        dx(27) = x(28);
        dx(28) = -1*x(28) - 2*x(27) - pi*x(25)*x(25)*dot(x(4:6),Rpo(:,3));
        dx(29:31) = Rph*x(32:34);
        dx(32:34) = S(7:9);
        dx(35:43) = reshape(Rphd', 9, 1);
        dx(44:52) = reshape( skewmat(S(10:12))', 9, 1);
    end
end