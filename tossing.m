function [t,x] = tossing()

t0 = 0; tf = 2;
x0 = zeros(6,1); % States for translation
x0 = [x0; reshape(eye(3), 9, 1); zeros(9,1)]; % States for rotation
x0 = [x0; 0.1; 0]; % States for the radial deformation
x0 = [x0; 0.0; 0]; % States for the vertical deformation
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
par.m = 0.5;           % [kg] 
% par.r0 = 0.3048;   % [m]
par.nu = 0.35;       
par.E = 5e+06;       % [Pa] = [kg/m/s/s]
par.rho = 1.2e+03;   % [kg/m/m/m]
par.g = 9.81;        % [m/s/s]
par.I = diag([1/4*par.m*x0(25)^2, 1/4*par.m*x0(25)^2, 1/2*par.m*x0(25)^2]);    % [kg m^2]

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
        m = par.m;
        g = par.g;
        nu = par.nu;
        E = par.E;
        rho = par.rho;
        
        R = reshape(x(7:15), 3, 3)';
        what = reshape(x(16:24), 3, 3)';
        Rp = R*what;
        
        I = m/18*diag( [3*x(25)^2 + x(27)^2, 3*x(25)^2 + x(27)^2, 1/6*x(25)^2] );
        Ip = m/18*diag( [6*x(25)*x(26) + 2*x(27)*x(28), 6*x(25)*x(26) + ...
                        2*x(27)*x(28), 1/3*x(25)*x(26)] );
%         Ip = m*x(25)/2*x(26)*[1,0,0;0,1,0;0,0,2];
%         I = par.I;
        Iinv = inv(I);
        
%         f = m*[0;0;g - 6*(x(3)-1) - 5*x(6)];
%         tau = skewmat([0; 0; -5*what(2,1)-6*(atan2(R(2,1), R(1,1))-pi/3)]);

        if t < 0.5
            f = m*[0; 0; 2*g];
            tau = skewmat([0; 0.02; 0.02]);
        else
            f = 0;
            tau = skewmat(zeros(3,1));
        end
        
        dx = zeros(26, 1);
        dx(1:3) = x(4:6);
        dx(4:6) = 1/m*(f - m*[0;0;g]);
        dx(7:15) = reshape(Rp', 9, 1);
        dx(16:24) = reshape( ( Ip*Iinv*what - (Ip*Iinv*what)' - ...
                               trace(Iinv*Ip)*what + ...
                               (I*what*Iinv*what)' - I*what*Iinv*what + ...
                               1/det(I)*I*tau*I)', 9, 1);
        dx(25) = x(26);
        dx(26) = -100*x(26) + 0.1*x(25)*what(2,1)^2;
        dx(27) = x(28);
        dx(28) = -1*x(28) - 2*x(27) - pi*x(25)*x(25)*dot(x(4:6),R(:,3));
    end
end