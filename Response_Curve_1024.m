%% Numerical simulation for dynamic response curve
clear all;clc;
% Color Order: gray, purple, blue, green, orange, red
color_full = ["#5f5f5f",  "#7262ac","#2e7ebb" ,"#2e974e" ,"#e25508","#d92523"   ];
colors_2 = ["#cecece",  "#cfcfe5","#b7d4ea","#b8e3b2","#fdc38d" ,"#fcab8f"]; 

%% 10Hz, mathematical
clear all;clc;
% Order: gray, purple, blue, green, orange, red
color_full = ["#5f5f5f",  "#7262ac","#2e7ebb" ,"#2e974e" ,"#e25508","#d92523"   ];
colors_2 = ["#cecece",  "#cfcfe5","#b7d4ea","#b8e3b2","#fdc38d" ,"#fcab8f"]; 

a = 50; %mm -> m
b = 80; %mm
g = 9.8; 
x0 = 5;
phi_0 = asin((b-x0)/b);

k = (300/1E3*g)/(80); % N/mm

% x1 = [0:1:45]; % Displacement
phi_1 = [phi_0:-0.005:asin((b-50)/b)];

EF_0 = sqrt(3)*b*cos(phi_0)*ones(size(phi_1));
EF_1 = sqrt(3)*b*cos(phi_1);

% disp = 2*b*ones(size(phi_1))-2*b*sin(phi_1);
F_Ver = 9/2*k*b*sin(phi_1).*(cos(phi_1) - cos(phi_0)*ones(size(phi_1)) )./cos(phi_1);



%---------- Dynamic -----

M  = max(F_Ver)/g;
kh = k;

fa_square = 5*0.2^2; % constant, check article for proving stage

% Input single f_values value
f_values = 10; % Hz

% Calculate za for the single f_values value
za = sqrt(fa_square / f_values);

tspan_1 = [0 4];
phi0 = [phi_0; 0]; % IC: phi(0) = phi_0, phi'(0) = 0

% Pass parameters to myODE_dyn_1 using an anonymous function
[t, phi_dyn] = ode45(@(t, phi) myODE_dyn_1(t, phi, b, phi_0, g, f_values, za, kh, M), tspan_1, phi0);

% Reinitialize phi_all to match the size of phi_dyn
phi_all = zeros(length(phi_dyn), 3);

% 2nd order derivative of phi
phi_ddot = 1 ./ (2 * b * cos(phi_dyn(:, 1))) .* (2 * b * sin(phi_dyn(:, 1)) .* phi_dyn(:, 2).^2 + ...
    4 * pi^2 * f_values^2 * za * sin(2 * pi * f_values * t) - g + (9 * kh * b) / (2 * M) * tan(phi_dyn(:, 1)) .* (cos(phi_dyn(:, 1)) - cos(phi_0) * ones(size(phi_dyn(:, 1)))));

% Combine into one parameter
phi_all(:, 1) = phi_dyn(:, 1);  % phi
phi_all(:, 2) = phi_dyn(:, 2);  % phi_dot
phi_all(:, 3) = phi_ddot(:, 1); % phi_ddot

% Calculate z_ddot
z_ddot = -g + 1 ./ (2 * M * b * cos(phi_all(:, 1))) .* (3 / 2 * kh * b^2 * sin(phi_all(:, 1)) .* cos(phi_all(:, 1)) + ...
    4 * M * b^2 * sin(phi_all(:, 1)) .* cos(phi_all(:, 1)) .* phi_all(:, 2).^2 - 4 * M * b^2 * cos(phi_all(:, 1)).^2 .* phi_all(:, 3));

y_ddot = z_ddot + 2 * b * cos(phi_all(:, 1)) .* phi_all(:, 3) - 2 * b * sin(phi_all(:, 1)) .* phi_all(:, 2).^2;
% y_ddot = 2*b*cos(phi_all(:, 1)) .* phi_all(:, 3) - 2 * b * sin(phi_all(:, 1)) .* phi_all(:, 2).^2 - 4*pi^2*f_values^2*za*sin(2*pi*f_values*t);



T_math = t(1:end-150+1);
Ydd_math = y_ddot(90:end-60)/1E4;
Zdd_math = z_ddot(90:end-60)/1E4;

Ydd_math = Ydd_math - mean(Ydd_math);
% Plot
figure()
hold on
% plot(T_math, 5E4 * (Ydd_math), 'Color', color_full(6), 'LineWidth', 1.2)
plot(T_math(1:end-29), 1*Zdd_math(30:end), 'Color', color_full(6), 'LineWidth', 1.2)
plot(T_math, 10*Zdd_math, 'Color', color_full(4), 'LineWidth', 1.2)

hold off
title('Numerical Result, 10Hz', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('Output', 'Input', 'FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('Time [s]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Acceleration [m/s^2]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
axis([0 2 -2 2])
grid on
box on