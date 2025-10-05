clc;
clear;

% Student Number: 110210172 Last digit (n)
n = 2;
N = 6000; 

% Initialize angular velocities for all cases

w = zeros(3, N, 4);

% Derivatives of angular velocities

W = zeros(3, N, 4);

% Initialize Euler angles for all cases

euler = zeros(3, N, 4);

% Derivatives of attitude angles

EULER = zeros(3, N, 4);

% Initial Parameters

% (rad)
euler_initial = [-0.01 - 0.002 * n ;     % -0.01 - 0.02(2)
                  0.01 - 0.002 * n;      % 0.01 - 0.02(2)
                  -0.005 - 0.002 * n];   % -0.005 - 0.002(2)

% (rad/s)
w_initial = [-0.002 - 0.0001 * n;    % -0.002 - 0.0001(2)
             0.003 - 0.0001 * n;     % 0.003 - 0.0001(2) 
             -0.004 - 0.0001 * n];   % -0.004 - 0.0001(2)

for i = 1:4

    euler(:, 1, i) = euler_initial;
    w(:, 1, i) = w_initial;

end

% Mass moments of inertia of the satellite (kg*m^2)

I_1 = 2.1 * 10^-3;
I_2 = 2.0 * 10^-3;
I_3 = 1.9 * 10^-3;

I_min = min([I_1, I_2, I_3]);

% Disturbance torques (N * m)

L_D = 3.6 * 10^-10;       % L_D1 = L_D2 = L_D3

b = ([500; 900; 1200] + 100 * n) * 10^(-9);           % Bias Vector for Magnetometer (T)
v = zeros(3, 6000);

v_rgb = zeros(3, 6000);
v_rg = zeros(3, 6000);

b_rg = zeros(3, 6000);                                % Bias Vector for Gyroscope (rad/s)
b_rg(:, 1) = ([5; 1+n; 1-n] * 10^(-3)) * (pi/180);    % Initial bias Vector for Gyroscope (rad/s)

% Orbital Period

mu = 3.98601 * 10^14;                       % Gravitational parameter (m^3/s^2)
orbit_r = (6371.2 + 800 + 5 * n) * 10^3;    % Orbit radius (m)
T = 2 * pi * sqrt(orbit_r^3 / mu);          % Period (s)

e_m = 80 + 0.5*n + 10;                      % Inclination (deg)
k_B = (4*pi/T) * (1 + sin(e_m)) * I_min;    % Positive controller gain

temp = load('B_igrf.mat');                  % Loading magnetic field values (calculated with IGRF) from hw2
B_igrf = temp.B_ECI * 10^(-9);              % T

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Calculations

C = zeros(3,3,6000);

Bb = zeros(3, N, 4);
Bb2 = zeros(3, N, 4);
Bb3 = zeros(3, N, 4);

W_BN = zeros(3, N, 4);
W_BN2 = zeros(3, N, 4);

L_applied = zeros(3, N, 4);

L_commanded = zeros(3, N, 4);

m_gen = zeros(3, N, 4);

dif_L = zeros(4, N);

W_true = zeros(3, 1, 6000);

b_hat = zeros(3, 6000, 4);

for i = 1:N-1

    v(:, i) = 100e-9 * randn(3,1);
    v_rg(:, i) = 1e-2 * (pi/180) * randn(3,1);
    v_rgb(:, i) = 1e-3 * (pi/180) * randn(3,1);
    b_rg(:, i+1) = b_rg(:, i) + v_rgb(:, i);

    for j = 1:4

        roll = euler(1, i, j); 
        pitch = euler(2, i, j); 
        yaw = euler(3, i, j);

        % Calculating 321 DCM

        C = [cos(pitch)*cos(yaw), cos(pitch)*sin(yaw), -sin(pitch); 
             sin(roll)*sin(pitch)*cos(yaw)-cos(roll)*sin(yaw), ...
             sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw), sin(roll)*cos(pitch); 
             cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw), ...
             cos(roll)*sin(pitch)*sin(yaw)-sin(roll)*cos(yaw), cos(roll)*cos(pitch)];

        Bb(:, i, j) = C * B_igrf(:, i) + b + v(:, i);
        Bb2(:, i, j) = C * B_igrf(:, i) + v(:, i);
        Bb3(:, i, j) = C * B_igrf(:, i);

        W_true = w(:, i, j);
        W_BN(:, i, j) = W_true + b_rg(:, i) + v_rg(:, i);
        W_BN2(:, i, j) = W_true + v_rg(:, i);

        if j == 1

            B_used = Bb(:, i, j); W_used = W_BN(:, i, j);

        elseif j == 2

            B_used = Bb(:, i, j); W_used = W_BN2(:, i, j);

        elseif j == 3

            B_used = Bb2(:, i, j); W_used = W_BN(:, i, j);

        else

            B_used = Bb2(:, i, j); W_used = W_BN2(:, i, j);

        end

        % Control Torque Calculations

        b_hat(:, i, j) = B_used / norm(B_used);
        L_commanded(:, i, j) = -k_B * (eye(3) - b_hat(:, i, j) * b_hat(:, i, j)') * W_used;
        m_gen(:, i, j) = -(k_B / norm(B_used)) * cross(b_hat(:, i, j), W_used);
        L_applied(:, i, j) = cross(m_gen(:, i, j), Bb3(:, i, j));

        dif_L(j, i) = norm(L_commanded(:, i, j) - L_applied(:, i, j));

        % Euler integration
        W(1, i, j) = (-(I_3 - I_2) * w(2, i, j) * w(3, i, j) + L_D + L_applied(1, i, j))/I_1;
        W(2, i, j) = (-(I_1 - I_3) * w(3, i, j) * w(1, i, j) + L_D + L_applied(2, i, j))/I_2;
        W(3, i, j) = (-(I_2 - I_1) * w(1, i, j) * w(2, i, j) + L_D + L_applied(3, i, j))/I_3;

        w(:, i+1, j) = w(:, i, j) + W(:, i, j);

        % Euler angles

        EULER(1, i, j) = (cos(pitch) * w(1, i, j) + sin(roll) * sin(pitch) * w(2, i, j) + cos(roll)*sin(pitch) * w(3, i, j)) / cos(pitch);
        EULER(2, i, j) = (cos(roll) * cos(pitch) * w(2, i, j) - sin(roll) * cos(pitch) * w(3, i, j)) / cos(pitch);
        EULER(3, i, j) = (sin(roll) * w(2, i, j) + cos(roll) * w(3, i, j)) / cos(pitch);

        euler(:, i+1, j) = wrapToPi(euler(:, i, j) + EULER(:, i, j));

    end
end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Plotting

x = 1:6000;

% L Applied

% Case 1

figure;
plot(x, L_applied(1, :, 1), 'r');
hold on;
plot(x, L_applied(2, :, 1), 'g');
plot(x, L_applied(3, :, 1), 'b');
xlabel('Time (s)');
ylabel('Applied Torque (N·m)');
title('Applied Control Torque - Case 1 (Bb1 + Wb1)');
legend('Lx', 'Ly', 'Lz');
grid on;

% Case 2

figure;
plot(x, L_applied(1, :, 2), 'r');
hold on;
plot(x, L_applied(2, :, 2), 'g');
plot(x, L_applied(3, :, 2), 'b');
xlabel('Time (s)');
ylabel('Applied Torque (N·m)');
title('Applied Control Torque - Case 2 (Bb1 + Wb2)');
legend('Lx', 'Ly', 'Lz');
grid on;

% Case 3

figure;
plot(x, L_applied(1, :, 3), 'r');
hold on;
plot(x, L_applied(2, :, 3), 'g');
plot(x, L_applied(3, :, 3), 'b');
xlabel('Time (s)');
ylabel('Applied Torque (N·m)');
title('Applied Control Torque - Case 3 (Bb2 + Wb1)');
legend('Lx', 'Ly', 'Lz');
grid on;

% Case 4

figure;
plot(x, L_applied(1, :, 4), 'r');
hold on;
plot(x, L_applied(2, :, 4), 'g');
plot(x, L_applied(3, :, 4), 'b');
xlabel('Time (s)');
ylabel('Applied Torque (N·m)');
title('Applied Control Torque - Case 4 (Bb2 + Wb2)');
legend('Lx', 'Ly', 'Lz');
grid on;

% L Difference

figure;
plot(x, dif_L(1, :));
xlabel('Time (s)');
ylabel('|Lcommanded - Lapplied|');
title('Difference for Case 1');
grid on;

figure;
plot(x, dif_L(2, :));
xlabel('Time (s)');
ylabel('|Lcommanded - Lapplied|');
title('Difference for Case 2');
grid on;

figure;
plot(x, dif_L(3, :));
xlabel('Time (s)');
ylabel('|Lcommanded - Lapplied|');
title('Difference for Case 3');
grid on;

figure;
plot(x, dif_L(4, :));
xlabel('Time (s)');
ylabel('|Lcommanded - Lapplied|');
title('Difference for Case 4');
grid on;

% Angular Velocities

% Case 1

figure;
plot(x, w(1, :, 1) * (180/pi));
hold on;
plot(x, w(2, :, 1) * (180/pi));
plot(x, w(3, :, 1) * (180/pi));
xlabel('Time (s)');
ylabel('Angular Velocity (deg/s)');
title('Angular Velocities Case 1');
legend('w1', 'w2', 'w3');
grid on;

% Case 2

figure;
plot(x, w(1, :, 2) * (180/pi));
hold on;
plot(x, w(2, :, 2) * (180/pi));
plot(x, w(3, :, 2) * (180/pi));
xlabel('Time (s)');
ylabel('Angular Velocity (deg/s)');
title('Angular Velocities Case 2');
legend('w1', 'w2', 'w3');
grid on;

% Case 3

figure;
plot(x, w(1, :, 3) * (180/pi));
hold on;
plot(x, w(2, :, 3) * (180/pi));
plot(x, w(3, :, 3) * (180/pi));
xlabel('Time (s)');
ylabel('Angular Velocity (deg/s)');
title('Angular Velocities Case 3');
legend('w1', 'w2', 'w3');
grid on;

% Case 4

figure;
plot(x, w(1, :, 4) * (180/pi));
hold on;
plot(x, w(2, :, 4) * (180/pi));
plot(x, w(3, :, 4) * (180/pi));
xlabel('Time (s)');
ylabel('Angular Velocity (deg/s)');
title('Angular Velocities Case 4');
legend('w1', 'w2', 'w3');
grid on;

% Euler Angles

% Case 1

figure;
plot(x, euler(1, :, 1) * (180/pi));
hold on;
plot(x, euler(2, :, 1) * (180/pi));
plot(x, euler(3, :, 1) * (180/pi));
xlabel('Time (s)');
ylabel('Angle (deg)');
title('Euler Angles Case 1');
ylim([-180 180]);
legend('roll', 'pitch', 'yaw');
grid on;

% Case 2

figure;
plot(x, euler(1, :, 2) * (180/pi));
hold on;
plot(x, euler(2, :, 2) * (180/pi));
plot(x, euler(3, :, 2) * (180/pi));
xlabel('Time (s)');
ylabel('Angle (deg)');
title('Euler Angles Case 2');
ylim([-180 180]);
legend('roll', 'pitch', 'yaw');
grid on;

% Case 3

figure;
plot(x, euler(1, :, 3) * (180/pi));
hold on;
plot(x, euler(2, :, 3) * (180/pi));
plot(x, euler(3, :, 3) * (180/pi));
xlabel('Time (s)');
ylabel('Angle (deg)');
title('Euler Angles Case 3');
ylim([-180 180]);
legend('roll', 'pitch', 'yaw');
grid on;

% Case 4

figure;
plot(x, euler(1, :, 4) * (180/pi));
hold on;
plot(x, euler(2, :, 4) * (180/pi));
plot(x, euler(3, :, 4) * (180/pi));
xlabel('Time (s)');
ylabel('Angle (deg)');
title('Euler Angles Case 4');
ylim([-180 180]);
legend('roll', 'pitch', 'yaw');
grid on;






