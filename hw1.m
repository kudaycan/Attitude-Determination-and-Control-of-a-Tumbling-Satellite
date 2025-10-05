
clc;
clear;

% Student Number: 110210172 Last digit (n)
n = 2;
N = 5999;

% Angular velocities 

w_1 = zeros(1, 6000); 
w_2 = zeros(1, 6000); 
w_3 = zeros(1, 6000); 
w_1r = w_1;
w_2r = w_2;
w_3r = w_3;

% Initial angular velocities (rad/s)

w_1(1) = -0.002 - 0.0001 * n; % -0.002 - 0.0001(2)
w_2(1) = 0.003 - 0.0001 * n; % 0.003 - 0.0001(2)
w_3(1) = -0.004 - 0.0001 * n; % -0.004 - 0.0001(2)
w_1r(1) = w_1(1);
w_2r(1) = w_2(1);
w_3r(1) = w_3(1);


% Derivatives of angular velocities

W1 = zeros(1, 6000); 
W2 = zeros(1, 6000); 
W3 = zeros(1, 6000); 

% Attitude angles

r = zeros(1, 6000); 
p = zeros(1, 6000); 
y = zeros(1, 6000); 
rr = r; 
pr = p; 
yr = y; 

% Differences between angles

d_1 = zeros(1, 6000); % rr - r
d_2 = zeros(1, 6000); % pr - p
d_3 = zeros(1, 6000); % yr - y


% Initial attitude angles (rad)

r(1) = -0.01 - 0.02 * n; % -0.01 - 0.02(2)
p(1) = 0.01 - 0.02 * n; % 0.01 - 0.02(2)
y(1) = -0.005 - 0.002 * n;  % -0.005 - 0.002(2)
rr(1) = r(1); 
pr(1) = p(1); 
yr(1) = y(1);  

% Derivatives of attitude angles

R = zeros(1, 6000);
P = zeros(1, 6000); 
Y = zeros(1, 6000);

% Mass moments of inertia of the satellite (kg*m^2)

I_1 = 2.1 * 10^-3;
I_2 = 2.0 * 10^-3;
I_3 = 1.9 * 10^-3;

% Disturbance torques (N*m)

L_1 = 3.6 * 10^-10;
L_2 = 3.6 * 10^-10;
L_3 = 3.6 * 10^-10;

% --------------------------------------------------------------------------------------------------------------------------------------------------

% For loop for the calculation of angular velocities using Euler's Method

for i = 1:N

    W1(i) = (-(I_3 - I_2) * w_2(i) * w_3(i) + L_1) / I_1;
    W2(i) = (-(I_1 - I_3) * w_3(i) * w_1(i) + L_2) / I_2;
    W3(i) = (-(I_2 - I_1) * w_1(i) * w_2(i) + L_3) / I_3;

    w_1(i + 1) = w_1(i) + W1(i) * 1;
    w_2(i + 1) = w_2(i) + W2(i) * 1;
    w_3(i + 1) = w_3(i) + W3(i) * 1;

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% For loop for the calculation of attitude angles using Euler's Method

for i = 1:N

    R(i) = (cos(p(i)) * w_1(i) + sin(r(i)) * sin(p(i)) * w_2(i) + cos(r(i)) * sin(p(i)) * w_3(i)) / cos(p(i));
    P(i) = (cos(r(i)) * cos(p(i)) * w_2(i) - sin(r(i)) * cos(p(i)) * w_3(i)) / cos(p(i));
    Y(i) = (sin(r(i)) * w_2(i) + cos(r(i)) * w_3(i)) / cos(p(i));

    r(i + 1) = r(i) + R(i) * 1;
    p(i + 1) = p(i) + P(i) * 1;
    y(i + 1) = y(i) + Y(i) * 1;
   
end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% For loop for the calculation of attitude angles using 4th Order Runge Kutta Method

for i = 1:N
    
    f_roll = @(r, p, w1, w2, w3) (cos(p) * w1 + sin(r) * sin(p) * w2 + cos(r) * sin(p) * w3) / cos(p);
    f_pitch = @(r, p, w2, w3) (cos(r) * cos(p) * w2 - sin(r) * cos(p) * w3) / cos(p);
    f_yaw = @(r, p, w2, w3) (sin(r) * w2 + cos(r) * w3) / cos(p);  

    k1_r = f_roll(rr(i), pr(i), w_1(i), w_2(i), w_3(i));
    k1_p = f_pitch(rr(i), pr(i), w_2(i), w_3(i));
    k1_y = f_yaw(rr(i), pr(i), w_2(i), w_3(i));  

    k2_r = f_roll(rr(i) + 0.5 * k1_r, pr(i) + 0.5 * k1_p, w_1(i), w_2(i), w_3(i));
    k2_p = f_pitch(rr(i) + 0.5 * k1_r, pr(i) + 0.5 * k1_p, w_2(i), w_3(i));
    k2_y = f_yaw(rr(i) + 0.5 * k1_r, pr(i) + 0.5 * k1_p, w_2(i), w_3(i));  

    k3_r = f_roll(rr(i) + 0.5 * k2_r, pr(i) + 0.5 * k2_p, w_1(i), w_2(i), w_3(i));
    k3_p = f_pitch(rr(i) + 0.5 * k2_r, pr(i) + 0.5 * k2_p, w_2(i), w_3(i));
    k3_y = f_yaw(rr(i) + 0.5 * k2_r, pr(i) + 0.5 * k2_p, w_2(i), w_3(i));  

    k4_r = f_roll(rr(i) + k3_r, pr(i) + k3_p, w_1(i), w_2(i), w_3(i));
    k4_p = f_pitch(rr(i) + k3_r, pr(i) + k3_p, w_2(i), w_3(i));
    k4_y = f_yaw(rr(i) + k3_r, pr(i) + k3_p, w_2(i), w_3(i));  

    rr(i+1) = rr(i) + (1/6) * (k1_r + 2*k2_r + 2*k3_r + k4_r);
    pr(i+1) = pr(i) + (1/6) * (k1_p + 2*k2_p + 2*k3_p + k4_p);
    yr(i+1) = yr(i) + (1/6) * (k1_y + 2*k2_y + 2*k3_y + k4_y);

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Calculating the differences between the methods for each Euler angles

for i = 1:N

    d_1(i) = rr(i) - r(i);  % roll
    d_2(i) = pr(i) - p(i);  % pitch
    d_3(i) = yr(i) - y(i);  % yaw

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Applying mod to get yaw angle between +180 and -180 degrees

for i = 1:N

    yr(i+1) = mod(yr(i+1) + pi, 2*pi) - pi;
    y(i+1) = mod(y(i+1) + pi, 2*pi) - pi;
    
end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Calculating vector norms (Euler)

n2 = norm([y(2000);p(2000);r(2000)]) * (180/pi);

n4 = norm([y(4000);p(4000);r(4000)]) * (180/pi);

n6 = norm([y(6000);p(6000);r(6000)]) * (180/pi);

% Calculating vector norms (Runge-Kutta)

n2r = norm([yr(2000);pr(2000);rr(2000)]) * (180/pi);

n4r = norm([yr(4000);pr(4000);rr(4000)]) * (180/pi);

n6r = norm([yr(6000);pr(6000);rr(6000)]) * (180/pi);

% Constructing DCM

B = zeros(3,3,6000);

for i = 1:6000

    % 321 sequence

    B(:,:,i) = [cos(p(i))*cos(y(i)), cos(p(i))*sin(y(i)), -sin(p(i)); 
         sin(r(i))*sin(p(i))*cos(y(i)) - cos(r(i))*sin(y(i)), sin(r(i))*sin(p(i))*sin(y(i)) + cos(r(i))*cos(y(i)), sin(r(i))*cos(p(i)); 
         cos(r(i))*sin(p(i))*cos(y(i)) + sin(r(i))*sin(y(i)), cos(r(i))*sin(p(i))*sin(y(i)) - sin(r(i))*cos(y(i)), cos(r(i))*cos(p(i))];

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Calculating 3-1-3 angles

r_313 = zeros(1, 6000);      % roll
y_313 = zeros(1, 6000);      % first yaw angle
y_2_313 = zeros(1, 6000);    % final yaw 

for i = 1:6000
    
    r_313(i) = acos(B(3,3,i));                        % Roll angle
    y_313(i) = atan(B(1,3,i) / -B(2,3,i));            % First yaw
    y_2_313(i) = atan(B(3,1,i) / B(3,2,i));           % Final yaw

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

x = 1:6000;

% Plotting angular velocities

figure;
plot(x, w_1 * (180/pi), 'b-', 'LineWidth', 1);
hold on;
plot(x, w_2 * (180/pi), 'r-', 'LineWidth', 1);
plot(x, w_3 * (180/pi), 'g-', 'LineWidth', 1);
hold off;
title('Euler Method');
xlabel('Time (s)');
ylabel('Angular Velocity (deg/s)');
legend('w_1', 'w_2', 'w_3');
grid on;

% Plotting attitude angles

% Roll

figure;

plot(x, r * (180/pi), 'b-', 'LineWidth', 1);
hold on;
plot(x, rr * (180/pi), 'g-', 'LineWidth', 1);
hold off;
title('Roll Angle');
xlabel('Time (s)');
ylabel('Angle (deg)');
ylim([-180 180]);
legend('Euler', 'Runge-Kutta');
grid on;

% Pitch

figure;

plot(x, p * (180/pi), 'g-', 'LineWidth', 1);
hold on;
plot(x, pr * (180/pi), 'r-', 'LineWidth', 1);
hold off;
title('Pitch Angle');
xlabel('Time (s)');
ylabel('Angle (deg)');
ylim([-180 180]);
legend('Euler', 'Runge-Kutta');
grid on;

% Yaw

figure;

plot(x, y * (180/pi), 'b-', 'LineWidth', 1);
hold on;
plot(x, yr * (180/pi), 'r-', 'LineWidth', 1);
hold off;
title('Yaw Angle');
xlabel('Time (s)');
ylabel('Angle (deg)');
ylim([-180 180]);
legend('Euler', 'Runge-Kutta');
grid on;

% Plotting differences

figure;

plot(x, d_1 * (180/pi), 'b', 'LineWidth', 1);
hold on;
plot(x, d_2 * (180/pi), 'r', 'LineWidth', 1);
plot(x, d_3 * (180/pi), 'g', 'LineWidth', 1);
hold off;
title('RK4 - Euler Angles');
xlabel('Time (s)');
ylabel('Angle (deg)');
legend('roll', 'pitch', 'yaw');
grid on;

% Displaying vector norms

disp('Vector norm at 2000s in Euler Method:');
disp(n2);

disp('Vector norm at 2000s in Runge-Kutta Method:');
disp(n2r);

disp('Vector norm at 4000s in Euler Method:');
disp(n4);

disp('Vector norm at 4000s in Runge-Kutta Method:');
disp(n4r);

disp('Vector norm at 6000s in Euler Method:');
disp(n6);

disp('Vector norm at 6000s in Runge-Kutta Method:');
disp(n6r);

% Displaying DCM

disp('DCM at 2000s in Euler Method:');
disp(B(:,:,2000));

disp('DCM at 4000s in Euler Method:');
disp(B(:,:,4000));

disp('DCM at 6000s in Euler Method:');
disp(B(:,:,6000));

% Plot 3-1-3 Euler Angles

figure;

plot(x, y_313 * (180/pi), 'r', 'LineWidth', 1);
hold on;
plot(x, r_313 * (180/pi), 'g', 'LineWidth', 1);
plot(x, y_2_313 * (180/pi), 'b', 'LineWidth', 1);
hold off;
title('3-1-3 Euler Angles Over Time');
xlabel('Time (s)');
ylabel('Angle (deg)');
ylim([-180 180]);
legend('Yaw', 'Roll', 'Yaw(2)');
grid on;




