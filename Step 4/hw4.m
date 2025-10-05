
clc;
clear;

% Student Number: 110210172 Last digit (n)

n = 2;
N = 6000;

% Initial and constant values 

b_rg = zeros(3, 6000);
b_rg(:, 1) = ([5; 1+n; 1-n] * 10^(-3));

temp = load('w_1.mat');                % Loading w in x direction from hw1
w_x = rad2deg(temp.w_1); 

temp2 = load('w_2.mat');               % Loading w in y direction from hw1
w_y = rad2deg(temp2.w_2); 

temp3 = load('w_3.mat');               % Loading w in z direction from hw1
w_z = rad2deg(temp3.w_3);

W_true = [w_x; w_y; w_z];

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Calculating Rate Gyro Bias Vector with Euler

v_rgb = zeros(3, 6000);

for i = 1:5999

    v_rgb(:, i) = 1 * 10^(-3) * randn(3,1);

    b_rg(:, i + 1) = b_rg(:, i) + v_rgb(:, i);

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Calculating Gyroscope Measurements

v_rg = zeros(3, 6000);

W_BN = zeros(3, 6000);
W_BN2 = zeros(3, 6000);
W_BN3 = zeros(3, 6000);

for i = 1:N

    v_rg(:, i) = 1 * 10^(-2) * randn(3,1);

    W_BN(:, i) = W_true(:, i) + b_rg(:, i) + v_rg(:, i);
    W_BN2(:, i) = W_true(:, i) + v_rg(:, i);
    W_BN3(:, i) = W_true(:, i);

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Plotting

dt = 1:6000;                           % Time  

% Plotting Rate Gyro Bias Vector

figure;
plot(dt, b_rg(1, :));
hold on;
plot(dt, b_rg(2, :));
plot(dt, b_rg(3, :));
title('Bias Vector');
xlabel('Time (s)');
ylabel('b (deg/s)');
legend('x', 'y', 'z');
grid on;

% Plotting Angular Velocity

figure;
plot(dt, W_BN(1, :));
hold on;
plot(dt, W_BN2(1, :));
plot(dt, W_BN3(1, :), 'k', 'LineWidth', 3);
hold off;
title('Angular Velocity x Direction');
xlabel('Time (s)');
ylabel('w (deg/s)');
legend('Case 1', 'Case 2', 'Case 3');
grid on;

figure;
plot(dt, W_BN(2, :));
hold on;
plot(dt, W_BN2(2, :));
plot(dt, W_BN3(2, :), 'k', 'LineWidth', 3);
hold off;
title('Angular Velocity y Direction');
xlabel('Time (s)');
ylabel('w (deg/s)');
legend('Case 1', 'Case 2', 'Case 3');
grid on;

figure;
plot(dt, W_BN(3, :));
hold on;
plot(dt, W_BN2(3, :));
plot(dt, W_BN3(3, :), 'k', 'LineWidth', 3);
hold off;
title('Angular Velocity z Direction');
xlabel('Time (s)');
ylabel('w (deg/s)');
legend('Case 1', 'Case 2', 'Case 3');
grid on;




