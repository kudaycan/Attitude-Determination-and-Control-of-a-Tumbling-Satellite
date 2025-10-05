
clc;
clear;

% Student Number: 110210172 Last digit (n)

n = 2;
N = 6000;

% Initial and constant values 

b = ([500; 900; 1200] + 100 * n);

temp = load('DCM.mat');                % Loading DCM from hw1
DCM = temp.B; 

temp2 = load('B_igrf.mat');            % Loading magnetic field values (calculated with IGRF) from hw2
B_igrf = temp2.B_ECI;

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Magnetometer Measurements

v = zeros(3, 6000);

Bb = zeros(3, 6000);
Bb2 = zeros(3, 6000);
Bb3 = zeros(3, 6000);

Bb_norm = zeros(3, 6000);
Bb_norm2 = zeros(3, 6000);
Bb_norm3 = zeros(3, 6000);

for i = 1:N

    v(:,i) = 100 * randn(3,1);

    Bb(:, i) = DCM(:, :, i) * B_igrf(:, i) + b + v(:, i);
    Bb2(:, i) = DCM(:, :, i) * B_igrf(:, i) + v(:, i);
    Bb3(:, i) = DCM(:, :, i) * B_igrf(:, i);

    Bb_norm(:, i) = Bb(:, i) / norm(Bb(:, i));
    Bb_norm2(:, i) = Bb2(:, i) / norm(Bb2(:, i));
    Bb_norm3(:, i) = Bb3(:, i) / norm(Bb3(:, i));

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Plotting

dt = 1:6000;                           % Time  

% Plotting Magnetic Field

figure;
plot(dt, Bb(1, :));
hold on;
plot(dt, Bb2(1, :));
plot(dt, Bb3(1, :));
hold off;
title('Magnetic Field in x Direction');
xlabel('Time (s)');
ylabel('B (nT)');
legend('Eq. 1', 'Eq. 2', 'Eq. 3');
grid on;

figure;
plot(dt, Bb(2, :));
hold on;
plot(dt, Bb2(2, :));
plot(dt, Bb3(2, :));
hold off;
title('Magnetic Field in y Direction');
xlabel('Time (s)');
ylabel('B (nT)');
legend('Eq. 1', 'Eq. 2', 'Eq. 3');
grid on;

figure;
plot(dt, Bb(3, :));
hold on;
plot(dt, Bb2(3, :));
plot(dt, Bb3(3, :));
hold off;
title('Magnetic Field in z Direction');
xlabel('Time (s)');
ylabel('B (nT)');
legend('Eq. 1', 'Eq. 2', 'Eq. 3');
grid on;

% Plotting (Magnetic Field) / (Vector Norm)

figure;
plot(dt, Bb_norm(1, :));
hold on;
plot(dt, Bb_norm2(1, :));
plot(dt, Bb_norm3(1, :));
hold off;
title('(B_b / Norm) in x Direction');
xlabel('Time (s)');
ylabel('B');
legend('Eq. 1', 'Eq. 2', 'Eq. 3');
grid on;

figure;
plot(dt, Bb_norm(2, :));
hold on;
plot(dt, Bb_norm2(2, :));
plot(dt, Bb_norm3(2, :));
hold off;
title('(B_b / Norm) in y Direction');
xlabel('Time (s)');
ylabel('B');
legend('Eq. 1', 'Eq. 2', 'Eq. 3');
grid on;

figure;
plot(dt, Bb_norm(3, :));
hold on;
plot(dt, Bb_norm2(3, :));
plot(dt, Bb_norm3(3, :));
hold off;
title('(B_b / Norm) in z Direction');
xlabel('Time (s)');
ylabel('B');
legend('Eq. 1', 'Eq. 2', 'Eq. 3');
grid on;



