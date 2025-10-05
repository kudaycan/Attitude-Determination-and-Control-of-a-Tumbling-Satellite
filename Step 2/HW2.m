
clc;
clear;

% Student Number: 110210172 Last digit (n)

n = 2;
N = 6000;

% Initial and constant values 

r_e = 6.3712 * 10^6;                  % Earth's radius (m)
r = (6371.2 + 800 + 5 * n) * 10^3;    % Orbit radius (m) 
mu = 3.98601 * 10^14;                 % Gravitational parameter (m^3/s^2)
inc = 80 + 0.5 * n;                   % Inclination (deg)
omg = 2 * n;                          % RAAN (deg)
w_Earth = 7.29212 * 10^(-5);          % Earth's rotational velocity (rad/s)
deltat = 1;                           % Sampling time (s)

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Period and sample calculation

T = 2 * pi * sqrt(r^3 / mu);          % Period (s)
N_o = ceil(T / deltat);               % Number of samples for 1 orbit

disp('Number of samples for a circular Earth orbiting satellite:');
fprintf('%d\n', N_o);

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Calculating the Position of the Satellite

x = zeros(1, 6000);
y = zeros(1, 6000);
z = zeros(1, 6000);
v_t = zeros(1, 6000);                 % True anomaly
pos = zeros(3, 6000);                 % Position vector

for i = 1:N

    v_t(i) = (360 / T) * i;

    x(i) = r * (cosd(0 + v_t(i)) * cosd(omg) - sind(0 + v_t(i)) * sind(omg) * cosd(inc));
    y(i) = r * (cosd(0 + v_t(i)) * sind(omg) + sind(0 + v_t(i)) * cosd(omg) * cosd(inc));
    z(i) = r * (sind(0 + v_t(i)) * sind(inc));

    % Assigning values to position vector

    pos(:, i) = [x(i), y(i), z(i)];

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Transforming Position Vector (ECI to ECEF) and Calculating Specific Parameters

% Date, time

utc1 = [2024, 4, 13, 20, 59, 0];      % Due date of this homework

jd1 = juliandate(utc1);

jd = zeros(1, 6000);
utc = NaT(1, 6000);

pos_ecef = zeros(3, 6000);            % Position vector
alt = zeros(1, 6000);
long = zeros(1, 6000);
lat = zeros(1, 6000);
theta = zeros(1, 6000);
phi = zeros(1, 6000);

for i = 1:N

    jd(i) = jd1 + i/24/60/60 ;

    utc(i) = datetime(jd(i), 'convertfrom', 'juliandate'); % Since we calculate utc(i) values in this for loop we can use them in other for loops without any other calculations

    pos_ecef(:, i) = eci2ecef(utc(i), pos(:, i));

    a = pos_ecef(1, i);
    b = pos_ecef(2, i);
    c = pos_ecef(3, i);

    theta(i) = atan2(b, a);           % Azimuthal angle in radians     
    phi(i) = acos(c / r);             % Polar angle in radians            

    % Altitude calculation

    alt(i) = sqrt(x(i)^2 + y(i)^2 + z(i)^2) - r_e;

    % Longitude and latitude calculation

    long(i) = (theta(i)) * (180/pi);   
    lat(i) = 90 - (phi(i) * (180/pi)); 

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Magnetic Field Calculations

h = zeros(1, N);
d = zeros(1, N);
I = zeros(1, N);
f = zeros(1, N);

% ...Using "igrfmagm" 

B_NED = zeros(3, N);   

Bx_ECEF = zeros(1, N);  
By_ECEF = zeros(1, N);
Bz_ECEF = zeros(1, N); 

B_ECI = zeros(3, N);

% ...Using "wrldmagm"

B_NED2 = zeros(3, N);   

Bx_ECEF2 = zeros(1, N);  
By_ECEF2 = zeros(1, N);
Bz_ECEF2 = zeros(1, N); 

B_ECI2 = zeros(3, N);

% Differences

dif_x = zeros(1, N);
dif_y = zeros(1, N);
dif_z = zeros(1, N);


for i = 1:N

    % igrfmagm

    [B_NED(:, i), h(i), d(i), I(i), f(i),] = igrfmagm(alt(i), lat(i), long(i), decyear(2024, 4, 13), 13);

    % NED to ECEF

    [Bx_ECEF(i), By_ECEF(i), Bz_ECEF(i)] = ned2ecefv(B_NED(1, i), B_NED(2, i), B_NED(3, i), lat(1), long(1), "degrees");

    % ECEF to ECI

    B_ECI(:, i) = ecef2eci(utc(i), [Bx_ECEF(i), By_ECEF(i), Bz_ECEF(i)]);

    % wrldmagm

    [B_NED2(:, i), h(i), d(i), I(i), f(i),] = wrldmagm(alt(i), lat(i), long(i), decyear(2024, 4, 13), '2020');

    % NED to ECEF

    [Bx_ECEF2(i), By_ECEF2(i), Bz_ECEF2(i)] = ned2ecefv(B_NED2(1, i), B_NED2(2, i), B_NED2(3, i), lat(1), long(1), "degrees");

    % ECEF to ECI

    B_ECI2(:, i) = ecef2eci(utc(i), [Bx_ECEF2(i), By_ECEF2(i), Bz_ECEF2(i)]);

    % Differences in two methods

    dif_x(i) = B_ECI(1, i) - B_ECI2(1, i);
    dif_y(i) = B_ECI(2, i) - B_ECI2(2, i);
    dif_z(i) = B_ECI(3, i) - B_ECI2(3, i);

end

% --------------------------------------------------------------------------------------------------------------------------------------------------

% Plotting

dt = 1:6000;                           % Time  

% 3D Plot

figure;
plot3(x, y, z);
grid on;
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('Position of the Satellite (3D)');

% XY Plot

figure;
plot(x, y);
xlabel('x (m)');
ylabel('y (m)');
title('Position of the Satellite (x-y)');
grid on;

% XZ Plot

figure;
plot(x, z);
xlabel('x (m)');
ylabel('z (m)');
title('Position of the Satellite (x-z)');
grid on;

% YZ Plot

figure;
plot(y, z);
xlabel('y (m)');
ylabel('z (m)');
title('Position of the Satellite (y-z)');
grid on;

% Altitude Plot

figure;
plot(dt, alt)
xlabel('Time (s)');
ylabel('Altitude (m)');
title('Altitude of the Satellite');
grid on;

% Longitude-Latitude Plot

figure;
plot(long, lat)
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees)');
title('Longitude-Latitude');
grid on;

% Magnetic Field "igrfmagm" ICE Plot 

figure;
plot(dt, B_ECI(1, :));
hold on;
plot(dt, B_ECI(2, :));
plot(dt, B_ECI(3, :));
hold off;
title('Magnetic Field Components ECI "igrfmagm"');
xlabel('Time (s)');
ylabel('B (nT)');
legend('Magnetic Field (x)', 'Magnetic Field (y)', 'Magnetic Field (z)');
grid on;

% Magnetic Field "wrldmagm" ICE Plot 

figure;
plot(dt, B_ECI2(1, :));
hold on;
plot(dt, B_ECI2(2, :));
plot(dt, B_ECI2(3, :));
hold off;
title('Magnetic Field Components ECI "wrldmagm"');
xlabel('Time (s)');
ylabel('B (nT)');
legend('Magnetic Field (x)', 'Magnetic Field (y)', 'Magnetic Field (z)');
grid on;

% Magnetic Field Differences Plot

figure;
plot(dt, dif_x);
hold on;
plot(dt, dif_y);
plot(dt, dif_z);
hold off;
title('Difference in Two Methods');
xlabel('Time (s)');
ylabel('B (nT)');
legend('Magnetic Field (x)', 'Magnetic Field (y)', 'Magnetic Field (z)');
grid on;


















