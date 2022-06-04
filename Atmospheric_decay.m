%% Pre-run Cleaning

clear 
clear global x
clc

%% Variable Definition

mu = 3.986004418*10^14;
Re = 6371000; % [m]
h = 380000; % [m]
J2 = 1082.63*10^-6; %[-] J2 effect parameter
a = Re + h;
e = 0;
i = deg2rad(97);
RAAN = 0;
true_anomaly = 0;
w = 0;
T = 2*pi*sqrt(a^3/mu); % [s] orbital period
Cd = 2.2; % [-] 
A = 0.0216; % [m^2] cross-sectional area
m = 8; % [m] mass


%% Initial Conditions

initial_time = [2017 3 21 11 00 00]; %real time at which the mission starts
%dyear = decyear(initial_time(1), initial_time(2), initial_time(3)); %decimal year at which the mission starts

position = kepler_to_cartesian(a, e, w, true_anomaly, i, RAAN);
velocity = calc_vel(a, e, w, true_anomaly, i, RAAN, position, mu);

sim_days = 30;
t_span = [0 sim_days*3600*24];
y0 = [position' velocity'];


%% Simulation

opts = odeset('RelTol', 1e-8);
[t, y] = ode78(@(t, y) odefunc(t, y, initial_time, Cd, A, mu, Re, J2, m), t_span, y0, opts);


%% Plotting

radii = vecnorm(y(:, 1:3), 2, 2);
heights = radii - Re;

plot(t/24/3600, heights/1000, t/24/3600, 100*ones(length(heights)), 'r');
xlabel('Time [days]')
ylabel('Height [km]')
grid on
yticks(0:50:h/1000 + 100)
ylim([0, h/1000 + 100])


%% Function Definition

function set_global_value(val)
% sets a global variable (x) a value (val)
global x
x = val;
end

function [val] = get_global_value
% gets the value of a global variable
global x
val = x;
end

function [position] = kepler_to_cartesian(a, e, w, true_anomaly, i, RAAN)
% Convert a position in the Keplerian system to a cartesian system
p = a * (1 - e^2);
r = p / (1 + e * cos(true_anomaly)); % radius

% Compute the Cartesian position vector
X = r * (cos(RAAN) * cos(w + true_anomaly) - sin(RAAN) * sin(w + true_anomaly) * cos(i));
Y = r * (sin(RAAN) * cos(w + true_anomaly) + cos(RAAN) * sin(w + true_anomaly) * cos(i));
Z = r * (sin(i) * sin(w + true_anomaly));

position = [X; Y; Z];
end

function [velocity] = calc_vel(a, e, w, true_anomaly, i, RAAN, position, mu)
% Get the velocity in cartesian coordinates
p = a * (1 - e^2);
r = p / (1 + e * cos(true_anomaly)); % radius
h = sqrt(mu * p);

V_X = (position(1) * h * e / (r * p)) * sin(true_anomaly) - (h / r) * (cos(RAAN) * sin(w + true_anomaly) + sin(RAAN) * cos(w + true_anomaly) * cos(i));
V_Y = (position(2) * h * e / (r * p)) * sin(true_anomaly) - (h / r) * (sin(RAAN) * sin(w + true_anomaly) - cos(RAAN) * cos(w + true_anomaly) * cos(i));
V_Z = (position(3) * h * e / (r * p)) * sin(true_anomaly) + (h / r) * (cos(w + true_anomaly) * sin(i));

velocity = [V_X; V_Y; V_Z];
end

function [time] = real_time(t, initial_time)
%computes the real Earth time from the mission time in seconds
time = initial_time;
remainder = time(4) + time(5)*60 + time(6)*3600;
doy = time(3) + remainder/24 + t/3600/24; % decimal day of the year

[yy, mm, dd, hh, min, ss] = datevec(datenum(time(1), time(2), doy));
time(1) = yy;
time(2) = mm;
time(3) = dd;
time(4) = hh;
time(5) = min;
time(6) = ss;
end

function [latitude, longitude, altitude] = lla(t, initial_time, position)
%computes the latitude, longitude and altitude of the spacecraft from
%orbital values and the time
%ECI position in meters
utc = real_time(t, initial_time); %Universal Coordinated Time
lla = eci2lla(position, utc);
latitude = lla(1);
longitude = lla(2);
altitude = lla(3);
end

function [rho] = density_calc(t, initial_time, position)
%computes the local density using the atmosnrlmsise00 model
%obtain real time in the form requested by atmosnrlmsise00
time = real_time(t, initial_time);
year = time(1);
d = datetime(time(1), time(2), time(3));
day_of_year = day(d, 'dayofyear');
UTseconds = (time(4)*60 + time(5))*60 + time(6);

% obtain the location
[latitude, longitude, altitude] = lla(t, initial_time, position);

% get solar flux and magnetic index data
[f107a, f107d] = getf107_func(year, day_of_year, false);
[magnetic_index] = getAPH_func(year, day_of_year, UTseconds, false);

% estimate mass density, disregarding the temperatures
[~, rhos] = atmosnrlmsise00(altitude, latitude, longitude, year, day_of_year, ...
                       UTseconds, f107a, f107d, magnetic_index); 
rho = rhos(6); % get the total mass density
end

function [a] = a_earth(position, mu)
a = -mu/(vecnorm(position))^2 *position/vecnorm(position);
end

function [a] = a_J2(position, mu, Re, J2)
x = position(1);
y = position(2);
z = position(3);
a_x = -mu*x/norm(position)^3*J2*3/2*(Re/norm(position))^2*(5*z^2/norm(position)^2-1);
a_y = -mu*y/norm(position)^3*J2*3/2*(Re/norm(position))^2*(5*z^2/norm(position)^2-1);
a_z = -mu*z/norm(position)^3*J2*3/2*(Re/norm(position))^2*(3-5*z^2/norm(position)^2);
a = [a_x; a_y; a_z];
end

function [a] = a_drag(Cd, v, A, m, rho)
a = -Cd/2*rho*norm(v)*A*v/m;
end

function [dydt] = odefunc(t, y, initial_time, Cd, A, mu, Re, J2, m)
position = y(1:3);
v = y(4:6);

% Density
if mod(floor(t), 60) == 0
    rho = density_calc(t, initial_time, position');
    set_global_value(rho);
end

rho = get_global_value;

a1 = a_earth(position, mu);
a2 = a_J2(position, mu, Re, J2);
a3 = a_drag(Cd, v, A, m, rho);
a = a1 + a3;

dydt = zeros(6, 1);

dydt(1) = v(1);
dydt(2) = v(2);
dydt(3) = v(3);
dydt(4) = a(1);
dydt(5) = a(2);
dydt(6) = a(3);

% Display time in console at an appropriate speed
if mod(floor(t), 1200) == 0
    disp(t/3600/24)
end
end
