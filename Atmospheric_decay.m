%% Pre-run Cleaning

clear
clear global x
clc

%% Variable Definition

mu = 3.986004418*10^14;
Re = 6371000; % [m]
J2 = 1082.63*10^-6; %[-] J2 effect parameter#


%% Simulation

% Stationkeeping
%{
global initial_time;
initial_time = [2014 3 21 11 00 00]; %real time at which the mission starts
h_p = 380*10^3;
h_a = 380*10^3;

a = (2*Re + h_p + h_a) / 2;
e = (h_a - h_p) / (2*a);
i = deg2rad(97);
RAAN = 0;
true_anomaly = 0;
w = 0;
T = 2*pi*sqrt(a^3/mu); % [s] orbital period
Cd = 3; % [-] 
A_m = 0.022; % area-to-mass ratio of the object
n_orbits = 3;
t_span = [0 n_orbits*T]; % 1 orbit simulated, the more orbits, the worse this approximation
position = kepler_to_cartesian(a, e, w, true_anomaly, i, RAAN);
velocity = calc_vel(a, e, w, true_anomaly, i, RAAN, position, mu);
y0 = [position' velocity'];
opts = odeset('RelTol', 1e-12);
[t1, y1] = ode78(@(t, y) odefunc(t, y, initial_time, Cd, A_m, mu, Re, J2, true, false), t_span, y0, opts); % With J2, without drag
[t2, y2] = ode78(@(t, y) odefunc(t, y, initial_time, Cd, A_m, mu, Re, J2, true, true), t_span, y0, opts); % With J2, with drag
deltav_stationkeeping = abs(norm(y1(end)) - norm(y2(end)))
deltav_over_one_year = deltav_stationkeeping * 265 * 24 * 3600 / T
%}

% Debris decay
global initial_time;
initial_time = [2006 8 21 11 00 00]; %real time at which the mission starts
h_a = 1000*10^3;
i = deg2rad(97);
RAAN = 0;
true_anomaly = 0;
w = 0;
Cd = 2.2; % [-]
A_m = 0.07946; % area-to-mass ratio of the object
sim_days = 730;
t_span = [0 sim_days*3600*24];
h_p_range = 350; % [km]
opts = odeset('RelTol', 1e-6, 'Events', @events);
decay_times = zeros(size(h_p_range));
global min_altitudes;
global t_min_altitudes;
t_min_altitudes = 0;
for idx = 1:numel(h_p_range)
    h_p = h_p_range(idx) * 10^3;
    min_altitudes = h_p;
    a = (2*Re + h_p + h_a) / 2;
    e = (h_a - h_p) / (2*a);
    position = kepler_to_cartesian(a, e, w, true_anomaly, i, RAAN);
    velocity = calc_vel(a, e, w, true_anomaly, i, RAAN, position, mu);
    y0 = [position' velocity'];
    [t, y] = ode78(@(t, y) odefunc(t, y, initial_time, Cd, A_m, mu, Re, J2, false, true), t_span, y0, opts); % Without J2
    decay_times(idx) = t(end);
    T = table(t_min_altitudes', min_altitudes');
    writetable(T, string(h_p) + ".txt")
end
    


%% Plotting

%{
figure(1)
plot(h_p_range, decay_times/24/3600);
xlabel('Periapsis Altitude [km]')
ylabel('Decay time [days]')
grid on
[max_time, mmax_idx] = max(decay_times);
%ylim([0, min(sim_days,max_time) + 10])
%saveas(gcf,'Peri-vs-Alt.png')
%}


%figure(2)
radii = vecnorm(y(:, 1:3), 2, 2);
heights = radii - Re;
stride = ceil(min(sim_days, t(end)/24/3600));
plot(t(1:stride:end)/24/3600, heights(1:stride:end)/1000);
xlabel('Time [days]')
ylabel('Height [km]')
grid on
yticks(0:50:h_a/1000 + 100)
ylim([0, h_a/1000 + 100])
savefig(string(h_p_range) + ".txt")


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

function [rho] = density_calc(t, latitude, longitude, altitude)
%computes the local density using the atmosnrlmsise00 model
%obtain real time in the form requested by atmosnrlmsise00
global initial_time;
time = real_time(t, initial_time);
year = time(1);
d = datetime(time(1), time(2), time(3));
day_of_year = day(d, 'dayofyear');
UTseconds = (time(4)*60 + time(5))*60 + time(6);

% obtain the location
%[latitude, longitude, altitude] = lla(t, initial_time, position');

% get solar flux and magnetic index data
[f107a, f107d] = getf107_func(year, day_of_year, false);
[magnetic_index] = getAPH_func(year, day_of_year, UTseconds, false);

if (altitude < 0) | (altitude > 1000*10^3)
    altitude
end
clamp_alt = max(min(1000*10^3,altitude),0);

% estimate mass density, disregarding the temperatures
[~, rhos] = atmosnrlmsise00(clamp_alt, latitude, longitude, year, day_of_year, ...
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

function [a] = a_drag(Cd, v, A_m, rho)
a = -Cd/2*rho*norm(v)*A_m*v;
end

function [value,isterminal,direction] = events(t,y)
% Locate the time when y passes through 0.111 in all 
% directions and stop integration.
global initial_time;
position = y(1:3);

[latitude, longitude, altitude] = lla(t, initial_time, position');
value = (altitude > 0);
isterminal = 1; % Stop the integration
direction = -1; % When altitude decreases and crosses 0
end


function [dydt] = odefunc(t, y, initial_time, Cd, A_m, mu, Re, J2, j2_enabled, drag_enabled)
position = y(1:3);
v = y(4:6);

% obtain the location
[latitude, longitude, altitude] = lla(t, initial_time, position');

global min_altitudes
global t_min_altitudes
if (t > t_min_altitudes(end))
    min_altitudes(end+1) = min(min_altitudes(end), altitude);
    t_min_altitudes(end+1) = t;
end

% Density
if mod(floor(t), 5) == 0
    rho = density_calc(t, latitude, longitude, altitude);
    set_global_value(rho);
end

rho = get_global_value;

a1 = a_earth(position, mu);
a2 = a_J2(position, mu, Re, J2);
a3 = a_drag(Cd, v, A_m, rho);
a = a1 + a2 * j2_enabled + a3 * drag_enabled;

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