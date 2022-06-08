import matplotlib.pyplot as plt
import numpy as np

# ---------------------------- CONSTANTS & GLOBAL PARAMETERS -------------------------------
# Debris characteristics
D_debris = 0.01  # [m] debris diameter
A_debris = np.pi*((D_debris/2)**2)  # [m^2] debris frontal area, assuming a circle
refl = 0.20  # Lambertian reflection of debris

# Lidar characteristics
lambda_laser = 532e-9  # [m] laser wavelength
B_sensor = 0.17e-9  # [m] sensor bandwidth]
tau_laser = 1e-9  # [s] laser pulse duration
E_pulse = 10e-3  # [J] total pulse energy
f_rep = 10000  # [Hz] repetition rate
P_pulse_max = E_pulse/tau_laser  # [W] maximum possible pulse power, occurring during the time of tau_laser
P_avg = E_pulse*f_rep  # [W] average power
eff_t = 0.70  # transmitter optics efficiency
eff_r = 0.70  # receiver optics efficiency
eff_q = 0.40  # quantum efficiency of detector
overlap = 0.3  # percentage of overlap between the images

# Creating a fov array, and calculating the background noise
step = 1e-4  # step size for fov [deg]
fov = np.arange(0.01, 20+step, step)*np.pi/180  # [rad] fov, using [deg] as input in np.arange(), for easier use
I_background = (B_sensor/lambda_laser)*10e-9*(2*np.pi*(1-np.cos(fov/2)))  # [W/m^2] extragalactic background irradiance
SNR = 2  # Signal-to-noise ratio

# ----------------------- POWER VERSUS FOV, FOR A SET OF DISTANCES --------------------------
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.set_xlabel('FoV [deg]')
ax1.set_ylabel('Pulse Power [W]')

# Use different linestyles to make it easier to read in black and white.
linestyles = [(0,(1,1)), (0,(1,4)), (0,(5,6)), (0,(5,2)), (0,(3,5,1,5)), (0,(3,4,1,4,1,4))]

x = 0  # Use a counter to iterate over the different linestyles
for distance in [50000, 100000, 200000, 300000, 400000, 500000]:
    # Explanation of calculations is shown in the Final Report.
    A_hemi = 2*np.pi*(distance**2)  # [m^2] area of hemisphere, caused by light reflection off debris
    P_debris = SNR*I_background*A_hemi/(eff_r*eff_q)  # [W] required total power over the hemisphere
    I_debris = (P_debris/A_debris)/refl  # [W/m^2] required power per area at the debris
    A_cone = 2*np.pi*(distance**2)*(1-np.cos(fov/2))  # [m^2] cone area illuminated by FoV
    P_pulse = I_debris*A_cone/eff_t   # [W] total power sent by the laser

    ax1.loglog(fov*180/np.pi, P_pulse, label=f'Distance = {int(distance/1000)} km', linestyle=linestyles[x])

    x += 1

ax1.hlines(P_pulse_max, xmin=fov[0]*180/np.pi, xmax=fov[-1]*180/np.pi, colors='r', label=r'$P_{pulse, max}$')
ax1.grid()
ax1.legend()

# --------------- FOV VERSUS DISTANCE AND SCANNING TIME, FOR A GIVEN MAX POWER ------------------
# Explicit relation for the distance, as a function of FoV, for the given max power
distance = ((P_pulse_max*eff_t*A_debris*refl*eff_r*eff_q)/(4*SNR*np.pi*np.pi*I_background*(1-np.cos(fov/2))))**(1/4)
A_cone = 2*np.pi*(distance**2)*(1-np.cos(fov/2))  # should replace A_search

fig2 = plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)

# Including a scanning pattern to increase the FoV.
x = 0
for t_pattern in [1/f_rep, 0.1, 1, 5, 10, 20]:  # [s] time for scanning pattern (first element = stationary FoV)
    '''If the total scanning area is larger than the entire sphere at the distance, the argument to be entered in 
    arccos exceeds the limit. So, only when abs(cos_arg) <= 1, they are entered in arccos. The other values
    are left as Not-a-Number.'''
    A_pattern = t_pattern*f_rep*A_cone*(1-overlap)  # more overlap between images reduces possible FoV
    cos_arg = 1-A_pattern/(2*np.pi*(distance**2))  # argument to be entered in arccos
    fov_pattern = np.array([np.nan]*len(distance))  # initiate empty array
    fov_pattern[np.abs(cos_arg) <= 1] = 2*np.arccos(cos_arg[np.abs(cos_arg) <= 1])  # check for compliance with arccos

    if x == 0:  # separate 'stationary' FoV from scanning FoV, to create a more logical label for the legend
        ax2.semilogx(distance/1000, fov_pattern*180/np.pi, linestyle=linestyles[x], label=f'Stationary (no scanning)')
    else:
        ax2.semilogx(distance/1000, fov_pattern*180/np.pi, linestyle=linestyles[x], label=f'{t_pattern:.2f} [s] Scanning')

    x += 1

ax2.set_ylabel('Total FoV [deg]')
ax2.set_xlabel('Distance [km]')
ax2.legend()
ax2.grid()

plt.show()
