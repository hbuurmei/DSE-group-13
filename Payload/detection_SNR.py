import matplotlib.pyplot as plt
import numpy as np

# Nomenclature:
# I = irradiance [W/m^2]
# P = power [W]

D_debris = 0.01  # [m]
A_debris = np.pi*((D_debris/2)**2)  # [m^2]

step = 1e-4
fov = np.arange(0.1, 50+step, step) * np.pi / 180  # [rad]

lambda_laser = 532e-9  # laser wavelength [m]
B_laser = 1e-9  # laser bandwidth [m]
tau_laser = 20e-15  # laser pulse duration [s]
E_pulse = 2  # total pulse energy [J]
f_max = 250  # maximum repetition rate (based on EM wave travel time for 2x500 km) [Hz]
P_pulse_max = E_pulse/tau_laser  # maximum possible pulse power, occurring during the time of tau_laser [W]
P_avg = E_pulse*f_max  # average power [W]

I_background = (B_laser/lambda_laser)*10e-9*(2*np.pi*(1-np.cos(fov/2)))  # [W/m^2] extragalactic background irradiance
SNR = 2
refl = 0.20  # reflection of debris
eff_t = 0.70  # transmitter optics efficiency
eff_r = 0.70  # receiver optics efficiency
eff_q = 0.40  # quantum efficiency of detector

# ----------------------- POWER VERSUS FOV, FOR A SET OF DISTANCES --------------------------
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.set_title(f'Active Detection (Debris Diameter: {D_debris*100:.1f} cm)')
ax1.set_xlabel('FoV [deg]')
ax1.set_ylabel('Pulse Power [W]')

for distance in [50000, 100000, 200000, 300000, 400000, 500000]:
    A_hemi = 2*np.pi*(distance**2)  # [m^2] area of hemisphere, caused by light reflection off debris
    P_debris = SNR*I_background*A_hemi/(eff_r*eff_q)  # [W] required total power over the hemisphere at the S/C, for SNR 2
    I_debris = (P_debris/A_debris)/refl  # [W/m^2] required power per area at the debris
    A_cone = 2*np.pi*(distance**2)*(1-np.cos(fov/2))  # [m^2] cone area illuminated by FoV
    P_pulse = I_debris*A_cone/eff_t   # [W] total power sent

    ax1.semilogy(fov*180/np.pi, P_pulse, label=f'Distance = {int(distance/1000)} km')

    # Calculate maximum FoV for the given max pulse power, and max scanning diameter at the given distance.
    fov_max = fov[P_pulse < P_pulse_max][-1]  # [rad]
    D_cone_max = 2*np.tan(fov_max/2)*distance  # [m]

    print(f'---------- Distance = {int(distance/1000)} km ----------')
    print(f'FoV at max power: {fov_max*180/np.pi} deg')
    print(f'Searching diameter at the debris: {D_cone_max/1000} km\n')

ax1.hlines(P_pulse_max, xmin=fov[0]*180/np.pi, xmax=fov[-1]*180/np.pi, colors='r', linestyles='dashed',
           label=r'$P_{pulse, max}$')
ax1.grid()
ax1.legend()

# -------------------- FOV AND SEARCHING DIAMETER VERSUS DISTANCE, FOR A GIVEN MAX POWER -----------------------
# Derivation of the formula is given in the Final Report of the project.
distance = ((P_pulse_max*eff_t*A_debris*refl*eff_r*eff_q)/(4*SNR*np.pi*np.pi*(B_laser/lambda_laser)*10e-9*2*np.pi*(1-np.cos(fov/2))**2))**(1/4)
D_search = 2*distance*np.tan(fov/2)  # [m^2] searching area illuminated by FoV
fig2 = plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)
ax2.semilogx(distance/1000, fov*180/np.pi, color='tab:blue')
ax2.set_title(f'Active Detection (Debris Diameter: {D_debris*100:.1f} cm)')
ax2.set_ylabel('FoV [deg]', color='tab:blue')
ax2.set_xlabel('Distance [km]')
ax2.tick_params(axis='y', labelcolor='tab:blue')

ax3 = ax2.twinx()
ax3.semilogx(distance/1000, D_search/1000, color='tab:orange')
ax3.set_ylabel('Searching Diameter [km]', color='tab:orange')
ax3.tick_params(axis='y', labelcolor='tab:orange')

plt.show()
