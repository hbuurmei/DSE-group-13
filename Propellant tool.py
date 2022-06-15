import math
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean

mission_duration = 1 #year
h = 380 #CHANGE#
G = 6.67*10**-11
m_E = 5.972*10**24 
r_E = 6371 #km
mu = G * m_E
a = (r_E + h)*1000
g0 = 9.80665
v_orbit = np.sqrt(mu/a)
period = 2* np.pi * np.sqrt((a**3)/(mu))
burn_year = (mission_duration * 365.25 * 24 * 3600) / period #number of orbits per year, thus number of burns per year

d_inclination = 121.5*np.pi/180 #deg
dV_inclination_worst_case = np.sqrt(2* v_orbit**2 * (1 - np.cos(d_inclination)))

#Inputs for orbit insertion and EOL
dV_insert = 105 #m/s
dV_EOL = 81.8 #m/s

#Inputs for drag compensation (from tool)
dV_tot = 335.8869 #m/s
dV_burn = 0.0587 #m/s per burn (per orbit)

#Thruster characteristics
M_thruster = 1.12 #kg
Isp = 234 #s
m_dot = 0.0506  #kg/s
density_prop = 982 #kg/m3 for Hydrazine
p_min = 1620270 #Pa

#Pressurant characteristics
B = 3.5 #blow-down ratio
pressure_req_for_thrusting = 3.9 #MPa
density_pressurant_i =  48.9 #kg/m^3 for N2

#margins:
margin_trapped = 0.03
margin_prop_loading_uncert = 0.005
margin_tank = 0.2
margin_pressure = 0.2
margin_feed = 0.2

#1st order propellant mass estimation
M_dry = 2880 #kg #CHANGE#
M_prop = M_dry * ((np.exp(dV_tot/(Isp*g0)) - 1))
M_wet = M_dry + M_prop

#detailed propellant mass estimation
propellant_mass = []
counter = []
mass_tot = []
burntime = []


start = M_dry
M_fin_EOL = start * np.exp((dV_EOL)/(Isp * g0))

M_prop_EOL_burn = M_fin_EOL - M_dry
burntime_EOL = M_prop_EOL_burn / m_dot

burntime.append(burntime_EOL)
counter.append(0)
counter.append(1)

propellant_mass.append(M_prop_EOL_burn)
mass_tot.append(M_dry)
mass_tot.append(M_fin_EOL)


start = M_fin_EOL
for i in range(2, int(burn_year)+2, 1):
    #print(i)
    counter.append(i)
    fin = start * np.exp((dV_burn)/(Isp * g0))
    mass_tot.append(fin)
    m_prop_burn = fin - start
    propellant_mass.append(m_prop_burn)
    burntime_pulse = m_prop_burn / m_dot
    burntime.append(burntime_pulse)
    start = fin

M_fin_insertion = start * np.exp((dV_insert)/(Isp * g0))
M_prop_insertion = M_fin_insertion - start
burntime_insert = M_prop_insertion / m_dot

burntime.append(burntime_insert)
propellant_mass.append(M_prop_insertion)
counter.append(int(burn_year)+3)
mass_tot.append(M_fin_insertion)

mass_tot.reverse()
propellant_mass.reverse()
burntime.reverse()


#print(mean(mass_tot)) average mass of spacecraft


# plt.plot(counter[2:-1], burntime[1:-1], label = "MONARC-445")
# plt.xlabel("Pulses")
# plt.ylabel("Burntime [s]")
# plt.legend()
# plt.show()

plt.plot(counter, mass_tot, label = "MONARC-90HT")
plt.xlabel("Pulses")
plt.ylabel("Total Spacecraft Mass [kg]")
plt.legend()
plt.show()

total_propellant_mass_required = mass_tot[0]- M_dry #total propellant mass in kg

Mp_loaded = total_propellant_mass_required * (1 + margin_trapped + margin_prop_loading_uncert)
Vp_loaded = Mp_loaded / density_prop
Vp_usable = total_propellant_mass_required / density_prop

print("Propellant Volume (including margins for trapped loading uncertainty) = ", round(Vp_loaded,2), "m3")

#pressurant tank for blowdown system:
V_ullage_i = Vp_usable / (B - 1)

#final tank volume
V_tank = (1 + margin_tank) * (Vp_loaded + V_ullage_i)
V_tank_liter = 1000* V_tank



#Tank mass relationships:
counter1 = []
pmd_tank = []
prop_volume_lst = []
diaphragm_tank = []
for j in range(0, 750, 1):
    counter1.append(j)
    prop_volume_lst.append(V_tank_liter)
    PMD = (2.7086 * (10**(-8)) * j**3) - (6.1703 * (10**(-5)) * j**2 )+ (6.6290 * (10**(-2)) * j) + 1.3192
    Diaphragm = (2.36 * (10**(-7)) * j**3) - (2.32 * (10**(-4))* j**2) + (0.131* j) + 0.264
    pmd_tank.append(PMD)
    diaphragm_tank.append(Diaphragm)
M_tank_PMD = (2.7086 * (10**(-8)) * V_tank_liter**3) - (6.1703 * (10**(-5)) * V_tank_liter**2 )+ (6.6290 * (10**(-2)) * V_tank_liter) + 1.3192
M_tank_Diaphragm = (2.36 * (10**(-7)) * V_tank_liter**3) - (2.32 * (10**(-4))* V_tank_liter**2) + (0.131* V_tank_liter) + 0.264

# plt.plot(prop_volume_lst[:100], counter1[:100],label = "Required Propellant Volume")
# plt.plot(counter1, pmd_tank, label = "PMD Tank" )
# plt.plot(counter1, diaphragm_tank, label = "Diaphragm Tank")
# plt.xlabel("Volume tank [liters]")
# plt.ylabel("Tank Mass [kg]")
# plt.legend()
# plt.show()
#print("Mass PMD tank = ", M_tank_PMD)
#print("Tank mass Diaphragm tank = ", M_tank_Diaphragm)


M_pressurant = V_ullage_i * density_pressurant_i


M_feed = ((M_thruster + M_tank_PMD) / (1 - margin_feed)) * margin_feed



print("========================================================================")
print("Thruster mass = ", M_thruster, "kg")
print("Tank mass = ", round(M_tank_PMD,3), "kg")
print("Feed system = ", round(M_feed,3), "kg")
print("========================================================================")
print("Total Dry mass = ", M_thruster + M_tank_PMD + M_feed, "kg")
print("========================================================================")
print("Propellant Mass = ", round(Mp_loaded,3), "kg")
print("Pressurant Mass = ", round(M_pressurant,3), "kg")
print("========================================================================")
print("Total Mass of Propulsion Subsystem = ", M_thruster + M_tank_Diaphragm + M_feed + Mp_loaded + M_pressurant, "kg")
print("========================================================================")

#Pressurant tank sizing
V_pressuranttank = 0.1
r = ((3*V_pressuranttank)/(4*np.pi))**(1/3)
V_tot = V_tank + V_pressuranttank
p_pressure_tank_BOL = (p_min * V_tot) / V_pressuranttank

sigma_material_pressuretank = 503 #MPa ~Aluminum 7075-T6
density_material_pressuretank = 2810 #kg/m3 ~Aluminum 7075-T6
sigma_c = []
sigma_sph = []
sigma_mat = []
counter2 = []

    #creating thickness vs stress plot and calculating thickness for spherical tank
for i in range(1,11,1):
    t = i/1000
    counter2.append(t)
    sigmat = sigma_material_pressuretank * 1000000
    sigma_mat.append(sigmat)
    sigma_cyl = (p_pressure_tank_BOL * r) / (t)
    sigma_spher = (p_pressure_tank_BOL * r) / (2*t)
    sigma_c.append(sigma_cyl)
    sigma_sph.append(sigma_spher)
    t_c = (p_pressure_tank_BOL * r) / (sigmat)
    t_s = (p_pressure_tank_BOL * r) / (2* sigmat)
    
V_s = (4/3) * np.pi * ((r+t_s)**3 - r**3)
M_s = density_material_pressuretank * V_s

print(M_s)
    
# plt.plot(counter2, sigma_c, label = "Cylinder")
# plt.plot(counter2, sigma_sph, label = "Spherical")
# plt.plot(counter2, sigma_mat, label = "Chosen Material")
# plt.xlabel("thickness [m]")
# plt.ylabel("Stress [Pa]")
# plt.legend()
# plt.show()









