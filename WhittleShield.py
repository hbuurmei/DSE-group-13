import numpy as np

# Based on: https://ntrs.nasa.gov/api/citations/19920010826/downloads/19920010826.pdf?attachment=true
# More detail can be added on later by looking at: https://ntrs.nasa.gov/api/citations/19920010785/downloads/19920010785.pdf?attachment=true



# Assuming a square based triangle, with the base as the entrance to the net:
# length of a side [m]
side_length = 1000
# height of triangle [m]
depth = 50
tan_side = np.sqrt((side_length/2)**2 + depth**2)
# Interior surface area
int_area = (0.5 * side_length * tan_side) * 4



# Projectile Diamater [cm], loses accuracy around 2.4 [cm] #
d = 2.4
# Density of debris projectile [g/cc], range given between 1.14 and 2.8 in report mentioned above.
delta_p = 2.8
# Density of bumper material [g/cc], range between 2.7 and 2.8
delta_b = 2.8
# Debris projectile mass [g]
M = 10
# Relative velocity [km/s]
v = 0.5
# Geometry angle between horizontal and net [rad]
theta = np.arctan(((side_length/2) / depth) * np.pi/180)
# Spacing between bumper and rear wall [cm], S/d between 15 and 96
S = 15 * d
# Yield strength of material used [ksi], range between 18 and 70
sigma = 70
# C constant
C = 0.16



# Bumper thickness
t_b = 0.25 * d * delta_p / delta_b

# Rear wall thickness
t_w = C * (d**0.5) * ((delta_p * delta_b)**(1/6)) * (M**(1/3)) * (((v*np.cos(theta))/S)**0.5) * ((70/sigma)**0.5)

# Critical particle diameter causing failure
d_crit = 3.918 * (t_w**(2/3)) * ((delta_p)**(-1/3)) * ((delta_b)**(-1/9)) * (v**(-2/3)) * ((np.cos(theta))**(-2/3)) * (S**(1/3)) * ((sigma/70)**(1/3))

# Mass of shield = combined thicknesses * internal surface area * density [kg]
mass = ((t_b + t_w) * delta_b * (int_area * 10000))/1000

print('Debris size is', d, 'cm')
print('Bumper thickness is', t_b, 'cm')
print('Real wall thickness is', t_w, 'cm')
print('Critical particle diameter is', d_crit, 'cm')
print('The mass of the shield is', mass, 'kg')

if t_b/d < 0.15:
    print('The t_b/d ratio is less than 0.15 (t_b/d =', t_b/d, '). The results are being underestimated.' )

if v*np.cos(theta) < 7:
    print('The normal velocity is less than 7 km/s (v_n =', v * np.cos(theta), '). The results are being underestimated.')

if S/d < 15:
    print('The ratio of S/d (spacing between bumper and rear wall to size of debris) is less than 15 (S/d =', S/d, '). The results are being underestimated.')












