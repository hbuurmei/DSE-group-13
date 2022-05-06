import numpy as np

n_particles = 45000  # evenly distributed across the orbit
r_earth = 6378
altitude = 800
mu = 369800.44
period = 2*np.pi*np.sqrt((r_earth+altitude)**3/mu)
i_scatter = 1*np.pi/180
h_scatter = 50
v_scatter = np.sin(i_scatter)*(r_earth+altitude)
area = h_scatter*v_scatter
frequency = n_particles/period
print('Frequency: ' + str(frequency) + ' debris fragments per second')
flux = frequency/(area*10**6)*10  # Assume factor of 10 to account for concentration of debris around the original orbit
print('The flux at the pinch points would be of: ' + str(flux) + ' debris fragments per second per square meter')
