import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from mpl_toolkits import mplot3d
plt.ioff()  # not plot without plt.show()

main = False
MRN_Bool = True
CI_Bool = False

# For plotting in 3D
fig = plt.figure(facecolor="white")
ax = plt.axes(projection='3d')
colours = ["red", "blue", "green", "black", "pink"]


# Constants
c = 299792458  # m/s
k_bolt = 1.380649 * 1e-23  # J/K
AU = 149597871  # km
J2_val = 0.00108263  # [-]
R_Earth = 6371  # km
Earth_day_duration = 24 * 60 * 60  # s
mu_Earth = 398600.4418  # km^3/s^2
T_end = 1*24*60*60  # s mission duration
Earth_Rotation_Rate = 2 * np.pi / Earth_day_duration
elevDeg = 1
f_trans = 3 * 1e9  # S Band

# Relay satellites
#     [a, e, i, argument of pericenter, right ascension of ascending node, name]
CI_TTW = [R_Earth+6000, 0, 0, 0, 0]

# right ascension of ascending node is the longitude of the ascending node basically
class planet:
    def __init__(self, radius, mass, day_duration, mu):
        self.mass = mass
        self.radius = radius
        self.T_revolution = day_duration
        self.mu = mu
        self.omega = 2 * np.pi / day_duration

    def drawPlanet(self):  # To be used only for illustration only in the report
        '''
        Draw the planet wanted. Here in orange for Mars.
        '''
        r = self.radius
        phi = np.linspace(-np.pi/2, np.pi/2, 100)
        lamb = np.linspace(-np.pi, np.pi, 100)
        x = r * np.outer(np.cos(lamb), np.cos(phi))
        y = r * np.outer(np.sin(lamb), np.cos(phi))
        z = r * np.outer(np.ones(np.size(lamb)), np.sin(phi))
        ax.plot_surface(x,y,z, rstride=5, cstride=5, color="orange")  # plotting

class GS:
    def __init__(self, Name, planetRadius=R_Earth, latitudeDeg=0.001, longitudeDeg=0.001):
        self.lat = np.deg2rad(latitudeDeg)
        self.long = np.deg2rad(longitudeDeg)
        self.height = planetRadius
        self.name = Name
        self.dataTransmittedMRN = 0
        self.dataTransmittedTTW = 0
        self.timeAvailableMRN = 0
        self.timeAvailableTTW = 0

    def getPosition(self, centerLatitudeDeg, SDLatitude): # This should be done in the setup of the code
        '''
        Give a random position to the rover. Uniform distribution in longitude
        and normal distribution on latitude.
        @param: centerLatitudeDeg, center of the normal distribution in latitude.
        @param: SDLatitude, standard deviation of the normal distribition in latitude.
        @return: long, the longitude attributed.
        @return: lat, the latitude attributed.
        '''
        # Uniform random distribution in longitude
        longOptions = np.linspace(-np.pi, np.pi, 1000)
        longSuggested = np.random.choice(longOptions)

        # Truncated normal random distribution for latitude by removing the rovers below 20deg S and higher than 90deg N
        latSuggested = 100
        while abs(latSuggested) > np.pi/2 or latSuggested < np.radians(-20) or latSuggested > np.radians(70):
            latSuggested = np.random.normal(np.radians(centerLatitudeDeg), SDLatitude)

        # Store result - for testing, put 0 0
        self.long = longSuggested
        self.lat = latSuggested
        return self.long, self.lat

    def printInfo(self, boolName, boolPos):
        if boolName:
            print(self.name)
        if boolPos:
            print(f"latitude = {np.degrees(self.lat)}, longitude = {np.degrees(self.long)}")

    def cartesianCoords(self):
        '''
        Get Cartesian coordinates of the rover (by converting the coords).
        @return: x - y - z, the coordinates of the rover in the Cartesian coordinate system
        '''
        h = self.height
        phi = self.lat
        lamb = self.long

        # Conversion from spherical coords to cartesian coords
        x = h * np.cos(phi) * np.cos(lamb)
        y = h * np.cos(phi) * np.sin(lamb)
        z = h * np.sin(phi)

        # self.xpos, self.ypos, self.zpos = x, y, z
        return x, y, z

    def updatePosition(self, RotationRate, timeStep):
        '''
        Update the rover position in an Earth centred but non-fixed system
        '''
        self.long = self.long + timeStep * RotationRate
        self.lat = self.lat
        return self.long, self.lat

    def results(self, MarsRotationRate):  # ATTENTION 0 !!!!!!!!
        print(f"{np.degrees(self.long - T_end*MarsRotationRate)} & {np.degrees(self.lat)} & {100* self.timeAvailableTTW/(T_end)} & {self.dataTransmittedTTW/1e9} & {100 * self.timeAvailableMRN/(T_end)} & {self.dataTransmittedMRN/1e9}")

    def drawRover(self):  # This is only used once for the rendering in the report
        '''
        Draw rovers on the surface of Mars. Illustration for the report
        '''
        xpos, ypos, zpos = GS.cartesianCoords()  # Get the cartesian coordinates of each rover
        r = 50 # Make sure that it can be seen

        # Define draw surface
        phi = np.linspace(-np.pi / 2, np.pi / 2, 100)
        lamb = np.linspace(-np.pi, np.pi, 100)
        x = r * np.outer(np.cos(lamb), np.cos(phi))
        y = r * np.outer(np.sin(lamb), np.cos(phi))
        z = r * np.outer(np.ones(np.size(lamb)), np.sin(phi))

        # Plot
        ax.plot_surface(xpos-x, ypos-y, zpos-z, color="blue")  # For on planet

class relay:
    def __init__(self, semiMajorAxis, eccentricity, inclinationRad, argumentOfPericenterRad, ascendingNodeRad, timeOfPericenterPassage, sc_name, mu = mu_Earth):
        self.e = eccentricity
        self.a = semiMajorAxis
        self.t_p = timeOfPericenterPassage
        self.i = inclinationRad
        self.Omega = ascendingNodeRad
        self.w = argumentOfPericenterRad
        self.mu = mu
        self.name = sc_name
        self.commTime = 0
        self.inCommLatRec = []
        self.inCommLongRec = []
        self.noCommLatRec = []
        self.noCommLongRec = []
        self.commDuration = []

    def orbitParameters(self):
        '''
        Determine various orbital parameters
        @param: no param
        @return: T, the orbital period in seconds
                 r_p, the pericenter position in kilometers
                 r_a, the apocenter position in kilometers
        '''
        T = 2 * np.pi * np.sqrt(self.a**3 / self.mu)  # Orbital period (s)
        r_p = self.a * (1 - self.e)  # Pericenter (km)
        r_a = self.a * (1 + self.e)  # Apocenter (km)
        return T, r_p, r_a

    def getPosition(self, t):
        '''
        Find the position of the spacecraft in the Keplerian system.
        @param: t, the time variable in seconds
        @return: true_anomaly, the true anomaly of the spacecraft at time t
        '''
        e = self.e
        n = np.sqrt(self.mu/self.a**3) # Mean motion
        M = n * (t-self.t_p) # Mean anomaly

        # Newton's method to solve the equation
        def eq(E, e, M): return E - e *np.sin(E) - M
        def d_eq(E, e): return 1 - e * np.cos(E)
        err, E_old = 100, M
        while abs(err) > 1e-13: # Make smaller late if possible
            #print(E_old)
            E_new = E_old - (eq(E_old, e, M)/d_eq(E_old, e))
            err = E_new - E_old
            E_old = E_new
        E = E_new  # Choose final value

        # Final equation
        true_anomaly = 2 * np.arctan(np.sqrt((1+e)/(1-e)) * np.tan(E/2))
        return true_anomaly

    def KeplerToCartesian(self, true_anomaly):
        '''
        Convert the a position in the Keplerian system to a cartesian system
        @param: true_anomaly, the true anomaly of the spacecraft at time t
        '''
        p = self.a * (1-self.e**2)
        r = (p)/(1 + self.e * np.cos(true_anomaly))  # radius
        h = np.sqrt(self.mu * p)  # Specific angular momentum

        # Compute the Cartesian position vector
        X = r * (np.cos(self.Omega) * np.cos(self.w + true_anomaly) - np.sin(self.Omega)*np.sin(self.w + true_anomaly) * np.cos(self.i))
        Y = r * (np.sin(self.Omega) * np.cos(self.w + true_anomaly) + np.cos(self.Omega) * np.sin(self.w + true_anomaly) * np.cos(self.i))
        Z = r * (np.sin(self.i) * np.sin(self.w + true_anomaly))

        # Compute the Cartesian velocity vector
        X_dot = ((X * h * self.e)/(r * p)) * np.sin(true_anomaly) - (h/r) * (np.cos(self.Omega) * np.sin(self.w + true_anomaly) + np.sin(self.Omega) * np.cos(self.w + true_anomaly) * np.cos(self.i))
        Y_dot = ((Y * h * self.e)/(r * p)) * np.sin(true_anomaly) - (h/r) * (np.sin(self.Omega) * np.sin(self.w + true_anomaly) - np.cos(self.Omega) * np.cos(self.w + true_anomaly) * np.cos(self.i))
        Z_dot = ((Z * h * self.e)/(r * p)) * np.sin(true_anomaly) + (h/r) * np.sin(self.i) * np.cos(self.w + true_anomaly)

        # # Mars Centred and Fixed coordinate system. THIS
        # X_MCF = X_dot + rotationRate * Y
        # Y_MCF = Y_dot - rotationRate * X
        # Z_MCF = Z_dot

        return X, Y, Z

    def groundTrackData(self, X, Y, Z, t, EarthRot, bool):
        t0 = 0
        if t*EarthRot >= 2*np.pi:
            t0 = t
        lat_t = np.arcsin(Z/np.sqrt(X**2 + Y**2 + Z**2))
        long_t = np.arctan2(Y,X) - (t-t0)*EarthRot
        if long_t < -np.pi:
            long_t = 2 * np.pi + long_t
        if bool:
            self.inCommLatRec.append(np.rad2deg(lat_t))
            self.inCommLongRec.append(np.rad2deg(long_t))
            self.commDuration.append(1)
        else:
            self.noCommLatRec.append(np.rad2deg(lat_t))
            self.noCommLongRec.append(np.rad2deg(long_t))
            self.commDuration.append(0)
        return np.rad2deg(lat_t), np.rad2deg(long_t)

    def J2(self, dt, mu=mu_Earth, R_e=R_Earth, J_2= J2_val):
        n = np.sqrt(mu / self.a ** 3)
        RAAN_dot = -1.5 * n * R_e ** 2 * J_2 * np.cos(np.deg2rad(self.i)) / self.a ** 2 / (1 - self.e ** 2) ** 2
        w_dot = 0.75 * n * R_e ** 2 * J_2 * (4 - 5 * (np.sin(np.deg2rad(self.i))) ** 2) / self.a ** 2 / (1 - self.e ** 2) ** 2
        self.Omega = self.Omega + RAAN_dot*dt
        self.w = self.w + w_dot*dt
        return False


    def inCommunicationCone(self, relayPosition, roverPosition, elevationAngleDeg):  # THIS ONE WILL STILL NEED SOME DEBUGGING
        '''
        Determines whether the satellite can be in communication with the rover or not.
        @param: relayPosition, position of the relay satellite under consideration [x, y, z]
        @param: roverPosition, position of the rover under consideration [x, y, z]
        @return: Comm, boolean stating whether telecommunication is possible
        @return: r_relayRovSys, the distance between the rover and the relay satellite
        '''
        Comm = False
        # Decompose basic information
        relayPosition = np.array(relayPosition)
        X_relay, Y_relay, Z_relay = relayPosition[0], relayPosition[1], relayPosition[2]
        X_rover, Y_rover, Z_rover = roverPosition[0], roverPosition[1], roverPosition[2]
        r_rover = np.sqrt(X_rover**2 + Y_rover**2 + Z_rover**2)
        k = np.array([0, 0, 1])  # Unit vector along Z-axis
        # print(np.sqrt((X_relay-X_rover)**2 + (Y_relay-Y_rover)**2 + (Z_relay-Z_rover)**2))  # OK but diverges a little bit somehow
        # Determine rotation angle
        theta_rot = np.arccos(np.dot(roverPosition, k)/r_rover)

        # Unit axis of rotation
        #print(roverPosition, k)
        #print(np.sqrt(X_rover**2 + Y_rover**2+ Z_rover**2))
        U = np.cross(roverPosition, k)#/np.abs(np.cross(roverPosition, k))
        u = U /np.sqrt(U[0]**2 + U[1]**2 + U[2]**2)
        ux, uy, uz = u[0], u[1], u[2]
        #print(u)
        # Rotation Matrix by theta_rot along u.
        R_theta = np.matrix([[np.cos(theta_rot) + ux**2 * (1-np.cos(theta_rot)), ux * uy * (1-np.cos(theta_rot))-uz*np.sin(theta_rot), ux*uz*(1-np.cos(theta_rot)) + uy * np.sin(theta_rot)],
                   [uy*ux*(1-np.cos(theta_rot)) + uz * np.sin(theta_rot), np.cos(theta_rot) + uy**2 * (1-np.cos(theta_rot)), uy*uz*(1-np.cos(theta_rot))-uz*np.sin(theta_rot)],
                   [uz*ux*(1-np.cos(theta_rot))-uy*np.sin(theta_rot), uz*uy*(1-np.cos(theta_rot))+ux*np.sin(theta_rot), np.cos(theta_rot)+uz**2 * (1-np.cos(theta_rot))]])

        # Convert to new coordinate system
        rotatedPosition = np.matmul(R_theta, relayPosition)
        #print(rotatedPosition)
        roverSystemPosition = (rotatedPosition-np.matmul(R_theta, roverPosition))
        #print(roverSystemPosition)
        r_relayRovSys = np.sqrt(roverSystemPosition[0,0] ** 2 + roverSystemPosition[0,1] ** 2 + roverSystemPosition[0,2] ** 2)
        #print(r_relayRovSys)
        lambda_relayRovSys = np.arctan2(roverSystemPosition[0,1], roverSystemPosition[0,0])
        phi_relayRovSys = np.arctan2(roverSystemPosition[0,2], np.sqrt(roverSystemPosition[0,0] ** 2 + roverSystemPosition[0,1] ** 2))

        if phi_relayRovSys >= np.radians(elevationAngleDeg):
            Comm = True
        return Comm, r_relayRovSys

    def plotOrbit(self, colorOrbit="blue"):  # To be used only for illustration only in the report
        '''
        Plot the complete orbit of the spacecraft
        @param: colorOrbit, the colour to use for the plot
        '''

        T = relay.orbitParameters(self)[0]
        time, dt = 0, 100
        xPlot, yPlot, zPlot = [], [], []
        while time < T * 1.01:
            theta = relay.getPosition(self, t=time)
            cartCoordsX = relay.KeplerToCartesian(self, true_anomaly=theta)[0]; xPlot.append(cartCoordsX)
            cartCoordsY = relay.KeplerToCartesian(self, true_anomaly=theta)[1]; yPlot.append(cartCoordsY)
            cartCoordsZ = relay.KeplerToCartesian(self, true_anomaly=theta)[2]; zPlot.append(cartCoordsZ)
            time += dt
        ax.plot(xPlot, yPlot, zPlot, color=colorOrbit, label=self.name)

    def plotSC(self, time, colorSC):  # To be used only for illustration only in the report
        '''
        Plot a sphere at the position of the spacecraft at the indicated time
        @param: time, the time variable defining the position of the spacecraft in the orbit
        @param: colorSC, the color of the spacecraft for plotting
        '''
        theta = relay.getPosition(self, t=time)
        X, Y, Z = relay.KeplerToCartesian(self, true_anomaly=theta)
        r = 500
        phi = np.linspace(-np.pi / 2, np.pi / 2, 100)
        lamb = np.linspace(-np.pi, np.pi, 100)
        x = r * np.outer(np.cos(lamb), np.cos(phi))
        y = r * np.outer(np.sin(lamb), np.cos(phi))
        z = r * np.outer(np.ones(np.size(lamb)), np.sin(phi))

        # Plot
        ax.plot_surface(X-x, Y-y, Z-z, color=colorSC)

    def performance(self, contactTime = 0):
        self.commTime += contactTime
        return True

    def results(self): # ATTENTION 0 !!!!!!!!
        print(f"{self.name} availability: {100* self.commTime/(T_end)} %")
        return self.inCommLatRec, self.inCommLongRec, self.noCommLatRec, self.noCommLongRec, self.commDuration, 100* self.commTime/(T_end)


# here the rover is the GS
KIR = GS(Name="KIR", planetRadius=R_Earth, latitudeDeg=67.85713, longitudeDeg=20.96432)
RED = GS(Name="RED", planetRadius=R_Earth, latitudeDeg=50.00046, longitudeDeg=5.145344)
KRU = GS(Name="KRU", planetRadius=R_Earth, latitudeDeg=5.251606, longitudeDeg=-52.80466)
SMA = GS(Name="SMA", planetRadius=R_Earth, latitudeDeg=36.99725, longitudeDeg=-25.13572) # Check if this one is necessary
AGO = GS(Name="AGO", planetRadius=R_Earth, latitudeDeg=-33.133333, longitudeDeg=-70.666667) # Check if this one is necessary
MAL = GS(Name="MAL", planetRadius=R_Earth, latitudeDeg=-2.995556, longitudeDeg=40.194511) # Check if this one is necessary
SG =  GS(Name="SG",  planetRadius=R_Earth, latitudeDeg=78.229772, longitudeDeg=15.407786) # Check if this one is necessary

ground = [KIR, RED, KRU, SMA, AGO, MAL, SG]
#180*satID/len(satNums)
satNums = np.arange(1, 30, 1)
satsA = [relay(semiMajorAxis=350+6371, eccentricity=0, inclinationRad=np.deg2rad(180*satID/len(satNums)), argumentOfPericenterRad=0,
               ascendingNodeRad=0, timeOfPericenterPassage=0, sc_name="SCA_" + str(satID), mu=mu_Earth) for satID in satNums]

satsB = [relay(semiMajorAxis=350+6371, eccentricity=0, inclinationRad=np.deg2rad(180*satID/len(satNums)), argumentOfPericenterRad=0,
               ascendingNodeRad=60, timeOfPericenterPassage=0, sc_name="SCB_" + str(satID), mu=mu_Earth) for satID in satNums]

satsC = [relay(semiMajorAxis=350+6371, eccentricity=0, inclinationRad=np.deg2rad(180*satID/len(satNums)), argumentOfPericenterRad=0,
               ascendingNodeRad=120, timeOfPericenterPassage=0, sc_name="SCC_" + str(satID), mu=mu_Earth) for satID in satNums]

satList = satsA + satsB + satsC
for sc in satList:
    sc.plotOrbit()

t, dt = 0, 10
while t < T_end:
    if t % 100:
        print(f"Percentage Done = {np.round(100 * t / T_end, 3)}%")

    for sat in satList:
        j = False
        theta = sat.getPosition(t)  # Get Keplerian position of the satellite
        xSAT, ySAT, zSAT = sat.KeplerToCartesian(theta)  # Convert to cartesian coordinates
        for gs in ground:
            xGS, yGS, zGS = gs.cartesianCoords()  # Get the current coordinates of all ground stations
            Comm_bool, distance_SC_rover = sat.inCommunicationCone([xSAT, ySAT, zSAT], [xGS, yGS, zGS],
                                                                        elevationAngleDeg=elevDeg)
            if Comm_bool:  # Determine if relay communication is possible
                j = True
                sat.performance(contactTime=dt)
        sat.groundTrackData(xSAT, ySAT, zSAT, t, Earth_Rotation_Rate, j)
        sat.J2(dt)
    for gs in ground:
        gs.updatePosition(Earth_Rotation_Rate, dt)
    t += dt

colors = ["red", "blue"]
bg = plt.imread("land_ocean_ice_350.jpg")
color_map = matplotlib.colors.ListedColormap(colors)
ext = [-180, 180, -90, 90]
fig, ax = plt.subplots()

for sat in satList:
    ax.imshow(bg, extent=ext)
    inlatrec, inlongrec, nolatrec, nolongrec, commDuration, perAvailability = sat.results()
    commDuration = np.array(commDuration)
    ax.scatter(nolongrec, nolatrec, color='blue', s=1)
    ax.scatter(inlongrec, inlatrec, color='red', s=1)

    # Contact window duration
    indices = (np.arange(1, len(commDuration)+1, 1))
    commIntervals = indices[commDuration == 1]
    commIntervals = np.append(commIntervals, np.array([10**4]))
    i, index0 = 0, 0
    while i <= len(commIntervals)-2:
        if commIntervals[i]+1 == commIntervals[i+1]:
            pass
        else:
            window = (commIntervals[i] - commIntervals[index0])*dt/60  # minutes
            # print(window, "min")
            index0 = i + 1
        i += 1
GS_lat = [67.85713, 50.00046, 5.251606, 36.99725, -33.133333, -2.995556, 78.229772]
GS_long = [20.96432, 5.145344, -52.80466, -25.13572, -70.666667, 40.194511, 15.407786]

AGO = GS(Name="AGO", planetRadius=R_Earth, latitudeDeg=-33.133333, longitudeDeg=-70.666667) # Check if this one is necessary
MAL = GS(Name="MAL", planetRadius=R_Earth, latitudeDeg=-2.995556, longitudeDeg=40.194511) # Check if this one is necessary
SG = GS(Name="SG",  planetRadius=R_Earth, latitudeDeg=78.229772, longitudeDeg=15.407786) # Check if this one is necessary
ax.scatter(GS_long, GS_lat, color='green', s=10)
ax.set_xlabel(r'Longitude $\lambda$ [deg]')
ax.set_ylabel(r'Latitude $\phi$ [deg]')
plt.show()
