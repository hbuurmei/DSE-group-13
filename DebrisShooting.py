import numpy as np
import matplotlib.pyplot as plt

mu_Earth = 3.986004418e14 / 10**9  # km^3/s^2
R_Earth = 6371  # km
M_Earth = 5.972e24  # kg
Earth_day_duration = 24 * 60 *60 # s
# For plotting in 3D
fig = plt.figure(facecolor="white")
ax = plt.axes(projection='3d')

class planet:
    def __init__(self, radius, mass, day_duration, mu):
        self.mass = mass
        self.radius = radius
        self.T_revolution = day_duration
        self.mu = mu
        self.omega = 2 * np.pi/ day_duration

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
        ax.plot_surface(x, y, z, rstride=5, cstride=5, color="blue")  # plotting

class satellite:
    def __init__(self, semiMajorAxis, eccentricity, inclinationRad, argumentOfPericenterRad, ascendingNodeRad, timeOfPericenterPassage, sc_name, mu = mu_Earth):
        self.e = eccentricity
        self.a = semiMajorAxis
        self.t_p= timeOfPericenterPassage
        self.i = inclinationRad
        self.Omega = ascendingNodeRad
        self.w = argumentOfPericenterRad
        self.mu = mu
        self.name = sc_name

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

    def plotOrbit(self, colorOrbit="red"):  # To be used only for illustration only in the report
        '''
        Plot the complete orbit of the spacecraft
        @param: colorOrbit, the colour to use for the plot
        '''
        T = satellite.orbitParameters(self)[0]
        time, dt = 0, 1
        xPlot, yPlot, zPlot = [], [], []
        while time < T * 1.01:
            theta = satellite.getPosition(self, t=time)
            cartCoordsX = satellite.KeplerToCartesian(self, true_anomaly=theta)[0]; xPlot.append(cartCoordsX)
            cartCoordsY = satellite.KeplerToCartesian(self, true_anomaly=theta)[1]; yPlot.append(cartCoordsY)
            cartCoordsZ = satellite.KeplerToCartesian(self, true_anomaly=theta)[2]; zPlot.append(cartCoordsZ)
            time += dt
            print(cartCoordsX, cartCoordsY, cartCoordsZ)
        ax.plot(xPlot, yPlot, zPlot, color=colorOrbit, label=self.name)

    def plotSC(self, time, colorSC):  # To be used only for illustration only in the report
        '''
        Plot a sphere at the position of the spacecraft at the indicated time
        @param: time, the time variable defining the position of the spacecraft in the orbit
        @param: colorSC, the color of the spacecraft for plotting
        '''
        theta = satellite.getPosition(self, t=time)
        X, Y, Z = satellite.KeplerToCartesian(self, true_anomaly=theta)
        r = 500
        phi = np.linspace(-np.pi / 2, np.pi / 2, 100)
        lamb = np.linspace(-np.pi, np.pi, 100)
        x = r * np.outer(np.cos(lamb), np.cos(phi))
        y = r * np.outer(np.sin(lamb), np.cos(phi))
        z = r * np.outer(np.ones(np.size(lamb)), np.sin(phi))

        # Plot
        ax.plot_surface(X-x, Y-y, Z-z, color=colorSC)

Mars = planet(R_Earth, M_Earth, Earth_day_duration, mu_Earth)
debrisList = ...
# debris = [satellite(semiMajorAxis=, eccentricity=, inclinationRad=, argumentOfPericenterRad=, ascendingNodeRad=, timeOfPericenterPassage=, sc_name= "") for debrisID in debrisList]
