import tradeoff_class as tc
import numpy as np
# Define paramters in tradeoff
# list of arguments:
# Name is the display name used in table printing
# weight is the assigned weight for the parameter
# func (default "LRTS") can either be LRTS, IRTS or DRTS see tutorial page in notion for description
# direc (default "HB") either "HB" for high is best or "LB" for Low is best
# Limitype (default "minmax") is how the range of values that will give map top [0,1]. there are three options "fixed"
# the rnage is manually defined,
# "minmax" the range is defined such that the worst option is 0 and the best is 1, "SD" the range is given as
# mean+-x*standard deviation
# Limit_val is only nesscarry for "fixed" where it defines the manual rangee, and "SD" where it defines x
# EXAMPLE: avi = tc.param(name="Availability", weight=1/3, Limitype="fixed", Limit_val=[1, 5])

risk = tc.param(name="Risk", weight=0.35, var=0.01, Limitype="minmax", direc="LB")
speed = tc.param(name="Debris removal speed", weight=0.324, var=0.01, Limitype="minmax", direc="HB")
power = tc.param(name="Power", weight=0.02, var=0.00002, Limitype="minmax", direc="LB")
complex = tc.param(name="Complexity", weight=0.0837, var=0.000016, Limitype="minmax", direc="LB")
volume = tc.param(name="Volume", weight=0.0439, var=0.000041, Limitype="minmax", direc="LB")
wet_mass = tc.param(name="Wet mass", weight=0.0975, var=0.000087, Limitype="minmax", direc="LB")
dry_mass = tc.param(name="Dry mass", weight=0.0810, var=0.000068, Limitype="minmax", direc="LB")
param_list = [complex, risk, speed, power, dry_mass, wet_mass, volume]

# Define the designs to be traded off source list being the values assigned
# EXAMPLE: gyro = tc.design(name="Gyroscope", sourcelist=[5, 3, 5])

scoring_matrix = np.array([[108.6, 45.34,	3, 49950, 3684.55, 3774.55, 96.6],
                  [108.6, 75.56, 5, 99900, 7369.1, 7579.1, 193.2],
                  [74.3,  65.52, 1, 10890, 868.56, 1475.5, 8.25],
                  [62.2, 80.86, 1, 8750, 4873, 12338.3, 13.55]])

# for i in range(0, 7):
    # column = scoring_matrix[:, i]
    # column = 4 * (column - np.amin(column)) / (np.amax(column) - np.amin(column)) + 1
    # scoring_matrix[:, i] = column

laser = tc.design(name="Laser", sourcelist=scoring_matrix[0])
laser_constellation = tc.design(name="Laser constellation", sourcelist=scoring_matrix[1])
foam = tc.design(name="Foam", sourcelist=scoring_matrix[2])
thruster = tc.design(name="Thruster", sourcelist=scoring_matrix[3])
design_list = [laser, laser_constellation, foam, thruster]


# colors for Latex table generation
colors = [tc.color("EF5350", "red"), tc.color("FB8C00", "orange"), tc.color("FFEB3B", "yellow"),
          tc.color("8BC34A", "green")]

# EXAMPLE
# tradeoff_att = tc.tradeoff(design_list=[gyro, sun, star, maf], param_list=[avi, acc, sam])

tradeoff_att = tc.tradeoff(design_list=design_list, param_list=param_list)
tradeoff_att.get_tradeoff()
tradeoff_att.get_output(language="latex", color_list=colors, width=10, rot="hor", caption="Tradeoff ASDR Concepts")

sensitivity_analysis = tc.sensitivity(tradeoff_att, samples=20000)
sensitivity_analysis.addto_weights()
sensitivity_analysis.get_sens()
print(sensitivity_analysis.per)
sensitivity_analysis.get_RMS()
print(sensitivity_analysis.RMS)





