from matplotlib import pyplot as plt
import numpy as np

threshold = 0.01/10
masses = np.arange(1, 3100, 1)
number = 0.1*masses**0.75*threshold**(-1.71)

plt.figure()
plt.plot(masses/2, number/2)
plt.xlabel('Mass of each of the counterparts [kg] (equivalent to the average of the two)')
plt.ylabel('Number of debris particles larger than 1 cm in each orbit')
plt.xlim([0, 1500])
plt.show()
