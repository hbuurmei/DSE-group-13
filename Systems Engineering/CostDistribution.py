""""This program visualizes the uncertainties in the mission cost estimation"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

plt.rcParams['text.usetex'] = True

# mu = 543771.6094423*1e-3
# sigma = 0.268355861325*mu
# print(sigma)
# x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
# plt.plot(x, stats.norm.pdf(x, mu, sigma))
# plt.xlabel(r"Cost estimation \euro")
# plt.ylabel("Probability density")
# plt.axvline(x=mu, color='k', linestyle='--')
# plt.show()




# import libraries
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

# mean
mu = 543771.6094423*1e-3
# standard deviation
sigma = 0.268355861325*mu
# number of observations in our dataset
n = 1_000_000
# construct a random Generator
s = np.random.default_rng(0)
# generate an array of 1_000_000 random values that follow a normal distribution with mean mu and st.dev. sigma
r = s.normal(mu, sigma, n)

# transform the array into a dataframe
df = pd.DataFrame(r, columns=["value"])
# view descriptive statistics on the data
dstats = df.describe()
# save stats above into individual variables
mean = dstats.loc["mean"].value
std = dstats.loc["std"].value
p25 = dstats.loc["25%"].value
median = dstats.loc["50%"].value
p75 = dstats.loc["75%"].value
# minv = dstats.loc["min"].value
# maxv = dstats.loc["max"].value

# plot distribution
plt.figure(figsize=(10, 5))

plt.hist(df["value"], bins=100, color="teal")

# plt.axvline(x=mean, color="red", ls="-")
plt.axvline(x=mean+std, color="orange", ls="-")
plt.axvline(x=mean-std, color="salmon", ls="-")
plt.axvline(x=p25, color="purple", ls=":")
plt.axvline(x=median, color="black", ls=":")
plt.axvline(x=p75, color="lime", ls=":")
# plt.axvline(x=minv, color="blue", ls="-.")
# plt.axvline(x=maxv, color="navy", ls="-.")

# plt.title("Normal Distribution", size=16, color="blue", pad=20)
plt.xlabel("Mission cost [$\euro$]") #, color="blue")
plt.ylabel("Frequency Density") #, color="blue")

plt.legend(["+Standard Deviation",
"-Standard Deviation", "25th Percentile",
"50th Percentile or Median", "75th Percentile"])
# ,"Minimum Value", "Maximum Value"])

plt.savefig("blog_fig1.jpg", dpi=200)

plt.tight_layout()
plt.show()
