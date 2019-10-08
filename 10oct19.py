# 10/10/19

import matplotlib.pyplot as plt

import readdata

# import arrays from data file
exoplanet_array = readdata.exoplanet_array
tele_array = readdata.tele_array

# plot
rank = []
metric = []
radius = []
mass = []
semi_major_axis = []

for i in range(len(exoplanet_array)):
    rank.append(i)
    metric.append(exoplanet_array[i].decision_metric)
    if exoplanet_array[i].radius == '': radius.append(0)
    else: radius.append(exoplanet_array[i].radius)
    #radius.append(exoplanet_array[i].radius)
    if exoplanet_array[i].mass == '': mass.append(0)
    else: mass.append(exoplanet_array[i].mass)
    semi_major_axis.append(exoplanet_array[i].semi_major_axis)

x = rank
y = metric

plt.scatter(x, y)
plt.xlabel(x)
plt.ylabel(y)
plt.show()