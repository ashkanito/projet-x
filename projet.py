import numpy as np
import csv
import matplotlib.pyplot as plt

# Reading data from the C file in Python
file_path = 'output.txt'

# Read the data from the file
with open(file_path, 'r') as file:
    data = file.readlines()

altitude = np.genfromtxt('output.txt', delimiter = ';', dtype='float', usecols=0)
mass = np.genfromtxt('output.txt', delimiter = ';', dtype='float', usecols=1)
speed = np.genfromtxt('output.txt', delimiter = ';', dtype='float', usecols=2)
time = np.genfromtxt('output.txt', delimiter = ';', dtype='float', usecols=3)


# Creating Figure 1
plt.figure(1)
plt.plot(time, altitude)
plt.title('Figure 1 : altitude in function of time')
plt.xlabel('time')
plt.ylabel('altitude')
plt.legend()

# Creating Figure 2
plt.figure(2)
plt.plot(time, mass)
plt.title('Figure 2 : mass in function of time')
plt.xlabel('time')
plt.ylabel('mass')
plt.legend()

# Creating Figure 3
plt.figure(3)
plt.plot(time, speed)
plt.title('Figure 3 : speed in function of time')
plt.xlabel('time')
plt.ylabel('speed')
plt.legend()



plt.show()