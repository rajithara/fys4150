import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

## Problem 7
###################################################

# Extract the binary files obtained in Problem7.cpp
eigvec_9 = pa.mat()
eigvec_99 = pa.mat()
eigvec_9.load("eigvec_9.bin")
eigvec_99.load("eigvec_99.bin")

# Convert to numpy
eigvec_9_arr = np.zeros([9, 3])
eigvec_99_arr = np.zeros([99, 3])
for j in range(3):
    for i in range(9):
        eigvec_9_arr[i, j] = eigvec_9[i, j]
    for i in range(99):
        eigvec_99_arr[i, j] = eigvec_99[i, j]

# Plot
x_9 = np.arange(9)/9
x_99 = np.arange(99)/99
u_9_list = [eigvec_9_arr[:, col] for col in range(3)]
u_99_list = [eigvec_99_arr[:, col] for col in range(3)]

fig, ax = plt.subplots()
for i in range(3):
    ax.plot(x_9, u_9_list[i])
ax.set_title("n=10")
plt.xlabel(r'\hat{x}')
plt.ylabel(r'\hat{u}', rotation=0)
fig.savefig("problem7_n10_plot.pdf")

fig, ax = plt.subplots()
for i in range(3):
    ax.plot(x_99, u_99_list[i])
ax.set_title("n=10")
plt.xlabel(r'\hat{x}')
plt.ylabel(r'\hat{u}', rotation=0)
fig.savefig("problem7_n100_plot.pdf")