import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

## Problem 6
###################################################

# Extract the binary files obtained in Problem6.cpp
I_mat = pa.mat()
N_info = pa.mat()
I_mat.load("iterations_vector.bin")
N_info.load("N_information.bin")

# Convert to numpy
N_min, N_max = N_info[0,0], N_info[1,0]
N = np.arange(N_min, N_max+1)
I = np.zeros(N.size)
for idx in range(N.size):
    I[idx] = I_mat[idx, 0]

# Plot
fig, ax = plt.subplots()
ax.plot(N, I)
ax.set_title(f'Plot of the number of iterations before convergence\nNmin={int(N_min)}. Nmax={(int(N_max))}')
plt.xlabel('N')
plt.ylabel('I(N)', rotation=0)
fig.savefig("problem6a_plot.pdf")