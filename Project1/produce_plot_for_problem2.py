import numpy as np
import matplotlib.pyplot as plt

# Read file and extraxt values for x and u
filename = "problem2_output.txt"
with open(filename) as f:
    lines = f.readlines()

x_values = []
u_values = []
for line in lines:
    x, u = [float(scientific_notation) for scientific_notation in line.strip("\n").split(" ")]
    x_values.append(x)
    u_values.append(u)
    
# Plot
fig, ax = plt.subplots()
ax.plot(x_values, u_values)
ax.set_title('Plot of the exact solution for u(x)')
plt.xlabel('x')
plt.ylabel('u(x)', rotation=0)
fig.savefig("problem2_plot.pdf")