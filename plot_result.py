from re import I
import matplotlib.pyplot as plt
import numpy as np

hist = np.loadtxt("susceptible.txt", dtype=int, delimiter=' ')
bin_edges = np.loadtxt("interval_boundaries.txt", dtype=int, delimiter=' ')

n_bins = 20
min_value = min(list(bin_edges))
max_value = max(list(bin_edges))
width = (max_value - min_value) / (n_bins+1)
range_of_bins = (max_value-min_value)/n_bins

print("Histogram bin intervals:")
for i in range(len(bin_edges)-1):
    print(f"[{bin_edges[i]}, {bin_edges[i+1]}]")

plt.bar(bin_edges[:-1], hist, width=width)
plt.xlim(min_value, max_value)
plt.xlabel('Susceptible humans')
plt.ylabel('Amount of simulations')
plt.title(f'Histogram showing the result of {sum(list(hist))} Monte Carlo computations')
plt.show()