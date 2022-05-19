import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("susceptible.txt", dtype=int, delimiter=' ')

n_bins = 20
min_value = min(list(data))
max_value = max(list(data))
range_of_bins = (max_value-min_value)/n_bins

current_min = min_value
current_max = min_value + range_of_bins

print("Histogram bin intervals:")
for i in range(n_bins):
    print(f"[{round(current_min, 1)}, {round(current_max, 1)}]")
    current_min = current_min + range_of_bins
    current_max = current_max + range_of_bins

plt.hist(data, n_bins)
plt.xlabel('Susceptible humans')
plt.ylabel('Amount of simulations')
plt.title(f'Histogram showing the result of {len(list(data))} Monte Carlo computations')
plt.show()