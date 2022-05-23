import matplotlib.pyplot as plt

n_cores = [1, 2, 4, 8, 16]
speedup = [1, 1.978, 3.954, 7.565, 14.206]
speedup_ideal = [1, 2, 4, 8, 16]

plt.plot(n_cores, speedup, label='Actual speedup')
plt.plot(n_cores, speedup_ideal, label='Ideal speedup')
plt.xlabel('Number of processes')
plt.ylabel('Speedup')
plt.title('Speedup as number of processes increase')
plt.legend()
plt.show()


time1 = [142.096, 142.934, 146.188, 148.743, 160.115]

plt.plot(n_cores, time1)
plt.xlabel('Number of processes')
plt.ylabel('Time (s)')
plt.title('Execution time as #simulations and #processes increase')
plt.legend()
plt.show()

"""
362.159494
81.985210
50.131371
37.377782
7488

1 40.391694
4 8.888543
9 5.504375
16 3.887872
25 1.717179
3600
"""