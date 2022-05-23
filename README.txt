In order to run the code, first use the Makefile to compile.

The program is then run by using the command:

mpirun --bind-to none -np num_processes ./malaria_sim num_simulations create_output

num_simulations is the amount of simulations performed. create_output needs to be either 0 or 1.

If create_output is 0, no output files will be created. If create_output is 1, three output files will
be created. susceptible.txt the amount of simulations that ended up in each bin.

interval_boundaries.txt contains the boundaries of the bin intervals.

processor_timings.txt contains the average time each processor took to reach t = 25, 50, 75 and 100.

The provided shell script run.sh can be used to both run the simulations as well
as plot a histogram of susceptible humans at the end of the simulations. The range of
the bins is also printed in the standard output. Use run.sh as follows:

./run.sh num_simulations num_processes