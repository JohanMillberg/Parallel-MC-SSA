#Use the following command to run simulation: ./run.sh num_simulations num_processes

make

echo "Running $1 simulations on $2 processes..."

mpirun -np $2 ./malaria_sim susceptible.txt $1

python plot_result.py

make clean