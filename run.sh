#Use the following command to run simulation: ./run.sh num_simulations num_processes

make

echo "Running $1 simulations on $2 processes..."

mpirun --bind-to none -np $2 ./malaria_sim $1 1

python plot_result.py

make clean