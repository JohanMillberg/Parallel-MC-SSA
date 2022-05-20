#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -A uppmax2022-2-11
#SBATCH -p core -n 16
#SBATCH -t 30:00

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1
make

echo "Strong scalability"
mpirun --bind-to none -np 1 ./malaria_sim output.txt 100000
mpirun --bind-to none -np 2 ./malaria_sim output.txt 100000
mpirun --bind-to none -np 4 ./malaria_sim output.txt 100000
mpirun --bind-to none -np 8 ./malaria_sim output.txt 100000
mpirun --bind-to none -np 16 ./malaria_sim output.txt 100000

echo "Weak scalability"
mpirun --bind-to none -np 1 ./malaria_sim output.txt 100000
mpirun --bind-to none -np 2 ./malaria_sim output.txt 200000
mpirun --bind-to none -np 4 ./malaria_sim output.txt 400000
mpirun --bind-to none -np 8 ./malaria_sim output.txt 800000
mpirun --bind-to none -np 16 ./malaria_sim output.txt 1600000