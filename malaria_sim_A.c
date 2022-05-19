#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "prop.c"
#include <sys/time.h>

void runGillespie(int xInitial[], double tInitial, double T, int P[], int* result);
double generateRandom();
int writeOutput(char *fileName, int *output, int outputSize);

int main(int argc, char* argv[]) {
    struct timeval time;
    gettimeofday(&time,NULL);
    // Seeding using time
    srand((time.tv_sec*1000)+(time.tv_usec / 1000));
    char* outputName = argv[1];    
    int N = atoi(argv[2]);
    int worldRank;
    MPI_Status status;
    int size;

    double T = 100;
    int x0[] = {900, 900, 30, 330, 50, 270, 20};

    int P[] = {1, 0, 0, 0, 0, 0, 0,
                -1, 0, 0, 0, 0, 0, 0,
                -1, 0, 1, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0, 0,
                0, -1, 0, 0, 0, 0, 0,
                0, -1, 0, 1, 0, 0, 0,
                0, 0, -1, 0, 0, 0, 0,
                0, 0, -1, 0, 1, 0, 0,
                0, 0, 0, -1, 0, 0, 0,
                0, 0, 0, -1, 0, 1, 0,
                0, 0, 0, 0, -1, 0, 0,
                0, 0, 0, 0, -1, 0, 1,
                0, 0, 0, 0, 0, -1, 0,
                1, 0, 0, 0, 0, 0, -1,
                0, 0, 0, 0, 0, 0, -1
                };

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    int seedingNumbers[size];
    if (worldRank == 0) {
        for (int i = 0; i < size; i++) {
            seedingNumbers[i] = rand()%10000;
        }
    }

    MPI_Bcast(seedingNumbers, size, MPI_INT, 0, MPI_COMM_WORLD);

    srand(seedingNumbers[worldRank]);

    if (N%size != 0) {
        printf("Number of experiments must be divisible with number of processes.\n");
        MPI_Finalize();
        return 0;
    }

    int simulationsPerProcess = N / size;

    // Create the datatype to store the values of X(1,1:N)
    MPI_Datatype susceptibleHumansVector, tempType;
    MPI_Type_vector(simulationsPerProcess, 1, 7, MPI_INT, &tempType);
    MPI_Type_create_resized(tempType, 0, sizeof(int), &susceptibleHumansVector);
    MPI_Type_commit(&susceptibleHumansVector);

    int* X = malloc(sizeof(int) * 7 * simulationsPerProcess);
    int* totalX;
    int* totalX1;

    // Start the timer
    double timeStart = MPI_Wtime();

    for (int i = 0; i < simulationsPerProcess; i++) {
        runGillespie(x0, 0, T, P, &(X[i*7]));
    }

    double maxTime;
    double executionTime = MPI_Wtime() - timeStart;

    if (worldRank == 0) {
        totalX1 = malloc(sizeof(int) * N);
        totalX = malloc(sizeof(int)*N*7);
    }
    MPI_Gather(X, simulationsPerProcess*7, MPI_INT, totalX, simulationsPerProcess*7, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(X, 1, susceptibleHumansVector, totalX1, simulationsPerProcess, MPI_INT, 0, MPI_COMM_WORLD);

    // Find the largest execution time
    MPI_Reduce(&executionTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (worldRank == 0) {
        printf("Final time: %lf\n", maxTime);
        writeOutput(outputName, totalX1, N);
        free(totalX1);
        free(totalX);
    }

    free(X);
    MPI_Type_free(&tempType);
    MPI_Type_free(&susceptibleHumansVector);

    MPI_Finalize();
}

// Run the Gillespie direct method
void runGillespie(int x[], double tInitial, double T, int P[], int* result) {
    double t = tInitial;
    int rowsInP = 15;
    int colsInP = 7;

    for (int i = 0; i < colsInP; i++) {
        result[i] = x[i];
    }

    double w[rowsInP];
    double u1, u2, tau;
    int r;

    while (t < T) {
        double a0 = 0;
        prop(result, w);
        for (int i = 0; i < rowsInP; i++) {
            a0 += w[i];
        }
    
        u1 = generateRandom();
        u2 = generateRandom();
        tau = -log(u1)/a0;

        double temp = a0*u2;
        r = 0;
        double sum = 0;
        for (int i = 0; i < rowsInP; i++) {
            sum += w[i];
            if (sum >= temp) {
                break;
            }
            r++;
        }

        for (int i = 0; i < colsInP; i++) {
            result[i] += P[r*colsInP + i];
        }
        t += tau;
    }

}

double generateRandom() {
    double randomNumber = (double)rand()/(double)(RAND_MAX);
    return randomNumber;
}

// Function for handling output. Inspired by the readInput function given in A1.
int writeOutput(char *fileName, int *output, int outputSize) {
    FILE *file;
    if (NULL == (file = fopen(fileName, "w"))) {
        perror("Couldn't open output file");
        return -1;
    }
    for (int j = 0; j < outputSize-1; j++) {
        if (0 > fprintf(file, "%d ", output[j])) {
            perror("Couldn't write to output file");
        }
    }

    if (0 > fprintf(file, "%d", output[outputSize-1])) {
        perror("Couldn't write to output file");
    }

    if (0 > fprintf(file, "\n")) {
        perror("Couldn't write to output file");
    }
    if (0 != fclose(file)) {
        perror("Warning: couldn't close output file");
    }
    return 0;
}