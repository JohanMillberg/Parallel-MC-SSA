#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "prop.c"
#include <sys/time.h>

void runGillespie(int xInitial[], double tInitial, double T, int P[], int* result, double checkpointTimings[]);
double gillespieIteration(int* result, int T, int P[], int rowsInP, int colsInP);
double generateRandom();
int writeOutput(char *fileName, int *output, int outputSize);
int writeProcessorTimings(double *output, int size);
double calculateMeanOfInterval(double* array, int simulationAmount, int timeInterval);

int main(int argc, char* argv[]) {

    char* outputName = argv[1];    
    int N = atoi(argv[2]);

    struct timeval time;
    gettimeofday(&time,NULL);
    // Seeding using time
    srand((time.tv_sec*1000)+(time.tv_usec / 1000));
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

    // Used for random seeding
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
    MPI_Datatype susceptibleHumansVector, tempType1;
    MPI_Type_vector(simulationsPerProcess, 1, 7, MPI_INT, &tempType1);
    MPI_Type_create_resized(tempType1, 0, sizeof(int), &susceptibleHumansVector);
    MPI_Type_commit(&susceptibleHumansVector);

    // Initialize the local result matrix
    int* X = malloc(sizeof(int) * 7 * simulationsPerProcess);

    // Global result matrix
    int* totalX;

    // Global vector for values of X(1,1:N)
    int* totalX1;

    double* simulationTimings = malloc(sizeof(double)*simulationsPerProcess*4);
    double* averageTimings = malloc(sizeof(double)*4);

    // Start the timer
    double timeStart = MPI_Wtime();

    // Run the Gillespie simulation N / size times
    for (int i = 0; i < simulationsPerProcess; i++) {
        double checkPointTimings[4];
        runGillespie(x0, 0, T, P, &(X[i*7]), checkPointTimings);

        // Store the times to reach each checkpoint (t = 25, 50, 75, 100)
        memcpy(&(simulationTimings[i*4]), checkPointTimings, sizeof(double)*4);
    }

    // Stop the timer
    double maxTime;
    double executionTime = MPI_Wtime() - timeStart;

    // Calculate the average time to reach each time
    for (int i = 0; i < 4; i++) {
        averageTimings[i] = calculateMeanOfInterval(simulationTimings, simulationsPerProcess, i);
    }

    double* allAverageTimings;
    MPI_Win win;

    // Declare accessible memory area
    MPI_Alloc_mem(sizeof(double)*size*4, MPI_INFO_NULL, &allAverageTimings);
    MPI_Win_create(allAverageTimings, sizeof(double)*4*size, sizeof(double)*4, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);

    // Write the average timings of each process to the window
    MPI_Put(averageTimings, 4, MPI_DOUBLE, 0, worldRank, 4, MPI_DOUBLE, win);

    if (worldRank == 0) {
        totalX1 = malloc(sizeof(int) * N);
        totalX = malloc(sizeof(int)*N*7);
    }

    //MPI_Gather(X, simulationsPerProcess*7, MPI_INT, totalX, simulationsPerProcess*7, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(X, 1, susceptibleHumansVector, totalX1, simulationsPerProcess, MPI_INT, 0, MPI_COMM_WORLD);

    // Find the largest execution time
    MPI_Reduce(&executionTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Win_fence(0, win);

    // Print final time and write to output files
    if (worldRank == 0) {
        printf("Final time: %lf\n", maxTime);
        writeOutput(outputName, totalX1, N);
        writeProcessorTimings(allAverageTimings, size);

        free(totalX1);
        free(totalX);
    }
    MPI_Win_free(&win);
    MPI_Free_mem(allAverageTimings);

    free(X);
    free(simulationTimings);
    free(averageTimings);
    MPI_Type_free(&tempType1);
    MPI_Type_free(&susceptibleHumansVector);

    MPI_Finalize();
}

// Function used to calculate the mean execution time for each time interval
double calculateMeanOfInterval(double* array, int simulationAmount, int timeInterval) {
    double mean = 0;
    for (int i = 0; i < simulationAmount; i++) {
        mean += array[i*4 + timeInterval];
    }
    return (mean / simulationAmount);
}

// Run the Gillespie direct method
void runGillespie(int x[], double tInitial, double T, int P[], int* result, double checkpointTimings[]) {
    double t = tInitial;
    int rowsInP = 15;
    int colsInP = 7;

    for (int i = 0; i < colsInP; i++) {
        result[i] = x[i];
    }

    double tau;
    int r;
    double timerStart = MPI_Wtime();

    // Define the timing intervals
    double timeIntervals[4] = {T/4, T/2, 3*T/4, T};

    while (t < timeIntervals[0]) {
        tau = gillespieIteration(result, T, P, rowsInP, colsInP);
        t += tau;
    }
    checkpointTimings[0] = MPI_Wtime() - timerStart;
    while (t < timeIntervals[1]) {
        tau = gillespieIteration(result, T, P, rowsInP, colsInP);
        t += tau;
    }
    checkpointTimings[1] = MPI_Wtime() - timerStart;
    while (t < timeIntervals[2]) {
        tau = gillespieIteration(result, T, P, rowsInP, colsInP);
        t += tau;
    }
    checkpointTimings[2] = MPI_Wtime() - timerStart;
    while (t < timeIntervals[3]) {
        tau = gillespieIteration(result, T, P, rowsInP, colsInP);
        t += tau;
    }
    checkpointTimings[3] = MPI_Wtime() - timerStart;

}

// Helper function executing one Gillespie iteration
double gillespieIteration(int* result, int T, int P[], int rowsInP, int colsInP) {
    double w[rowsInP];
    double u1, u2, tau;
    int r;

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

    // Store the result
    for (int i = 0; i < colsInP; i++) {
        result[i] += P[r*colsInP + i];
    }
    return tau;
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

int writeProcessorTimings(double *output, int size) {
    FILE *file;
    if (NULL == (file = fopen("processor_timings.txt", "w"))) {
        perror("Couldn't open output file");
        return -1;
    }
    fprintf(file, "Average time it takes for each processor to reach a certain time interval:\n\n");
    for (int i = 0; i < 4; i++) {
        fprintf(file, "Average time(s) for time interval %d: ", i);
        for (int j = 0; j < size; j++) {
            fprintf(file, "%lf ", output[j*4 + i]);
        }
        if (i != 3)
            fprintf(file, "\n");
    }

    return 0;
}