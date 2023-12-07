/*
 * Trabalho PCD - Rainbow Game of Life (Versão MPI)
 * Nome: Rodrigo Peixe Oliveira
 * RA: 147873
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "mpi.h"


char version[] = "mpi";
const int N = 2048;
const int maxGenerations = 2000;
double **grid, **newGrid;
int rank, start, end, numProcs;


void printGrid(double** grid, int generation, int n);
void setInitialGeneration(double** grid);
void synchronizeGrids(double** grid);
double getNewValue(double** grid, int i, int j);
int countAlive(double** grid);
int getNeighbors(double** grid, int i, int j);
double getNeighborsMean(double** grid, int i, int j);
int enforceBorders(int a);
double **createSquareMatrix(int size);
void swap(double*** a, double*** b);


int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    struct timeval timeStart, timeEnd;
    int tmili;
    gettimeofday(&timeStart, NULL);

    start = ceil(rank * N / numProcs);
    end = ceil((rank+1) * N / numProcs);
    
    int currentGeneration, i, j;
    

    grid = createSquareMatrix(N);
    newGrid = createSquareMatrix(N);

    setInitialGeneration(grid);
    

    for (currentGeneration = 1; currentGeneration <= maxGenerations; currentGeneration++) {
        if (numProcs > 1) {
            synchronizeGrids(grid);
        }

        for(i = start; i < end; i++) {
            for(j = 0; j < N; j++) {
                newGrid[i][j] = getNewValue(grid, i, j);
            }
        }

        swap(&grid, &newGrid);
    }


    int localAlive = countAlive(grid);
    int totalAlive;
    MPI_Reduce(&localAlive, &totalAlive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    
    gettimeofday(&timeEnd, NULL);
    tmili = (int) (1000 * (timeEnd.tv_sec - timeStart.tv_sec) + (timeEnd.tv_usec - timeStart.tv_usec) / 1000);

    if (rank == 0) {
        printf("Processes: %d\n", numProcs);
        printf("Generation %d: %d alive\n", currentGeneration-1, totalAlive);
        printf("Total time: %d ms\n", tmili);
        printf("----------\n");
    }
    
    MPI_Finalize();
    
    return 0;
}


void setInitialGeneration(double** grid) {
    //GLIDER
    int lin = 1, col = 1;
    if (lin >= start && lin < end) {
        grid[lin  ][col+1] = 1.0;
    }
    if (lin+1 >= start && lin < end) {
        grid[lin+1][col+2] = 1.0;
    }
    if (lin+2 >= start && lin < end) {
        grid[lin+2][col  ] = 1.0;
        grid[lin+2][col+1] = 1.0;
        grid[lin+2][col+2] = 1.0;
    }
    
    //R-pentomino
    lin =10; col = 30;
    if (lin >= start && lin < end) {
        grid[lin  ][col+1] = 1.0;
        grid[lin  ][col+2] = 1.0;
    }
    if (lin+1 >= start && lin < end) {
        grid[lin+1][col  ] = 1.0;
        grid[lin+1][col+1] = 1.0;
    }
    if (lin+2 >= start && lin < end) {
        grid[lin+2][col+1] = 1.0;
    }
}

void synchronizeGrids(double** grid) {
    MPI_Request request1;
    MPI_Request request2;
    MPI_Request request3;
    MPI_Request request4;
    MPI_Status status;

    int previous = rank == 0 ? numProcs-1 : rank-1;
    int next = rank == numProcs-1 ? 0 : rank+1;

    int previousRow = start == 0 ? N-1 : start-1;
    int nextRow = end == N? 0 : end;

    MPI_Isend(grid[start], N, MPI_DOUBLE, previous, 10, MPI_COMM_WORLD, &request1);
    MPI_Irecv(grid[previousRow], N, MPI_DOUBLE, previous, 10, MPI_COMM_WORLD, &request2);
    MPI_Isend(grid[end-1], N, MPI_DOUBLE, next, 10, MPI_COMM_WORLD, &request3);
    MPI_Irecv(grid[nextRow], N, MPI_DOUBLE, next, 10, MPI_COMM_WORLD, &request4);

    MPI_Wait(&request1, &status);
    MPI_Wait(&request2, &status);
    MPI_Wait(&request3, &status);
    MPI_Wait(&request4, &status);
}


double getNewValue(double** grid, int i, int j) {
    double currentValue = grid[i][j];
    int cellsAlive = getNeighbors(grid, i, j);
    
    double newValue;
    
    if (currentValue > 0 && (cellsAlive == 2 || cellsAlive == 3)) {
        newValue = currentValue;
    }
    else if (cellsAlive == 3) {
        newValue = getNeighborsMean(grid, i, j);
    }
    else {
        newValue = 0.0;
    }

    return newValue;
}


int countAlive(double** grid) {
    int i, j, total;

    total = 0;

    for (i = start; i < end; i++) {
        for (j = 0; j < N; j++) {
            if (grid[i][j] > 0) {
                total++;
            }
        }
    }

    return total;
}


int getNeighbors(double** grid, int i, int j) {
    int k, l, total;
    
    total = 0;
    
    for (k = i-1; k <= i+1; k++) {
        for (l = j-1; l <= j+1; l++) {
            if (k == i && l == j) {
                continue;
            }
            
            int a, b;
            a = enforceBorders(k);
            b = enforceBorders(l);

            if (grid[a][b] > 0) {
                total++;
            }
        }
    }

    return total;
}


double getNeighborsMean(double** grid, int i, int j) {
    int k, l;
    double total, mean;
    
    total = 0.0;
    
    for (k = i-1; k <= i+1; k++) {
        for (l = j-1; l <= j+1; l++) {
            if (k == i && l == j) {
                continue;
            }
            
            int a, b;
            a = enforceBorders(k);
            b = enforceBorders(l);

            total += grid[a][b];
        }
    }

    mean = total / 8;

    return mean;
}


int enforceBorders(int a) {
    int b;
    if (a >= 0 && a <= N-1) {
        b = a;
    }
    else if (a == -1) {
        b = N-1;
    }
    else if (a == N) {
        b = 0;
    }
    else {
        printf("Unexpected value when enforcing borders: %d\n", a);
        exit(2);
    }
    return b;
}


double **createSquareMatrix(int size) {
    double **m;
    int i;
    
    m = calloc(size, sizeof(double *));
    for(i=0; i<size; i++) {
        m[i] = calloc(size, sizeof(double));
    }

    return m;
}


void swap(double*** a, double*** b) {
    double** aux = *a;
    *a = *b;
    *b = aux;
}