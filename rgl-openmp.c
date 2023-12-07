/*
 * Trabalho PCD - Rainbow Game of Life (Versão OpenMP)
 * Nome: Rodrigo Peixe Oliveira
 * RA: 147873
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>


char version[] = "openmp";
const int N = 2048;
const int maxGenerations = 2000;
double **grid, **newGrid;


void setInitialGeneration(double** grid);
double getNewValue(double** grid, int i, int j);
int countAlive(double** grid, int* time_elapsed);
int getNeighbors(double** grid, int i, int j);
double getNeighborsMean(double** grid, int i, int j);
int enforceBorders(int a);
double **createSquareMatrix(int size);
void swap(double*** a, double*** b);


int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Usage: %s <numThreads>\n", argv[0]);
        exit(1);
    }
    const int numThreads = atoi(argv[1]);
    omp_set_num_threads(numThreads);

    struct timeval timeStart, timeEnd;
    int tmili;
    gettimeofday(&timeStart, NULL);

    
    int currentGeneration, i, j;

    grid = createSquareMatrix(N);
    newGrid = createSquareMatrix(N);

    setInitialGeneration(grid);
    
    int threadID;
    #pragma omp parallel private(threadID, currentGeneration, i, j) shared(maxGenerations, grid, newGrid, N) default(none)
    {
        threadID = omp_get_thread_num();

        for (currentGeneration = 1; currentGeneration <= maxGenerations; currentGeneration++) {
            for(i = 0; i < N; i++) {
                #pragma omp for
                for(j = 0; j < N; j++) {
                    newGrid[i][j] = getNewValue(grid, i, j);
                }
            }
            #pragma omp barrier
            if (threadID == 0) {
                struct timeval timeStartNP, timeEndNP;
                gettimeofday(&timeStartNP, NULL);

                swap(&grid, &newGrid);

                gettimeofday(&timeEndNP, NULL);
            }
            #pragma omp barrier
        }
    }

    int tuCountAlive;
    int alive = countAlive(grid, &tuCountAlive);
    
    gettimeofday(&timeEnd, NULL);
    tmili = (int) (1000 * (timeEnd.tv_sec - timeStart.tv_sec) + (timeEnd.tv_usec - timeStart.tv_usec) / 1000);


    printf("Threads: %d\n", numThreads);
    printf("Generation 2000: %d alive\n", alive);
    printf("Total time: %d ms\n", tmili);
    printf("countAlive() time: %d us\n", tuCountAlive);
    printf("----------\n");
    
    return 0;
}


void setInitialGeneration(double** grid) {
    //GLIDER
    int lin = 1, col = 1;
    grid[lin  ][col+1] = 1.0;
    grid[lin+1][col+2] = 1.0;
    grid[lin+2][col  ] = 1.0;
    grid[lin+2][col+1] = 1.0;
    grid[lin+2][col+2] = 1.0;
    
    //R-pentomino
    lin =10; col = 30;
    grid[lin  ][col+1] = 1.0;
    grid[lin  ][col+2] = 1.0;
    grid[lin+1][col  ] = 1.0;
    grid[lin+1][col+1] = 1.0;
    grid[lin+2][col+1] = 1.0;
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


int countAlive(double** grid, int* time_elapsed) {
    int i, j, total;
    struct timeval timeStart, timeEnd;

    total = 0;
    gettimeofday(&timeStart, NULL);

    #pragma omp parallel for private(i, j) reduction(+:total)
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (grid[i][j] > 0) {
                total++;
            }
        }
    }

    gettimeofday(&timeEnd, NULL);

    *time_elapsed = (int) (timeEnd.tv_usec - timeStart.tv_usec);

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