# Projeto 3

Projeto 3 da disciplina Programação Concorrente e Distribuída 2023-1 da UNIFESP

Aluno:
- Rodrigo Peixe Oliveira (RA: 147873)

Comandos utilizados:

```
gcc -o rgl-openmp rgl-openmp.c -fopenmp
./rgl-openmp <numThreads>


mpicc -o rgl-mpi rgl-mpi.c
mpiexec -np <numProcs> --oversubscribe ./rgl-mpi
```
