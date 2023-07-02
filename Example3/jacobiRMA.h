#ifndef JACOBIRMA_H
#define JACOBIRMA_H

#include "poisson.h"

void sweep2d(double a[][maxn], double f[][maxn], int nx,
	int* xs, int* xe, int* ys, int* ye, double b[][maxn]);
void exchangrma2d(double x[][maxn], int nx, int* xs, int* xe, int* ys, int* ye, int rank, MPI_Comm comm,
    int nbrleft, int nbrright, int nbrup, int nbrdown, MPI_Win win);
void exchangrma_2d_pscw(double x[][maxn], int ny, int* xs, int* xe, int* ys, int* ye, MPI_Win win,
	int nbrleft, int nbrright, int nbrup, int nbrdown, MPI_Group origin[2], MPI_Group target[2]);
double griddiff2d(double a[][maxn], double b[][maxn], int* xs, int* xe, int* ys, int* ye);

#endif
