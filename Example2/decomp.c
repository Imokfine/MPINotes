#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int MPE_Decomp1d(int n, int nprocs, int rank, int* s, int* e){
    int epp = n / nprocs; //elements per process without remainder
    int r = n % nprocs; //remainder

    *s = rank * epp + 1;
    *s = *s + ((rank < r) ? rank : r);
    if (rank < r) epp++;
    *e = *s + epp - 1;
    if (*e > n || rank == nprocs - 1) *e = n;

    return MPI_SUCCESS;
}

int MPE_Decomp2d(int nx, int ny, int* dims, int rank, int* coords, int* xs, int* xe, int* ys, int* ye) {
    int xprocs = dims[0];
    int yprocs = dims[1];

    int xepp = nx / xprocs; //elements per process without remainder
    int xr = nx % xprocs; //remainder

    *xs = coords[0] * xepp + 1;
    *xs = *xs + ((coords[0] < xr) ? coords[0] : xr);
    if (coords[0] < xr) xepp++;
    *xe = *xs + xepp - 1;
    if (*xe > nx || coords[0] == xprocs - 1) *xe = nx; 

    int yepp = ny / yprocs; //elements per process without remainder
    int yr = ny % yprocs; //remainder

    *ys = coords[1] * yepp + 1;
    *ys = *ys + ((coords[1] < yr) ? coords[1] : yr);
    if (coords[1] < yr) yepp++;
    *ye = *ys + yepp - 1;
    if (*ye > ny || coords[1] == yprocs - 1) *ye = ny;

    return MPI_SUCCESS;
}
