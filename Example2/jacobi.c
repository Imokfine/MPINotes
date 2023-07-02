#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "poisson1d.h"
#include "jacobi.h"

void sweep1d(double a[][maxn], double f[][maxn], int nx,
    int s, int e, double b[][maxn])
{
    double h;
    int i, j;

    h = 1.0 / ((double)(nx + 1));

    for (i = s; i <= e; i++) {
        for (j = 1; j < nx + 1; j++) {
            b[i][j] = 0.25 * (a[i - 1][j] + a[i + 1][j] + a[i][j + 1] + a[i][j - 1] - h * h * f[i][j]);
        }
    }
}

void sweep2d(double a[][maxn], double f[][maxn], int nx,
    int* xs, int* xe, int* ys, int* ye, double b[][maxn])
{
    double h;
    int i, j;

    h = 1.0 / ((double)(nx + 1));

    for (i = *xs; i <= *xe; i++) {
        for (j = *ys; j <= *ye; j++) {
            b[i][j] = 0.25 * (a[i - 1][j] + a[i + 1][j] + a[i][j + 1] + a[i][j - 1] - h * h * f[i][j]);

	}
    }

}

/* sequentialized if there is no buffering */
void exchang1(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
    int nbrleft, int nbrright)
{
    int rank;
    int ny;

    MPI_Comm_rank(comm, &rank);

    ny = nx;

    MPI_Ssend(&x[e][1], ny, MPI_DOUBLE, nbrright, 0, comm);
    /* printf("(myid: %d) sent \"col\" %d with %d entries to nbr: %d\n",rank, e, ny, nbrright); */

    MPI_Recv(&x[s - 1][1], ny, MPI_DOUBLE, nbrleft, 0, comm, MPI_STATUS_IGNORE);
    /* printf("(myid: %d) recvd into \"col\" %d from %d entries from nbr: %d\n",rank, s-1, ny, nbrleft); */

    MPI_Ssend(&x[s][1], ny, MPI_DOUBLE, nbrleft, 1, comm);
    /* printf("(myid: %d) sent \"col\" %d with %d entries to nbr: %d\n",rank, s, ny, nbrleft); */
    MPI_Recv(&x[e + 1][1], ny, MPI_DOUBLE, nbrright, 1, comm, MPI_STATUS_IGNORE);
    /* printf("(myid: %d) recvd into \"col\" %d from %d entries from nbr: %d\n",rank, e+1, ny, nbrright); */

}

/* ordered sends / receives */
void exchang2(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
    int nbrleft, int nbrright)
{
    int coord;
    int rank;

    MPI_Comm_rank(comm, &rank);

    coord = rank;

    if (coord % 2 == 0) {

        MPI_Ssend(&x[e][1], nx, MPI_DOUBLE, nbrright, 0, comm);

        MPI_Recv(&x[s - 1][1], nx, MPI_DOUBLE, nbrleft, 0, comm, MPI_STATUS_IGNORE);

        MPI_Ssend(&x[s][1], nx, MPI_DOUBLE, nbrleft, 1, comm);

        MPI_Recv(&x[e + 1][1], nx, MPI_DOUBLE, nbrright, 1, comm, MPI_STATUS_IGNORE);
    }
    else {
        MPI_Recv(&x[s - 1][1], nx, MPI_DOUBLE, nbrleft, 0, comm, MPI_STATUS_IGNORE);

        MPI_Ssend(&x[e][1], nx, MPI_DOUBLE, nbrright, 0, comm);

        MPI_Recv(&x[e + 1][1], nx, MPI_DOUBLE, nbrright, 1, comm, MPI_STATUS_IGNORE);

        MPI_Ssend(&x[s][1], nx, MPI_DOUBLE, nbrleft, 1, comm);

    }

}

/* sendrecv */
void exchang3(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
    int nbrleft, int nbrright)
{

    MPI_Sendrecv(&x[e][1], nx, MPI_DOUBLE, nbrright, 0, &x[s - 1][1], nx, MPI_DOUBLE, nbrleft,
        0, comm, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&x[s][1], nx, MPI_DOUBLE, nbrleft, 1, &x[e + 1][1], nx, MPI_DOUBLE, nbrright,
        1, comm, MPI_STATUS_IGNORE);

}

void exchang3_2d(double x[][maxn], int nx, int ny, int nprocs, int* dims, int* coords, int* xs, int* xe, int* ys, int* ye, int rank, MPI_Comm comm,
    int nbrleft, int nbrright, int nbrup, int nbrdown)
{
    MPI_Sendrecv(&x[*xe][*ys], *ye - *ys + 1, MPI_DOUBLE, nbrright, 0, &x[*xs - 1][*ys], *ye - *ys + 1, MPI_DOUBLE, nbrleft,
        0, comm, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&x[*xs][*ys], *ye - *ys + 1, MPI_DOUBLE, nbrleft, 1, &x[*xe + 1][*ys], *ye - *ys + 1, MPI_DOUBLE, nbrright,
        1, comm, MPI_STATUS_IGNORE);

    MPI_Datatype coltype;
    MPI_Type_vector(*xe - *xs + 1, 1, ny + 2, MPI_DOUBLE, &coltype);
    MPI_Type_commit(&coltype);

    MPI_Sendrecv(&x[*xs][*ye], 1, coltype, nbrup, 3, &x[*xs][*ys - 1], 1, coltype, nbrdown,
            3, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&x[*xs][*ys], 1, coltype, nbrdown, 2, &x[*xs][*ye + 1], 1, coltype, nbrup,
            2, comm, MPI_STATUS_IGNORE);
}

void exchangi1(double x[][maxn], int nx, int s, int e, MPI_Comm comm,
    int nbrleft, int nbrright)
{
    MPI_Request reqs[4];

    MPI_Irecv(&x[s - 1][1], nx, MPI_DOUBLE, nbrleft, 0, comm, &reqs[0]);
    MPI_Irecv(&x[e + 1][1], nx, MPI_DOUBLE, nbrright, 0, comm, &reqs[1]);
    MPI_Isend(&x[e][1], nx, MPI_DOUBLE, nbrright, 0, comm, &reqs[2]);
    MPI_Isend(&x[s][1], nx, MPI_DOUBLE, nbrleft, 0, comm, &reqs[3]);
    /* not doing anything useful here */

      /* not doing anything useful here */

    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
}

void exchangi2d(double x[][maxn], int nx, int ny, int* coords, int* xs, int* xe, int* ys, int* ye, MPI_Comm comm,
    int nbrleft, int nbrright, int nbrup, int nbrdown)
{

    if (coords[0] % 2 == 0) {
        MPI_Request req1[4];

        MPI_Irecv(&x[*xs - 1][*ys], nx, MPI_DOUBLE, nbrleft, 0, comm, &req1[0]);
        MPI_Irecv(&x[*xe + 1][*ys], nx, MPI_DOUBLE, nbrright, 0, comm, &req1[1]);
        MPI_Isend(&x[*xe][*ys], nx, MPI_DOUBLE, nbrright, 0, comm, &req1[2]);
        MPI_Isend(&x[*xs][*ys], nx, MPI_DOUBLE, nbrleft, 0, comm, &req1[3]);

        MPI_Waitall(4, req1, MPI_STATUSES_IGNORE);
    }
    else {

        MPI_Request req2[4];

        MPI_Isend(&x[*xe][*ys], nx, MPI_DOUBLE, nbrright, 0, comm, &req2[0]);
        MPI_Isend(&x[*xs][*ys], nx, MPI_DOUBLE, nbrleft, 0, comm, &req2[1]);
        MPI_Irecv(&x[*xs - 1][*ys], nx, MPI_DOUBLE, nbrleft, 0, comm, &req2[2]);
        MPI_Irecv(&x[*xe + 1][*ys], nx, MPI_DOUBLE, nbrright, 0, comm, &req2[3]);

        MPI_Waitall(4, req2, MPI_STATUSES_IGNORE);

    }

    MPI_Datatype coltype;
    MPI_Type_vector(*xe - *xs + 1, 1, ny + 2, MPI_DOUBLE, &coltype);
    MPI_Type_commit(&coltype);

    if (coords[1] % 2 == 0) {
        MPI_Request req3[4];

        MPI_Isend(&x[*xs][*ys], 1, coltype, nbrdown, 2, comm, &req3[0]);
        MPI_Irecv(&x[*xs][*ye + 1], 1, coltype, nbrup, 2, comm, &req3[1]);
        MPI_Isend(&x[*xs][*ye], 1, coltype, nbrup, 3, comm, &req3[2]);
        MPI_Irecv(&x[*xs][*ys - 1], 1, coltype, nbrdown, 3, comm, &req3[3]);

        MPI_Waitall(4, req3, MPI_STATUSES_IGNORE);

    }
    else {
        MPI_Request req4[4];

        MPI_Irecv(&x[*xs][*ye + 1], 1, coltype, nbrup, 2, comm, &req4[0]);
        MPI_Isend(&x[*xs][*ys], 1, coltype, nbrdown, 2, comm, &req4[1]);
        MPI_Irecv(&x[*xs][*ys - 1], 1, coltype, nbrdown, 3, comm, &req4[2]);
        MPI_Isend(&x[*xs][*ye], 1, coltype, nbrup, 3, comm, &req4[3]);

        MPI_Waitall(4, req4, MPI_STATUSES_IGNORE);

    }
}


void nbxchange_and_sweep(double u[][maxn], double f[][maxn], int nx, int ny,
    int s, int e, double unew[][maxn], MPI_Comm comm,
    int nbrleft, int nbrright)
{
    MPI_Request req[4];
    MPI_Status status;
    int idx;
    double h;
    int i, j, k;

    int myid;
    MPI_Comm_rank(comm, &myid);

    h = 1.0 / ((double)(nx + 1));

    /* int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, */
    /*               int source, int tag, MPI_Comm comm, MPI_Request *request); */
    /* int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, */
    /*            int tag, MPI_Comm comm, MPI_Request *request); */

    MPI_Irecv(&u[s - 1][1], ny, MPI_DOUBLE, nbrleft, 1, comm, &req[0]);
    MPI_Irecv(&u[e + 1][1], ny, MPI_DOUBLE, nbrright, 2, comm, &req[1]);

    MPI_Isend(&u[e][1], ny, MPI_DOUBLE, nbrright, 1, comm, &req[2]);
    MPI_Isend(&u[s][1], ny, MPI_DOUBLE, nbrleft, 2, comm, &req[3]);

    /* perform purely local updates (that don't need ghosts) */
    /* 2 cols or less means all are on processor boundary */
    if (e - s + 1 > 2) {
        for (i = s + 1; i < e; i++) {
            for (j = 1; j < ny + 1; j++) {
                unew[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j + 1] + u[i][j - 1] - h * h * f[i][j]);
            }
        }
    }

    /* perform updates in j dir only for boundary cols */
    for (j = 1; j < ny + 1; j++) {
        unew[s][j] = 0.25 * (u[s][j + 1] + u[s][j - 1] - h * h * f[s][j]);
        unew[e][j] = 0.25 * (u[e][j + 1] + u[e][j - 1] - h * h * f[e][j]);
    }

    /* int MPI_Waitany(int count, MPI_Request array_of_requests[], */
    /*      int *index, MPI_Status *status) */
    for (k = 0; k < 4; k++) {
        MPI_Waitany(4, req, &idx, &status);

        /* idx 0, 1 are recvs */
        switch (idx) {
        case 0:
            /* printf("myid: %d case idx 0: status.MPI_TAG: %d; status.MPI_SOURCE: %d (idx: %d)\n",myid,status.MPI_TAG, status.MPI_SOURCE,idx); */
            if (nbrleft != MPI_PROC_NULL &&
                (status.MPI_TAG != 1 || status.MPI_SOURCE != nbrleft)) {
                fprintf(stderr, "Error: I don't understand the world: (tag %d; source %d)\n",
                    status.MPI_TAG, status.MPI_SOURCE);
                MPI_Abort(comm, 1);
            }

            /* left ghost update completed; update local leftmost column */
            for (j = 1; j < ny + 1; j++) {
                unew[s][j] += 0.25 * (u[s - 1][j]);
            }
            break;
        case 1:
            /* printf("myid: %d case idx 1: status.MPI_TAG: %d; status.MPI_SOURCE: %d (idx: %d)\n",myid, status.MPI_TAG, status.MPI_SOURCE,idx); */
            if (nbrright != MPI_PROC_NULL &&
                (status.MPI_TAG != 2 || status.MPI_SOURCE != nbrright)) {
                fprintf(stderr, "Error: I don't understand the world: (tag %d; source %d)\n",
                    status.MPI_TAG, status.MPI_SOURCE);
                MPI_Abort(comm, 1);
            }
            /* right ghost update completed; update local rightmost
 *                column */
            for (j = 1; j < ny + 1; j++) {
                unew[e][j] += 0.25 * (u[e + 1][j]);
            }
            break;
        default:
            break;
        }
    }
    /* splitting this off to take account of case of one column assigned
 *        to proc -- so left and right node neighbours are ghosts so both
 *               the recvs must be complete*/
    for (j = 1; j < ny + 1; j++) {
        unew[s][j] += 0.25 * (u[s + 1][j]);
        unew[e][j] += 0.25 * (u[e - 1][j]);
    }

}

double griddiff(double a[][maxn], double b[][maxn], int nx, int s, int e)
{
    double sum;
    double tmp;
    int i, j;

    sum = 0.0;

    for (i = s; i <= e; i++) {
        for (j = 1; j < nx + 1; j++) {
            tmp = (a[i][j] - b[i][j]);
            sum = sum + tmp * tmp;
        }
    }

    return sum;

}

double griddiff2d(double a[][maxn], double b[][maxn], int nx, int* xs, int* xe, int* ys, int* ye)
{
    double sum;
    double tmp;
    int i, j;

    sum = 0.0;

    for (i = *xs; i <= *xe; i++) {
        for (j = *ys; j < *ye; j++) {
            tmp = (a[i][j] - b[i][j]);
            sum = sum + tmp * tmp;
        }
    }

    return sum;

}

