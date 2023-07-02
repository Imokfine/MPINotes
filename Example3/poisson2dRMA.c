#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include <mpi.h>

#include "poisson.h"
#include "jacobiRMA.h"

#define maxit 10000

#include "decomp.h"

void init_full_grid(double g[][maxn]);
void init_full_grids(double a[][maxn], double b[][maxn], double f[][maxn]);

void twodinit_basic(double a[][maxn], double b[][maxn], double f[][maxn],
    int nx, int ny, int* xs, int* xe, int* ys, int* ye);

void print_full_grid(double x[][maxn]);
void print_in_order(double x[][maxn], MPI_Comm comm);
void print_grid_to_file(char* fname, double x[][maxn], int nx, int ny);

/* Definition of write_grid2d() and GatherGrid2d() */
void write_grid2d(double a[][maxn], int nx, int* xs, int* xe, int* ys, int* ye, int rank);
void GatherGrid2d(double a[][maxn], int nx, int ny, int nprocs, int* dims, int* coords, int rank, MPI_Comm comm);

int main(int argc, char** argv)
{
    double a[maxn][maxn], b[maxn][maxn], f[maxn][maxn];
    int nx, ny;
    int myid, nprocs;
    /* MPI_Status status; */
    int nbrleft, nbrright, nbrup, nbrdown;
    int* xs = (int*)malloc(sizeof(int));
    int* xe = (int*)malloc(sizeof(int));
    int* ys = (int*)malloc(sizeof(int));
    int* ye = (int*)malloc(sizeof(int));
    int it;
    double glob_diff;
    double ldiff;
    double t1, t2;
    double tol = 1.0E-11;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (myid == 0) {
        /* set the size of the problem */
        if (argc > 2) {
            fprintf(stderr, "---->Usage: mpirun -np <nproc> %s <nx>\n", argv[0]);
            fprintf(stderr, "---->(for this code nx=ny)\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (argc == 2) {
            nx = atoi(argv[1]);
        }
        if (argc == 1) {
            nx = 15;
        }

        if (nx > maxn - 2) {
            fprintf(stderr, "grid size too large\n");
            exit(1);
        }
    }

    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("(myid: %d) nx = %d\n", myid, nx);
    ny = nx;

    init_full_grids(a, b, f);

    /*q4: Use MPI Cartesian topology functions such as MPI_Cart_create()
          and MPI_Cart_shift to calculate the up, down, left and right
          neighbours of each process */

    MPI_Comm comm_cart;
    int ndim = 2;
    int dims[2] = { 0, 0 };
    int periods[2] = { 0, 0 };
    int coords[2];

    /* Set the dims of processors */
    MPI_Dims_create(nprocs, ndim, dims);

    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, 0, &comm_cart);
    MPI_Cart_shift(comm_cart, 0, 1, &nbrleft, &nbrright);
    MPI_Cart_shift(comm_cart, 1, 1, &nbrdown, &nbrup);

    /* Get the coordinates of rank */
    MPI_Cart_coords(comm_cart, myid, ndim, coords);

    MPE_Decomp2d(nx, ny, dims, myid, coords, xs, xe, ys, ye);

    printf("(myid: %d coord: %d, %d) nx: %d; ny: %d; xs: %d; xe: %d; ys: %d; ye: %d; nbrleft: %d; nbrright: %d;  nbrup: %d; nbrdown: %d\n",
        myid, coords[0], coords[1], nx, ny, *xs, *xe, *ys, *ye, nbrleft, nbrright, nbrup, nbrdown);

    twodinit_basic(a, b, f, nx, ny, xs, xe, ys, ye);

    MPI_Win wina, winb;

    MPI_Win_create(&a[*xs - 1][0], (maxn) * (*xe - *xs + 3) * sizeof(double), sizeof(double),
        MPI_INFO_NULL, MPI_COMM_WORLD, &wina);
    MPI_Win_create(&b[*xs - 1][0], (maxn) * (*xe - *xs + 3) * sizeof(double), sizeof(double),
        MPI_INFO_NULL, MPI_COMM_WORLD, &winb);

    t1 = MPI_Wtime();

    glob_diff = 1000;
    for (it = 0; it < maxit; it++) {


        exchangrma2d(a, ny, xs, xe, ys, ye, myid, MPI_COMM_WORLD, nbrleft, nbrright, nbrup, nbrdown, wina);
        sweep2d(a, f, nx, xs, xe, ys, ye, b);



        exchangrma2d(b, ny, xs, xe, ys, ye, myid, MPI_COMM_WORLD, nbrleft, nbrright, nbrup, nbrdown, winb);
        sweep2d(b, f, nx, xs, xe, ys, ye, a);

        ldiff = griddiff2d(a, b, xs, xe, ys, ye);
        MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (myid == 0 && it % 10 == 0) {
            printf("(myid %d) locdiff: %lf; glob_diff: %lf\n", myid, ldiff, glob_diff);
        }
        if (glob_diff < tol) {
            if (myid == 0) {
                printf("iterative solve converged\n");
            }
            break;
        }

    }

    t2 = MPI_Wtime();

    printf("DONE! (it: %d)\n", it);

    if (myid == 0) {
        if (it == maxit) {
            fprintf(stderr, "Failed to converge\n");
        }
        printf("Run took %lf s\n", t2 - t1);
    }

    print_in_order(a, MPI_COMM_WORLD);
    if (nprocs == 1) {
        print_grid_to_file("grid", a, nx, ny);
        print_full_grid(a);
    }

    GatherGrid2d(a, nx, ny, nprocs, dims, coords, myid, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myid == 0) {
        printf(" The grid on root is \n");
        print_full_grid(a);

        *xs = 1;
        *xe = nx;
        *ys = 1;
        *ye = ny;
        write_grid2d(a, nx, xs, xe, ys, ye, myid);
    }

    free(xs);
    free(xe);
    free(ys);
    free(ye);
    MPI_Win_free(&wina);
    MPI_Win_free(&winb);
    MPI_Comm_free(&comm_cart);
    MPI_Finalize();
    return 0;
}

void twodinit_basic(double a[][maxn], double b[][maxn], double f[][maxn],
    int nx, int ny, int* xs, int* xe, int* ys, int* ye)
{
    int i, j;


    double left, bottom, right, top;

    left = -1.0;
    bottom = 1.0;
    right = 2.0;
    top = 3.0;

    /* set everything to 0 first */
    for (i = *xs - 1; i <= *xe + 1; i++) {
        for (j = *ys - 1; j <= *ye + 1; j++) {
            a[i][j] = 0.0;
            b[i][j] = 0.0;
            f[i][j] = 0.0;
        }
    }

    /* deal with boundaries */
    /* Set the left boundary */
    if (*xs == 1) {
        for (i = *ys; i <= *ye; i++) {
            a[0][i] = left;
            b[0][i] = left;
        }
    }

    /* Set the right boundary */
    if (*xe == ny) {
        for (i = *ys; i <= *ye; i++) {
            a[ny + 1][i] = right;
            b[ny + 1][i] = right;
        }
    }
    /* Set the bottom boundary */
    if (*ys == 1) {
        for (i = *xs; i <= *xe; i++) {
            a[i][0] = bottom;
            b[i][0] = bottom;
        }
    }

    /* Set the top boundary */
    if (*ye == ny) {
        for (i = *xs; i <= *xe; i++) {
            a[i][nx + 1] = top;
            b[i][nx + 1] = top;
        }
    }
}

void init_full_grid(double g[][maxn])
{
    int i, j;
    const double junkval = -5;

    for (i = 0; i < maxn; i++) {
        for (j = 0; j < maxn; j++) {
            g[i][j] = junkval;
        }
    }
}

/* set global a,b,f to initial arbitrarily chosen junk value */
void init_full_grids(double a[][maxn], double b[][maxn], double f[][maxn])
{
    int i, j;
    const double junkval = -5;

    for (i = 0; i < maxn; i++) {
        for (j = 0; j < maxn; j++) {
            a[i][j] = junkval;
            b[i][j] = junkval;
            f[i][j] = junkval;
        }
    }

}

/* prints to stdout in GRID view */
void print_full_grid(double x[][maxn])
{
    int i, j;
    for (j = maxn - 1; j >= 0; j--) {
        for (i = 0; i < maxn; i++) {
            if (x[i][j] < 10000.0) {
                printf("|%2.6lf| ", x[i][j]);
            }
            else {
                printf("%9.2lf ", x[i][j]);
            }
        }
        printf("\n");
    }

}

void print_in_order(double x[][maxn], MPI_Comm comm)
{
    int myid, size;
    int i;

    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &size);
    MPI_Barrier(comm);
    printf("Attempting to print in order\n");
    sleep(1);
    MPI_Barrier(comm);

    for (i = 0; i < size; i++) {
        if (i == myid) {
            printf("proc %d\n", myid);
            print_full_grid(x);
        }
        fflush(stdout);
        usleep(500);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void print_grid_to_file(char* fname, double x[][maxn], int nx, int ny)
{
    FILE* fp;
    int i, j;

    fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "Error: can't open file %s\n", fname);
        exit(4);
    }

    for (j = ny + 1; j >= 0; j--) {
        for (i = 0; i < nx + 2; i++) {
            fprintf(fp, "%lf ", x[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}


/* This function allows each processor to write it’s part of the grid to file or stdout.
 *    The output should be in “mesh/grid” format*/
void write_grid2d(double a[][maxn], int nx, int* xs, int* xe, int* ys, int* ye, int rank) {

    char filename[50];
    sprintf(filename, "RMA2d_Grid_of_rank_%d.txt", rank);

    FILE* fp;
    fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: can't open file %s\n", filename);
        exit(-1);
    }

    for (int j = *ye; j >= 1; j--) {
        for (int i = *xs; i <= *xe; i++) {
            fprintf(fp, "%lf ", a[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/* This function gathers the parallel solution which
 *    is distributed across multiple processors onto one processor (rank 0)*/
void GatherGrid2d(double a[][maxn], int nx, int ny, int nprocs, int* dims, int* coords, int rank, MPI_Comm comm) {
    int xprocs = dims[0];
    int yprocs = dims[1];

    int* lxs = (int*)malloc(sizeof(int));
    int* lxe = (int*)malloc(sizeof(int));
    int* lys = (int*)malloc(sizeof(int));
    int* lye = (int*)malloc(sizeof(int));

    int n_cols[nprocs];
    int n_rows[nprocs];
    int ixs[nprocs];
    int ixe[nprocs];
    int iys[nprocs];
    int iye[nprocs];

    int rx = nx % xprocs;
    int ry = ny % yprocs;

    for (int i = 0; i < nprocs; i++) {

        n_cols[i] = nx / xprocs;
        n_rows[i] = ny / yprocs;

        if (coords[0] < rx) {
            n_cols[i]++;
        }
        if (coords[1] < ry) {
            n_rows[i]++;
        }

    }


    MPE_Decomp2d(nx, ny, dims, rank, coords, lxs, lxe, lys, lye);
    ixs[rank] = *lxs;
    ixe[rank] = *lxe;
    iys[rank] = *lys;
    iye[rank] = *lye;

    MPI_Barrier(comm);

    if (rank != 0) {
        MPI_Send(&ixs[rank], 1, MPI_INT, 0, rank * 20000 + 1, comm);
        MPI_Send(&ixe[rank], 1, MPI_INT, 0, rank * 20000 + 2, comm);
        MPI_Send(&iys[rank], 1, MPI_INT, 0, rank * 20000 + 3, comm);
        MPI_Send(&iye[rank], 1, MPI_INT, 0, rank * 20000 + 4, comm);
    }

    if (rank == 0) {
        for (int i = 1; i < nprocs; i++) {
            MPI_Recv(&ixs[i], 1, MPI_INT, i, i * 20000 + 1, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&ixe[i], 1, MPI_INT, i, i * 20000 + 2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&iys[i], 1, MPI_INT, i, i * 20000 + 3, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&iye[i], 1, MPI_INT, i, i * 20000 + 4, comm, MPI_STATUS_IGNORE);
        }
    }

    MPI_Barrier(comm);
    /*
        if (rank == 0){
        for(int i = 0; i < nprocs; i++)
               printf("xs[%d] = %d, xe[%d] = %d, ys[%d] = %d, ye[%d] = %d \n", i, ixs[i], i, ixe[i], i, iys[i], i, iye[i]);
        }
    */

    if (rank != 0) {
        for (int col = ixs[rank]; col <= ixe[rank]; col++) {
            MPI_Send(&a[col][iys[rank]], n_rows[rank], MPI_DOUBLE, 0, rank * 10000 + col, comm);
        }

    }

    if (rank == 0) {
        for (int i = 1; i < nprocs; i++) {
            for (int col = ixs[i]; col <= ixe[i]; col++) {
                MPI_Recv(&a[col][iys[i]], n_rows[i], MPI_DOUBLE, i, i * 10000 + col, comm, MPI_STATUS_IGNORE);
            }
        }
    }

    free(lxs);
    free(lxe);
    free(lys);
    free(lye);
}

