#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "jacobiRMA.h"

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

void exchangrma2d(double x[][maxn], int ny, int* xs, int* xe, int* ys, int* ye, int rank, MPI_Comm comm,
    int nbrleft, int nbrright, int nbrup, int nbrdown, MPI_Win win) {

    MPI_Datatype coltype;
    MPI_Type_vector(*xe - *xs + 1, 1, ny + 2, MPI_DOUBLE, &coltype);
    MPI_Type_commit(&coltype);

    MPI_Aint offset;

	MPI_Win_fence(0,win);

	// right
	offset = (MPI_Aint)(maxn + *ys);
	MPI_Get(&x[*xe + 1][*ys], *ye - *ys + 1, MPI_DOUBLE, nbrright, offset, *ye - *ys + 1, MPI_DOUBLE, win);

	// up
	offset = (MPI_Aint)(maxn + *ye + 1);//We want to move 1 column in hence an offset of max n. We then want to get from e[1] + 1 up the column.
	MPI_Get(&x[*xs][*ye + 1], 1, coltype, nbrup, offset, 1, coltype, win);

	// left
	offset = (MPI_Aint)(*ys);
	MPI_Put(&x[*xe][*ys], *ye - *ys + 1, MPI_DOUBLE, nbrright, offset, *ye - *ys + 1, MPI_DOUBLE, win);

	// down
	offset = (MPI_Aint)(maxn + *ye);
	MPI_Put(&x[*xs][*ye], 1, coltype, nbrup, offset, 1, coltype,win);

	MPI_Win_fence(0, win);

}

void exchangrma_2d_pscw(double x[][maxn], int ny, int* xs, int* xe, int* ys, int* ye, MPI_Win win,
	int nbrleft, int nbrright, int nbrup, int nbrdown, MPI_Group origin[2], MPI_Group target[2]) {

	// Origin process calls MPI_Win_start and MPI_Win_complete
	// Target process calls MPI_Win_post and MPI_Win_wait


	MPI_Datatype coltype;
	MPI_Type_vector(*xe - *xs + 1, 1, ny + 2, MPI_DOUBLE, &coltype);
	MPI_Type_commit(&coltype);

	MPI_Aint offset;

	if (nbrleft != MPI_PROC_NULL) {
		MPI_Win_post(target[0], 0, win);
		MPI_Win_wait(win);
	}

	if (nbrdown != MPI_PROC_NULL) {
		MPI_Win_post(target[1], 0, win);
		MPI_Win_wait(win);
	}

	if (nbrright != MPI_PROC_NULL) {
		MPI_Win_start(origin[0], 0, win);
		offset = (MPI_Aint)(maxn + *ys);
		MPI_Get(&x[*xe + 1][*ys], *ye - *ys + 1, MPI_DOUBLE, nbrright, offset, *ye - *ys + 1, MPI_DOUBLE, win);

		offset = (MPI_Aint)(*ys);
		MPI_Put(&x[*xe][*ys], *ye - *ys + 1, MPI_DOUBLE, nbrright, offset, *ye - *ys + 1, MPI_DOUBLE, win);
		MPI_Win_complete(win);
	}

	if (nbrup != MPI_PROC_NULL) {
		MPI_Win_start(origin[1], 0, win);
		offset = (MPI_Aint)(maxn + *ye + 1);
		MPI_Get(&x[*xs][*ye + 1], 1, coltype, nbrup, offset, 1, coltype, win);

		offset = (MPI_Aint)(maxn + *ye);
		MPI_Put(&x[*xs][*ye], 1, coltype, nbrup, offset, 1, coltype, win);
		MPI_Win_complete(win);
	}



}

double griddiff2d(double a[][maxn], double b[][maxn], int* xs, int* xe, int* ys, int* ye)
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
