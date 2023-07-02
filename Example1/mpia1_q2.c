#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
int decomp1d(int n, int p, int myid, int *s, int *e);

int main(int argc,char **argv){

    int n;
    int *s = malloc(sizeof(int));
    int *e = malloc(sizeof(int));
    n = atoi(argv[1]);

    MPI_Init(&argc, &argv);

    int p, myid;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    for (int i = 0; i < p; i++){
	if (myid == i){
	decomp1d(n,p,i,s,e);
	printf("rank = %d, start = %d, end = %d\n", i, *s, *e);
	}
    }

    MPI_Finalize();

    free(s);
    free(e);

    return 0;
}

int decomp1d(int n, int p, int myid, int *s, int *e){
    int epp = n/p; //elements per process without remainder
    int r = n%p; //remainder


    if (r == 0 ){
	*s = myid * epp;
	*e = *s + epp - 1;
    }

    else{
	if (myid < r){
	    *s = myid * (epp + 1);
	    *e = *s + epp;
	}
	else{
	    *s = myid * epp + r;
	    *e = *s + epp - 1;
	}
    }
}
