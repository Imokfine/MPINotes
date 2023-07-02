#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
int decomp1d(int n, int p, int myid, int *s, int *e);

int main(int argc, char** argv){
    char str[100];
    
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FILE *file = fopen(argv[1],"r");

    if (file == NULL){
	fprintf(stderr, "Error: unable to open file\n");
	return 1;
    }

    int N = atoi(fgets(str,100,file));
    double *vector = malloc(N * sizeof(double));

    int i = 0; 

    while(fgets(str,100,file)!=NULL){
	vector[i] = atoi(str);
	i++;
    }
	
    fclose(file);


    double *gathered_vector = malloc(N * sizeof(double));  // Allocate memory for the gathered vector on rank 0
    double *gathered_vector2 = malloc(N * sizeof(double));  // Allocate memory for the gathered vector on rank 2
    double *recvbuf = malloc((N/size+1) * sizeof(double));  // Allocate memory for the vector on all ranks
    double *sl2 = malloc(sizeof(double));
    double *sl2_2 = malloc(sizeof(double));
    int *s = malloc(sizeof(int));
    int *e = malloc(sizeof(int));
    int *epp = malloc(sizeof(int));


    // Send the vector elements to all ranks
    if (rank == 0){
	for (int i = 1; i < size; i++){
	    decomp1d(N, size, i, s, e);
	    *epp = *e - *s + 1;
	    MPI_Send(&vector[*s], *epp, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
	}
	
	decomp1d(N, size, 0, s, e);
	*epp = *e - *s + 1;
    
	// calculate the square of the Euclidean (l2) norm of the vector	
	*sl2 = 0;
 	for (int j = 0; j < *epp; j++){
	    recvbuf[j] = vector[j];
	    *sl2 += recvbuf[j] * recvbuf[j];
	}
	*sl2_2 = *sl2;
//	printf("rank = %d, sl2 = %f\n",rank, *sl2);
    }

    else{
	decomp1d(N, size, rank, s, e);
	*epp = *e - *s + 1;
	MPI_Recv(recvbuf, *epp, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    	
	// calculate the square of the Euclidean (l2) norm of the vector	
	*sl2 = 0;
	for (int j = 0; j < *epp; j++){
	    *sl2 += recvbuf[j] * recvbuf[j];
	}

	*sl2_2 = *sl2;
//	printf("rank = %d, sl2 = %f\n",rank, *sl2);
    }

    MPI_Barrier(MPI_COMM_WORLD);



    double sum1 = 0;
    // Gather the local vectors back into a vector on rank 0
    MPI_Gather(sl2, 1, MPI_DOUBLE, gathered_vector, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0){
	for (int i =0; i < N; i++)
	    sum1 += gathered_vector[i];
	printf("full sum = %d (rank 0)\n",(int)sum1);
    }


    // Gather the local vectors back into a vector on rank 2
    if (size < 3){
	printf("the number of processors is less than 3, so there is no rank 2\n");	
    }

    else{
    if (rank != 2){
	MPI_Send(sl2_2, 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD);
//	printf("rank = %d, sl2 = %f\n",rank, *sl2_2);
    }

    else{
	for (int i = 0; i < size; i++){
	    if (i != 2){
		MPI_Recv(&gathered_vector2[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    }
	}
    }

    double sum2 = 0;
    if (rank == 2){
	gathered_vector2[2] = *sl2_2;
	for (int i = 0; i < size; i++){
//	    printf("i = %d, gathered_v = %f\n",i ,gathered_vector2[i]);
	    sum2 += gathered_vector2[i];
	}
	printf("full sum = %d (rank 2)\n",(int)sum2);
    }    
    }

    free(vector);
    free(recvbuf);
    free(gathered_vector);
    free(sl2);
    free(epp);
    free(gathered_vector2);
    free(sl2_2);

    MPI_Finalize();


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
