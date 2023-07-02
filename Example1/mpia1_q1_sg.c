#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv){
    char str[100];
    
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size == 1){
	fprintf(stderr,"Error: only one process");
	return 1;
    }


    // Read the vector from the file
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


    int epp = N/size; //elements per process expect rank 0
    double *gathered_vector = malloc(N * sizeof(double));  // Allocate memory for the gathered vector on rank 0
    double *recvbuf = malloc(epp * sizeof(double));  // Allocate memory for the vector on all ranks

    // Send the vector elements to all ranks
    MPI_Scatter(vector, epp, MPI_DOUBLE, recvbuf, epp, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    for(int i = 0; i < size; i++){
	if(rank == i){
	    for(int k = 0; k < epp; k++){
		recvbuf[k] += 1.0;
		//printf("rank = %d, recvb = %f\n", i ,recvbuf[k]);
	    }
	}
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Gather the local vectors back into a vector on rank 0
    MPI_Gather(recvbuf, epp, MPI_DOUBLE, gathered_vector, epp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
    if (rank == 0){   
    // Write the gathered vector to a file
        FILE* fp = fopen("gathered_vector_sg.txt", "w");
        fprintf(fp, "%d\n",N);
        for (int i = 0; i < N; i++) {
	    fprintf(fp, "%d\n", (int)gathered_vector[i]);
        }
        fclose(fp);
    }

    free(vector);
    free(recvbuf);
    free(gathered_vector);

    MPI_Finalize();


    return 0;
}

