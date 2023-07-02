#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv){

    char *vector = argv[1];
    int n = atoi(argv[1]);

    char filename[] = "vector_";
    strcat(filename,vector);
    strcat(filename,".txt");


    FILE *fp = fopen(filename,"w");
    fprintf(fp,"%d\n",n);
    for (int i = 1; i < n+1; i++){
	fprintf(fp,"%d\n",i);
    }

    fclose(fp);

    return 0;
}
