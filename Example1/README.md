## Requirements
    (1)To reproduce the project codes, it is necessary to set mpi environment on chuck. Execute the code as follows:
    [ user@chuck ~] $ module load cports openmpi


## Usage
    (2) Use makefile to compile all the c files as follows:
    [ user@chuck ~]$ make

-------------------------------------------------------------------------------------------------------------------

    For question 1, there are two programs.
    One is mpia1_q1_sr, which uses MPI_Send/MPI_Recv to implement.
    Another is mpia1_q1_sg, which uses MPI collective calls to implement.

    (3) Create a vector of length n = 16 as follows:
    [ user@chuck ~]$ ./create_vector 16

    (4) Create a vector of length n = 10000 as follows:
    [ user@chuck ~]$ ./create_vector 10000

    (5) To run the code on 4 processors and for a vector of length n = 16, execute:
    [ user@chuck ~]$ mpirun -n 4 ./mpia1_q1_sr vector_16.txt
    or
    [ user@chuck ~]$ mpirun -n 4 ./mpia1_q1_sg vector_16.txt
    and then check the results in "gathered_vector_sr.txt" or "gathered_vector_sg.txt"

    (6) To run the code on at least 10 processors for a vector length at n = 10000, execute:
    [ user@chuck ~]$ mpirun -n 10 ./mpia1_q1_sr vector_10000.txt
    or
    [ user@chuck ~]$ mpirun -n 10 ./mpia1_q1_sg vector_10000.txt
    and then check the results in "gathered_vector_sr.txt" or "gathered_vector_sg.txt"

-------------------------------------------------------------------------------------------------------------------

    For question 2, the program is named mpia1_q2.
    
    (7) To demonstrate that the decomposition produces consistent results on up to 9 processors for n = 25, execute:
    [ user@chuck ~]$ mpirun -n 9 ./mpia1_q2 25

    (8) To demonstrate that the decomposition produces consistent results on up to 10 processors for n = 100000, execute:
    [ user@chuck ~]$ mpirun -n 10 ./mpia1_q2 100000

    (9) To demonstrate that the function also produces correct results when run on one processor, execute:
    [ user@chuck ~]$ mpirun -n 1 ./mpia1_q2 25

-------------------------------------------------------------------------------------------------------------------

    For question 3, the program is named mpia1_q3.

    (10) Run on – q3file_16.txt for 4 processors, execute:
    [ user@chuck ~]$ mpirun -n 4 ./mpia1_qe q3file_16.txt

    (11) Run on – q3file_16.txt for 4 processors, execute:
    [ user@chuck ~]$ mpirun -n 9 ./mpia1_qe q3file_16.txt

    (12) Run on – q3file_16.txt for 4 processors, execute:
    [ user@chuck ~]$ mpirun -n 1 ./mpia1_qe q3file_16.txt

    (13) Run on – q3file_16.txt for 4 processors, execute:
    [ user@chuck ~]$ mpirun -n 11 ./mpia1_qe q3file_97.txt

    (14) Run on – q3file_16.txt for 4 processors, execute:
    [ user@chuck ~]$ mpirun -n 6 ./mpia1_qe q3file_97.txt

