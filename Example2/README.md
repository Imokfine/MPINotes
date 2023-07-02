**To use MPI on chuck, run the code below**

```bash
module load cports openmpi
```

Run the Makefile

```bash
make
```



## Question 1

The functions were implemented in q12.c



## Question 2

**For 1D processor decomposition, when grid size = 15, run**

```bash
mpirun -np 4 ./poiss1dq12 15
```

The global difference is less than 1.0E-11, iterative solve converged

**For 1D processor decomposition, when grid size = 31, run**

```bash
mpirun -np 4 ./poiss1dq12 31
```

The global difference is less than 1.0E-11, iterative solve converged



## Question 3

The write_grid(), GatherGrid() and heatmap() functions were implemented in q3.c

**For 1D processor decomposition, when grid size = 31, to gather and write global solution to a file, run**

```bash
mpirun -np 4 ./poiss1dq3 31
```

The full solution will be write to the file "Grid_of_rank_0.txt"



## Question 4

**For 2D processor decomposition, which uses sendrecv(), run**

```bash
 mpirun -np 4 ./poiss2d 31
```

The global difference is less than 1.0E-11, iterative solve converged

The full solution will be write to the file "2d_Grid_of_rank_0.txt"

**For 2D processor decomposition, which uses non-blocking sends and receives, run**

```bash
 mpirun -np 4 ./poiss2dnb 31
```

The global difference is less than 1.0E-11, iterative solve converged

The full solution will be write to the file "2dnb_Grid_of_rank_0.txt"



**After comparison, the output of  "Grid_of_rank_0.txt", "2d_Grid_of_rank_0.txt" and "2dnb_Grid_of_rank_0.txt" are same.**

