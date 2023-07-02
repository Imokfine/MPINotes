**To use MPI on chuck, run the code below**

```bash
module load cports openmpi
```

Run the Makefile

```bash
make
```



## Question

**For Fence Version, run**

```bash
mpirun -np 4 ./poiss2dRMA
```



**For PSCW Version, run**

```bash
mpirun -np 4 ./poiss2dPSCW
```



**After comparison with the previous assignment, the output of  "2d_Grid_of_rank_0.txt", "RMA2d_Grid_of_rank_0.txt" and "PSCW2d_Grid_of_rank_0.txt" are same. Thus, RMA operations match the answers of my non-RMA version**

