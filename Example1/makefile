#
# @file Makefile
# @Makefile for Assignment 1 for HPC software MPI
# @author Xueni Fu
# @version 1.0
# @date 2022-12-23
#

all: mpia1_q1_sr mpia1_q1_sg mpia1_q2 mpia1_q3 create_vector

mpia1_q1_sr:mpia1_q1_sr.c
	mpicc mpia1_q1_sr.c -o mpia1_q1_sr

mpia1_q1_sg:mpia1_q1_sg.c
	mpicc mpia1_q1_sg.c -o mpia1_q1_sg
 	
mpia1_q2:mpia1_q2.c
	mpicc mpia1_q2.c -o mpia1_q2

mpia1_q3:mpia1_q3.c
	mpicc mpia1_q3.c -o mpia1_q3

create_vector:create_vector.c
	gcc create_vector.c -o create_vector

.Phony: 

clean:
	rm -f mpia1_q1_sr mpia1_q1_sg mpia1_q2 mpia1_q3 create_vector


