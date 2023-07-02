CC = mpicc

CFLAGS = -g -Wall #-D_CC_OVERLAP

LDFLAGS = -lm

SRCS = jacobi.c decomp.c
OBJS = $(SRCS:.c=.o)
DEPS = jacobi.h decomp.h poisson1d.h
EXECS = poiss1dq12 poiss1dq3 poiss2d poiss2dnb

all: $(EXECS)

poiss1dq12: q12.o $(OBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
poiss1dq3: q3.o $(OBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
poiss2d: q4.o $(OBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
poiss2dnb: q4nb.o $(OBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

q12.o:q12.c $(DEPS)
	$(CC) -c $< -o $@ $(CFLAGS)
q3.o:q3.c $(DEPS)
	$(CC) -c $< -o $@ $(CFLAGS)
q4.o:q4.c $(DEPS)
	$(CC) -c $< -o $@ $(CFLAGS)
q4nb.o:q4nb.c $(DEPS)
	$(CC) -c $< -o $@ $(CFLAGS)


.PHONY: clean tags tests

clean:
	rm *.o *.txt $(EXECS)
