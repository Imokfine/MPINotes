CC = mpicc

CFLAGS = -g -Wall #-D_CC_OVERLAP

LDFLAGS = -lm

SRCS = jacobiRMA.c decomp.c
OBJS = $(SRCS:.c=.o)
DEPS = jacobiRMA.h decomp.h poisson.h
EXECS = poiss2dPSCW poiss2dRMA 

all: $(EXECS)

poiss2dPSCW:poisson2dPSCW.o jacobiRMA.o decomp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

poisson2dPSCW.o:poisson2dPSCW.c jacobiRMA.h decomp.h poisson.h
	$(CC) -c $< -o $@ $(CFLAGS)

poiss2dRMA:poisson2dRMA.o jacobiRMA.o decomp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

poisson2dRMA.o:poisson2dRMA.c jacobiRMA.h decomp.h poisson.h
	$(CC) -c $< -o $@ $(CFLAGS)

jacobiRMA.o:jacobiRMA.c jacobiRMA.h poisson.h
	$(CC) -c $< -o $@ $(CFLAGS)

decomp.o:decomp.c
	$(CC) -c $< -o $@ $(CFLAGS)

.PHONY: clean tags tests

clean:
	rm *.o *.txt $(EXECS)
