EXE	= swtest
CC	= mpicc
SCC	= sw5cc
HFLAGS	= -host   -O2 -msimd 
SFLAGS	= -slave  -O2 -msimd 


$(EXE):swtest.o mypbuf.o setupNetwork.o runNetwork_sw.o runNetwork_sl_mpi.o my_slave.o dmaWait.o
	$(CC)  *.o -o $@

swtest.o: ./swtest.c
	$(CC) -c $(HFLAGS) ./swtest.c

mypbuf.o: ./mypbuf.c
	$(CC) -c $(HFLAGS) ./mypbuf.c

setupNetwork.o: ./setupNetwork.c
	$(CC) -c $(HFLAGS) ./setupNetwork.c
	
runNetwork_sw.o: ./runNetwork_sw.c
	$(CC) -c $(HFLAGS) ./runNetwork_sw.c

runNetwork_sl_mpi.o: ./runNetwork_sl_mpi.c
	$(SCC) -c $(SFLAGS) ./runNetwork_sl_mpi.c

my_slave.o: ./my_slave.c
	$(SCC) -c $(SFLAGS) ./my_slave.c

dmaWait.o: ./dmaWait.c
	$(SCC) -c $(SFLAGS) ./dmaWait.c

clean:
	rm *.o
	rm $(EXE)
