OBJ  = ./build
EXE	= main
CC	= mpicc
SCC	= sw5cc.new
HFLAGS	= -host   -O2 -msimd 
SFLAGS	= -slave  -O2 -msimd
objects = $(OBJ)/swtest.o $(OBJ)/mypbuf.o $(OBJ)/setupNetwork.o $(OBJ)/runNetwork_sw.o $(OBJ)/runNetwork_sl_mpi.o $(OBJ)/my_slave.o $(OBJ)/dmaWait.o $(OBJ)/setup_slave.o


$(OBJ)/$(EXE):$(objects)
	$(CC)  $(objects) -o $(OBJ)/$(EXE)

$(OBJ)/swtest.o: swtest.c
	$(CC) -c $(HFLAGS) swtest.c -o $(OBJ)/swtest.o

$(OBJ)/mypbuf.o: mypbuf.c
	$(CC) -c $(HFLAGS) mypbuf.c -o $(OBJ)/mypbuf.o
$(OBJ)/setupNetwork.o: setupNetwork.c
	$(CC) -c $(HFLAGS) setupNetwork.c -o $(OBJ)/setupNetwork.o
	
$(OBJ)/runNetwork_sw.o: runNetwork_sw.c 
	$(CC) -c $(HFLAGS) runNetwork_sw.c -o $(OBJ)/runNetwork_sw.o

$(OBJ)/runNetwork_sl_mpi.o: runNetwork_sl_mpi.c 
	$(SCC) -c $(SFLAGS) runNetwork_sl_mpi.c -o $(OBJ)/runNetwork_sl_mpi.o

$(OBJ)/my_slave.o: my_slave.c
	$(SCC) -c $(SFLAGS) my_slave.c -o $(OBJ)/my_slave.o

$(OBJ)/dmaWait.o: dmaWait.c
	$(SCC) -c $(SFLAGS) dmaWait.c -o $(OBJ)/dmaWait.o

$(OBJ)/setup_slave.o: setup_slave.c
	$(SCC) -c $(SFLAGS) setup_slave.c -o $(OBJ)/setup_slave.o
run: 
	bsub -b -I -q q_sw_expr  -n 2 -np 2 -cgsp 64 -host_stack 256 -share_size 4096 $(OBJ)/main 1000 20000 2 0.125 1000 5.5  
clean:
	rm $(OBJ)/*.o
	rm $(OBJ)/$(EXE)
