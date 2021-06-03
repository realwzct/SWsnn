
rm -f *.o swtest

sw5cc -host -c swtest.c -O2 -msimd -I ./ -I/usr/sw-mpp/mpi2/include
sw5cc -host -c mypbuf.c -O2 -msimd -I ./ -I/usr/sw-mpp/mpi2/include
sw5cc -host -c setupNetwork.c -O0 -msimd -I ./ -I/usr/sw-mpp/mpi2/include
sw5cc -host -c runNetwork_sw.c -O2 -msimd -I ./ -I/usr/sw-mpp/mpi2/include

sw5cc -slave -c runNetwork_sl_mpi.c -O2 -msimd -I ./ -I/usr/sw-mpp/mpi2/include
sw5cc -slave -c my_slave.c -O2 -msimd -I ./ -I/usr/sw-mpp/mpi2/include
sw5cc -slave -c dmaWait.c -O2 -msimd -I ./ -I/usr/sw-mpp/mpi2/include

mpicc swtest.o mypbuf.o setupNetwork.o runNetwork_sw.o runNetwork_sl_mpi.o my_slave.o dmaWait.o -o swtest

