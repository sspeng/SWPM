
mpicc -host -O2 testphoton.c -lm -c
sw5cc -slave -msimd -O2 photon-slave.c -c
mpicc -hybrid -O2 testphoton.o photon-slave.o -lm -o pm
