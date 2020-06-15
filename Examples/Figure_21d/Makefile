# Simple Makefile for SEP_propagator
# DTS & FE

comp=gfortran
#comp=ifort
flags=
#-O1

main: SEP_propagator.o
	$(comp) $(flags) -o SEP_Propagator SEP_propagator.o

SEP_propagator.o: SEP_propagator.f90
	$(comp) $(flags) -c SEP_propagator.f90

clean:
	rm SEP_Propagator
	rm SEP_propagator.o

