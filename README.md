# SEP_propagator
# For details see van den Berg et al., 2020: ArXiv link: xxx.xxx.xx

1D solar energetic particle (SEP) transport model

To compile the code use either the Intel Fortran compiler with

user@fskrdts-lap:~/Downloads/1D model/Test> ifort -O2 SEP_propagator.f90

and run the executable a.out or ./a.out depending on your system.

or

compile with the Gfortran compiler as

user@fskrdts-lap:~/Downloads/1D model/Test> gfortran -O2 SEP_propagator.f90

and run the executable a.out or ./a.out depending on your system.

Ater running the .f90 code we also provide a Python code that does some basic plotting of the output. It seems to work with

user@fskrdts-lap:~/Downloads/1D model/Test> python plot_output.py 

i.e. might not be Python3 compatable.
