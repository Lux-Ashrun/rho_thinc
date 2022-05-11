# ============================================================================
# Name        : Makefile
# Author      : 
# Version     :
# Copyright   : Your copyright notice
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran

all: src/main.f90
	$(FORTRAN_COMPILER) -O2 -g \
		-o bin/main \
		src/main.f90

clean:
	rm -f bin/main *.mod
