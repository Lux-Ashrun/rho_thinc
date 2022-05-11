# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= /home/sbagley/git/BoxLibNew/BoxLib
FAST_HOME ?= /home/sbagley/git/fast

NDEBUG    :=
MPI       :=
OMP       :=
PROF      :=
COMP      :=gfortran
MKVERBOSE := t
TARGETHOST := t

INC_DIAGNOSTICS := t


EXE_NAME = rhothinc
include $(FAST_HOME)/Tools/F_mk/GMakedefs.mak


include ./src/GPackage.mak
VPATH_LOCATIONS += ./src

#main.$(suf).exe: $(objects)
#	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

$(EXE_NAME).$(suf).exe: $(objects)
	$(LINK.f90) -o $(EXE_NAME).$(suf).exe $(objects) $(libraries)

include $(FAST_HOME)/Tools/F_mk/GMakerules.mak
#DEFINES := -Dtest
#-Dtest=$(WENO_ORDER)
print-%  : ; @echo $* = $($*)

F90=gfortran
FC=gfortran
CC=gcc
CXX=g++




