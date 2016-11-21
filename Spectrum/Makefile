SHELL=/bin/bash

CC       = mpicc
OPTIMIZE = -O2 -Wall -std=c99
MPICHLIB = -lmpich
MPICHINCL=

#---------target system
SYSTEM="Linux"
#SYSTEM="Cluster"

ifeq ($(SYSTEM), "Linux")
NCORE      :=$(grep -c ^processor /proc/cpuinfo)
GSL_INCL    = $(shell pkg-config --cflags gsl)
GSL_LIBS    = $(shell pkg-config --libs gsl)
LAPACK_INCL = -I/usr/include/
LAPACK_LIBS = -L/usr/lib/lapack -llapacke -lblas
endif

ifeq ($(SYSTEM), "Cluster")
GSL_INCL = -I/mbh/mbhd01/user/liyanrong/soft/gsl/include
GSL_LIBS = -L/mbh/mbhd01/user/liyanrong/soft/gsl/lib  -lgsl -lgslcblas -lm
MPICHLIB = -L/mbh/mbhd01/user/liyanrong/soft/mpich3/lib -lmpich
MPIINCL  = -I/mbh/mbhd01/user/liyanrong/soft/mpich3/include
LAPACK_INCL = -I/mbh/mbhd01/user/liyanrong/soft/lapack/include
LAPACK_LIBS = -L/mbh/mbhd01/user/liyanrong/soft/lapack/lib -llapacke -llapack -lblas -lgfortran
#CBLAS_INCL  = -I/mbh/mbhd01/user/liyanrong/soft/atlas/include
#CBLAS_LIBS  = -L/mbh/mbhd01/user/liyanrong/soft/atlas/lib -lcblas
endif

EXEC     = Radiation
SRC      = ./src
OBJS     = $(SRC)/main.o $(SRC)/mathfun.o
INCL     = Makefile $(SRC)/main.h  $(SRC)/mathfun.h
OPTIONS  = $(OPTIMIZE)

ifeq ($(SYSTEM), "Linux")
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(CBLAS_INCL)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(CBLAS_LIBS)
endif

ifeq ($(SYSTEM), "Cluster")
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(CBLAS_INCL) $(MPIINCL)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(CBLAS_LIBS) $(MPICHLIB)
endif

$(EXEC):$(OBJS)
	cd $(SRC)
	$(CC) $(OBJS) $(LIBS) -o $@

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC) $(SRC)/*~ *~
