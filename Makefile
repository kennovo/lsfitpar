# $Id: Makefile,v 1.4 2015/05/12 04:38:29 kenno Exp $
# This is the line that will most often be in need of customization; the
# libraries in question are at different locations and have different names
# depending on the distro, architecture, choice of libraries (ATLAS, ACML, MKL,
# OpenBLAS,...)
LDLIBS= -llapack -lblas -lm
#LDLIBS= -L/opt/acml/acml5.1.0/gfortran64/lib -lacml -L/opt/acml/amdlibm-3.1-lin64/lib/dynamic -lamdlibm

# Change compiler according to your liking
CC=gcc
DEBUG= -g -D DEBUG
#DEBUG= -g -D DEBUG -pedantic # sadly, gopt triggers too many pedantic warnings
#DEBUG= -O2 -g -D DEBUG

# Adding machine-specific optimizations is up to the user
CFLAGS= -O2 -ffast-math
#CFLAGS= -O2 -ffast-math $(DEBUG)

# Shouldn't be changed by user
LSFPOBJS= lsfitpar.o gopt.o

## Commented out because it's redundant; see http://mrbook.org/tutorials/make/
# default    : lsfitpar
all        : lsfitpar
lsfitpar   : $(LSFPOBJS)
clean      :
	rm $(LSFPOBJS)
debug      : CFLAGS=$(DEBUG)
debug      : all

## As long as you're compiling a native binary on a little-endian architecture
## such as x86, everything below this point can safely be ignored; see the
## comments in the source code for explanation.

# m32 and m64 respectively instruct the compiler to (cross-)compile for 32- and
# 64-bit architecture. Completely independently, -D FINT32 and -D FINT64
# specify whether the fortran integers of the blas and lapack libraries are 32-
# or 64-bit. The example below reflects the fact that gfortran currently uses
# 32-bit integers both on 32-bit and 64-bit architectures. If memory serves me
# well, this was not always the case, and is not true for all compilers. At any
# rate, it is not guaranteed to be so, as the Fortran 2003 standard only
# specifies a minimum size for integers, and most fortran compilers feature a
# flag to control this. Bottom line is that the integer size depends on how the
# blas and lapack libraries of choice were compiled...
gfortran32 : CFLAGS += -D FINT32
gfortran32 : m32

gfortran64 : CFLAGS += -D FINT32
gfortran64 : m64

m32        : CFLAGS += -m32
m32        : LDFLAGS += -m32
m32        : all

m64        : CFLAGS += -m64
m64        : LDFLAGS += -m64
m64        : all
