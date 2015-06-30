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

# Warning: distributing a static binary that was created this way counts as
# redistributing a modified binary form of all the libraries to which we're
# linking. For example, if we have "LDLIBS= -llapack -lblas -lm" above, and
# we're linking against the libm that comes with gnu libc and BSD
# implementations of libblas and liblapack that were linked against
# libgfortran, we're also taking libm, libc, libgcc_s, libpthread and
# libquadmath with us. We need to satisfy the conditions for distributing all
# of these. Most of them are LGPL which can be linked into this GPL program
# without worries. The exception is the 3-clause BSD license of libblas and
# liblapack . I've seen different opinions on whether a BSD acknowledgement
# and/or disclaimer need to be distributed in the binary and/or its
# documentation. Most people don't, but to play it really safe, one should
# either consult a lawyer or just precautionary include the acknowledgements
# and disclaimers anyway.
# In modern versions of the Netlib libraries and derivatives (e.g. Atlas),
# xerbla_ is defined both in BLAS and LAPACK, leading to the dreaded "multiple
# definition" error. Of all the solutions that can be found online (further
# discussed in private README.static ), the following is by far the most
# portable (albeit also the most dangerous, because if there are other
# unexpected symbol clashes, it might quietly do the wrong thing).
# The -lgfortran is necessary because its functions are called by blas . It is
# possible to link lapack and blas static and libgfortran , libm , libpthread ,
# libquadmath ,... dynamic, but then we get a build procedure that is either
# vastly more complex or vastly less portable (we've been doing the latter in
# the past).
static     : LDLIBS += -lgfortran
static     : LDFLAGS += -static -Wl,--allow-multiple-definition
static     : all

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
