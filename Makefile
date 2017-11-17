# Compiler and compiler flags
CC = gcc
CFLAGS = -g -Og -Wall -Wextra -Wstrict-prototypes -Wmissing-prototypes -Werror -Wno-sign-compare -march=native #-mtune=native

# additional header and library file directories
INSTALLPATH = /usr/local/intel
INCDIR := $(INCDIR) -I$(INSTALLPATH)/mkl/include/
LIBDIR := $(LIBDIR) -L$(INSTALLPATH)/mkl/lib/intel64
LIBDIR := $(LIBDIR) -L$(INSTALLPATH)/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64_lin/
#LIBDIR := $(LIBDIR) -L$(INSTALLPATH)/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin

# executable to be generated
EXE = bin/numerov2d
OBJ = numerov2d.o stencils.o cubic_spline.o InputFunction.o InputFunctionDipole.o Help.o spline_interpolate.o
LIB = `pkg-config --cflags --libs gsl` -lmkl_core -lmkl_intel_ilp64 -lmkl_intel_thread -liomp5 -lpthread -ldl -m64 -DMKL_ILP64

all: $(EXE)
# define rule to build object files out of C-source files
%.o : %.c
	$(CC) $(CFLAGS) $(INCDIR) $(LIB) -c $<

# link all objects to create the executable
$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(LIBDIR) $(LIB) $(OBJ) -o $@

clean:
	rm -f $(OBJ) $(EXE)
