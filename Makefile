# Compiler and compiler flags
CC = gcc
CFLAGS = -g -Og -Wall -Wextra -Wstrict-prototypes -Wmissing-prototypes -Werror -Wno-sign-compare -Wno-unused-parameter -march=native #-mtune=native

CFLAGS = -g -Og

# additional header and library file directories
#INCDIR = -I/usr/local/intel/mkl/include/
#LIBDIR = -L/usr/local/intel/mkl/lib/intel64 -L/usr/local/intel/composer_xe_2015.3.187/compiler/lib/intel64
#INCDIR = -I/opt/intel/mkl/include
#LIBDIR = -L/opt/intel/mkl/lib/intel64 -L/opt/intel/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin

# executable to be generated
EXE = bin/watson
OBJ = main.o InputComFile.o InputNormalMode.o InvertMatrix.o
#LIB = `pkg-config --cflags --libs gsl` -lmkl_core -lmkl_intel_ilp64 -lmkl_intel_thread -liomp5 -lpthread -ldl -m64 -DMKL_ILP64
LIB = `pkg-config --cflags --libs gsl`

all: $(EXE)
# define rule to build object files out of C-source files
%.o : %.c
	$(CC) $(CFLAGS) $(INCDIR) $(LIB) -c $<

# link all objects to create the executable
$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(LIBDIR) $(LIB) $(OBJ) -o $@

clean:
	rm -f $(OBJ) $(EXE)
