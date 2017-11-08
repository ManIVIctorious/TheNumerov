# Compiler and compiler flags
CC = gcc
CFLAGS = -g -Og -Wall -Wextra -Wstrict-prototypes -Wmissing-prototypes -Werror -Wno-unused-parameter -march=native #-mtune=native

# additional header and library file directories
#INCDIR = -I/usr/local/intel/mkl/include/
#LIBDIR = -L/usr/local/intel/mkl/lib/intel64 -L/usr/local/intel/composer_xe_2015.3.187/compiler/lib/intel64
#INCDIR = -I/opt/intel/mkl/include
#LIBDIR = -L/opt/intel/mkl/lib/intel64 -L/opt/intel/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin

# executable to be generated
EXE1 = bin/watson_correction
EXE2 = bin/coriolis_coefficients
OBJ1 = main_coriolis_correction.o   InputNormalMode.o InputComFile.o InvertMatrix.o
OBJ2 = main_coriolis_coefficients.o InputNormalMode.o
#LIB = `pkg-config --cflags --libs gsl` -lmkl_core -lmkl_intel_ilp64 -lmkl_intel_thread -liomp5 -lpthread -ldl -m64 -DMKL_ILP64
LIB = `pkg-config --cflags --libs gsl`

all: $(EXE1) $(EXE2)
# define rule to build object files out of C-source files
%.o : %.c
	$(CC) $(CFLAGS) $(INCDIR) $(LIB) -c $<

# link all objects to create the executable
$(EXE1): $(OBJ1)
	$(CC) $(CFLAGS) $(LIBDIR) $(LIB) $(OBJ1) -o $@

$(EXE2): $(OBJ2)
	$(CC) $(CFLAGS) $(LIBDIR) $(LIB) $(OBJ2) -o $@

clean:
	rm -f $(OBJ1) $(EXE1) $(OBJ2) $(EXE2)
