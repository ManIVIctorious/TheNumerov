# Compiler
CC = gcc -g -Wall -Werror
# Compilerflags
#CFLAGS = -g03 -Wall -Werror -Wstrict-prototypes -Wno-parentheses -Wmissing-prototypes -mtune=native
CFLAGS = -I/usr/local/intel/mkl/include/ -L/usr/local/intel/mkl/lib/intel64 -L/usr/local/intel/composer_xe_2015.3.187/compiler/lib/intel64 -DMKL_ILP64
#CFLAGS = -I/opt/intel/mkl/include -L/opt/intel/mkl/lib/intel64 -L/opt/intel/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin -DMKL_ILP64 -w
# Zu erstellende fertige Programme (mit relativem Pfad vorangestellt)
EXE = bin/numerov2d
# Objekts (.o bzw. .out), Libraries und Dependencies
# (wenn Abhängigkeit geändert -> make wird neu ausgeführt)
OBJ =          numerov2d.o stencils.o cubic_spline.o InputFunction.o InputFunctionDipole.o Help.o
DEP = Makefile numerov2d.c stencils.c cubic_spline.c InputFunction.c InputFunctionDipole.c Help.c
LIB = `pkg-config --cflags --libs gsl` -lmkl_core -lmkl_intel_ilp64 -lmkl_intel_thread -liomp5 -lpthread -ldl -m64
# Kommandoblock "all:" als Einsprungspunkt für make
all: $(EXE)
#
$(EXE): $(OBJ) $(DEP)
	$(CC) $(CFLAGS) $(OBJ) $(LIB) -o $@

clean:
	for i in $(OBJ) $(EXE); do if [ -e $$i ]; then rm $$i; fi; done
