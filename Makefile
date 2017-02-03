# Compiler
CC = icc
# Compilerflags
#CFLAGS = -g03 -Wall -Werror -Wstrict-prototypes -Wno-parentheses -Wmissing-prototypes -mtune=native
CFLAGS = -I/usr/local/intel/composer_xe_2015.3.187/mkl/include/ -L/usr/local/intel/composer_xe_2015.3.187/compiler/lib/intel64 -L/usr/local/intel/composer_xe_2015.3.187/mkl/lib/intel64 -lgsl -lgslcblas -lmkl_core -lmkl_intel_ilp64 -lmkl_intel_thread -liomp5 -lpthread -lm -ldl -m64  -w -DMKL_ILP64
# Zu erstellende fertige Programme (mit relativem Pfad vorangestellt)
EXE = bin/numerov2d
# Objektdateien (.o bzw. .out), Librarys und Abhängigkeiten
# (wenn Abhängigkeit geändert -> make wird neu ausgeführt)
OBJ =          numerov2d.o stencils.o cubic_spline.o
DEP = Makefile numerov2d.c stencils.c cubic_spline.c
LIB  = -lm `pkg-config --cflags --libs gsl`
# Kommandoblock "all:" als Einsprungspunkt für make
all: $(EXE)
#
$(EXE): $(OBJ) $(DEP)
	$(CC) $(CFLAGS) $(OBJ) $(LIB) -o $@

clean:
	for i in $(OBJ) $(EXE); do if [ -e $$i ]; then rm $$i; fi; done
