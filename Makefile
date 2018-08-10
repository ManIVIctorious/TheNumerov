# Compiler
  CC = gcc
# List of compiler flags
  CFLAGS = -O2 -Wall -Wextra -Werror -march=native

# Resulting executables
  EXEDIR = $(if ${MyLocalPath}, ${MyLocalPath}, bin)
  EXE1 = $(EXEDIR)/coriolis_coefficients
  EXE2 = $(EXEDIR)/watson_correction


# List of linked libraries
  LIB += `pkg-config --cflags --libs gsl`

# List of resulting object files
  OBJ1 += main_coriolis_coefficients.o
  OBJ1 += InputNormalMode.o

  OBJ2 += main_coriolis_correction.o
  OBJ2 += InputComFile.o
  OBJ2 += InputMasses.o
  OBJ2 += InputCoriolisCoefficients.o
  OBJ2 += InvertMatrix.o

  OBJ = $(OBJ1) $(OBJ2)

all: $(EXE1) $(EXE2)
# define rule to build object files out of C-source files
%.o : %.c
	$(CC) $(CFLAGS) $(INCDIR) $(LIB) -c $<

# link all objects to create the executable
$(EXE1): $(OBJ1) Makefile
	$(CC) $(CFLAGS) $(LIBDIR) $(LIB) $(OBJ1) -o $@

$(EXE2): $(OBJ2) Makefile
	$(CC) $(CFLAGS) $(LIBDIR) $(LIB) $(OBJ2) -o $@

clean:
	rm -f $(OBJ) $(EXE1) $(EXE2)
