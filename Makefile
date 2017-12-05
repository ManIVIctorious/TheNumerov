# Compiler
  CC = gcc
# List of compiler flags
  CFLAGS += -g                      # Enable debug symbols
  CFLAGS += -Og                     # Set optimisation level, should be g if debug symbols are enabled
  CFLAGS += -march=native           # Tune for current chipset, don't bother about backwards compatibility
 #CFLAGS += -mtune=native           # Tune for current chipset, remain backwards compatible

  CFLAGS += -Werror                 # Treat warnings as errors
  CFLAGS += -Wall                   # Enable base set of warnings
  CFLAGS += -Wextra                 # Enable additional warnings
  CFLAGS += -Wstrict-prototypes     # Enable strict-prototypes warning
  CFLAGS += -Wmissing-prototypes    # Enable missing-prototypes warning
  CFLAGS += -Wno-sign-compare       # Disable sign-compare warning
 #CFLAGS = -g -w                    # Disable all warnings
 #CFLAGS += -Wno-unused-parameter

# List of linked libraries
  LIB += `pkg-config --cflags --libs gsl`


# executables to be generated
  EXE1 = bin/coriolis_coefficients
  EXE2 = bin/watson_correction

  EXE = $(EXE1) $(EXE2)

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
	rm -f $(OBJ) $(EXE)
