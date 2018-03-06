# Compiler
  CC = gcc
# List of compiler flags
  CFLAGS += -g                      # Enable debug symbols
  CFLAGS += -Og                     # Set optimisation level, should be g if debug symbols are enabled
  CFLAGS += -march=native           # Tune for current chipset, don't bother about backwards compatibility
 #CFLAGS += -mtune=native           # Tune for current chipset, remain backwards compatible

 #CFLAGS += -Werror                 # Treat warnings as errors
  CFLAGS += -Wall                   # Enable base set of warnings
  CFLAGS += -Wextra                 # Enable additional warnings
  CFLAGS += -Wstrict-prototypes     # Enable strict-prototypes warning
  CFLAGS += -Wmissing-prototypes    # Enable missing-prototypes warning
  CFLAGS += -Wno-sign-compare       # Disable sign-compare warning
 #CFLAGS = -g -w                    # Disable all warnings





# additional header and library file directories
  INSTALLPATH = /usr/local/intel


    # MKL includes
      MKLINC += -I$(INSTALLPATH)/mkl/include
      MKLINC += -DMKL_ILP64
      MKLINC += -m64

    # additional libraries
      MKLLIBDIR += -L$(MKLPATH)/mkl/lib/intel64
      MKLLIBDIR += -L$(MKLPATH)/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin
      #MKLLIBDIR += -L$(MKLPATH)/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin
      LIB += -lmkl_core
      LIB += -lmkl_intel_ilp64
      LIB += -lmkl_intel_thread
      LIB += -liomp5
      LIB += -lpthread
      LIB += -ldl 


##oldfiles


#  INCDIR += -I/usr/local/intel/composer_xe_2015.3.187/mkl/include/
#  LIBDIR += -L/usr/local/intel/composer_xe_2015.3.187/compiler/lib/intel64
#  LIBDIR += -L/usr/local/intel/composer_xe_2015.3.187/mkl/lib/intel64
 #LIBDIR += -L$(INSTALLPATH)/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin

# List of linked libraries
  LIB += `pkg-config --cflags --libs gsl`
#  LIB += -lmkl_core
#  LIB += -lmkl_intel_ilp64
#  LIB += -lmkl_intel_ilp64
#  LIB += -lmkl_intel_thread
#  LIB += -liomp5
#  LIB += -lpthread
  LIB += -fPIC
  LIB += -m64
  LIB += -DMKL_ILP64


# Resulting executable
  EXE = bin/numerov2d

# List of resulting object files
  OBJ += numerov2d.o
  OBJ += InputFunction.o
  OBJ += InputFunctionDipole.o
  OBJ += InputCoriolisCoefficients.o
  OBJ += stencils.o
  OBJ += cubic_spline.o
  OBJ += spline_interpolate.o
  OBJ += Help.o

all: $(EXE)
# define rule to build object files out of C-source files
%.o : %.c
	$(CC) $(CFLAGS) $(INCDIR) $(LIB) -c $<

# link all objects to create the executable
$(EXE): $(OBJ) Makefile
	$(CC) $(CFLAGS) $(LIBDIR) $(LIB) $(OBJ) -o $@

clean:
	rm -f $(OBJ) $(EXE)
