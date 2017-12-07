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

# Preprocessor flags (compile time flags)
  PPF += -D HAVE_OPT_SPLINE         # add splining ability
  PPF += -D HAVE_MKL_INSTALLED      # intel math kernel library support


# Resulting executable
  EXE = bin/numerov2d


# additional header and library file directories for intel math kernel library
ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
  INSTALLPATH = /usr/local/intel
  INCDIR += -I$(INSTALLPATH)/mkl/include
  LIBDIR += -L$(INSTALLPATH)/mkl/lib/intel64
  LIBDIR += -L$(INSTALLPATH)/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64_lin
 #LIBDIR += -L$(INSTALLPATH)/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin

# additional libraries for intel math kernel library
  LIB += -lmkl_core
  LIB += -lmkl_intel_ilp64
  LIB += -lmkl_intel_ilp64
  LIB += -lmkl_intel_thread
  LIB += -liomp5
  LIB += -lpthread
  LIB += -ldl
  LIB += -m64
  LIB += -DMKL_ILP64
endif

# List of linked libraries
  LIB += `pkg-config --cflags --libs gsl`

# List of resulting object files
  OBJ += numerov2d.o
  OBJ += InputFunction.o
  OBJ += InputFunctionDipole.o
  OBJ += InputCoriolisCoefficients.o
  OBJ += MetaGetStencil.o
  OBJ += FillStencil2D.o
  OBJ += Help.o
ifeq ($(findstring HAVE_OPT_SPLINE, $(PPF)), HAVE_OPT_SPLINE)
  OBJ += cubic_spline.o
  OBJ += spline_interpolate.o
endif

all: $(EXE)
# define rule to build object files out of C-source files
%.o : %.c Makefile
	$(CC) $(CFLAGS) $(PPF) $(INCDIR) $(LIB) -c $<

# link all objects to create the executable
$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(LIBDIR) $(LIB) $(OBJ) -o $@

# allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# remove all generated binary files
clean:
	rm -f $(OBJ) $(EXE)
