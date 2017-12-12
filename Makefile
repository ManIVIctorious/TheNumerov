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
  PPF += -D HAVE_OPT_SPLINE         # add splining ability, requires HAVE_GSL_INSTALLED
  PPF += -D HAVE_GSL_INSTALLED      # GNU Scientific Library support
  PPF += -D HAVE_MKL_INSTALLED      # intel math kernel library support


# Resulting executable
  EXE = bin/numerov2d


# List of linked libraries
  LIB += -lm

# additional header and library file directories for intel math kernel library
ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
  MKLPATH    = /usr/local/intel
  MKLINCDIR += -I$(MKLPATH)/mkl/include
  MKLLIBDIR += -L$(MKLPATH)/mkl/lib/intel64
  MKLLIBDIR += -L$(MKLPATH)/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64_lin
 #MKLLIBDIR += -L$(MKLPATH)/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin

# additional libraries for intel math kernel library
  MKLLIB += -lmkl_core
  MKLLIB += -lmkl_intel_ilp64
  MKLLIB += -lmkl_intel_thread
  MKLLIB += -liomp5
  MKLLIB += -lpthread
  MKLLIB += -ldl
  MKLLIB += -m64
  MKLLIB += -DMKL_ILP64
endif

# GNU Scientific Library
ifeq ($(findstring HAVE_GSL_INSTALLED, $(PPF)), HAVE_GSL_INSTALLED)
  GSLLIB += `pkg-config --cflags --libs gsl`
endif

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
  GSLOBJ += spline_interpolate.o
endif
ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
  MKLOBJ += EigensolverFEAST_MKL_2D.o
endif

all: $(EXE)
# define rule to build object files out of C-source files
%.o : %.c Makefile
	$(CC) $(CFLAGS) $(PPF) $(LIB) -c $<

# define rule to build object files requiring the Intel MKL
ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
$(MKLOBJ): $(subst .o,.c,$(MKLOBJ))
	$(CC) $(CFLAGS) $(PPF) $(LIB) $(MKLINCDIR) $(MKLLIB) -c $<
endif

# define rule to build object files requiring GNU Scientific Library
ifeq ($(findstring HAVE_GSL_INSTALLED, $(PPF)), HAVE_GSL_INSTALLED)
$(GSLOBJ): $(subst .o,.c,$(GSLOBJ))
	$(CC) $(CFLAGS) $(PPF) $(LIB) $(GSLLIB) -c $<
endif

# link all objects to create the executable
$(EXE): $(OBJ) $(GSLOBJ) $(MKLOBJ)
	$(CC) $(CFLAGS) $(LIB) $(GSLLIB) $(MKLLIBDIR) $(MKLLIB) $(OBJ) $(GSLOBJ) $(MKLOBJ) -o $@

# allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# remove all generated binary files
clean:
	rm -f $(OBJ) $(GSLOBJ) $(MKLOBJ) $(EXE)
