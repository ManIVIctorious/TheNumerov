# Compiler
  CC   = gcc
  CppC = g++
# List of compiler flags
  CFLAGS += -g                      # Enable debug symbols
  CFLAGS += -Og                     # Set optimisation level, should be g if debug symbols are enabled
  CFLAGS += -march=native           # Tune for current chipset, don't bother about backwards compatibility
 #CFLAGS += -mtune=native           # Tune for current chipset, remain backwards compatible

  CFLAGS += -Werror                 # Treat warnings as errors
  CFLAGS += -Wall                   # Enable base set of warnings
  CFLAGS += -Wextra                 # Enable additional warnings
  CFLAGS += -Wno-sign-compare       # Disable sign-compare warning
  CFLAGS += -Wno-misleading-indentation
 #CFLAGS += -Wstrict-prototypes     # Enable strict-prototypes warning  | not allowed in C++
 #CFLAGS += -Wmissing-prototypes    # Enable missing-prototypes warning | not allowed in C++
  CFLAGS += -Wno-unused-but-set-variable
 #CFLAGS = -g -w                    # Disable all warnings

# Preprocessor flags (compile time flags)
  PPF += -D HAVE_OPT_SPLINE         # add splining ability, requires HAVE_GSL_INSTALLED
  PPF += -D HAVE_GSL_INSTALLED      # GNU Scientific Library support
 #PPF += -D HAVE_MKL_INSTALLED      # intel math kernel library support
  PPF += -D HAVE_ARMA_INSTALLED     # Armadillo C++ linear algebra library


# Resulting executable
  EXE = bin/numerov2d


# List of linked libraries
  LIB += -lm

# List of resulting object files
# Standard objects
    OBJ += numerov2d.o
    OBJ += InputFunction.o
    OBJ += InputFunctionDipole.o
    OBJ += InputCoriolisCoefficients.o
    OBJ += MetaGetStencil.o
    OBJ += FillStencil2D.o
    OBJ += Help.o
  # GSL objects
    ifeq ($(findstring HAVE_GSL_INSTALLED, $(PPF)), HAVE_GSL_INSTALLED)
      ifeq ($(findstring HAVE_OPT_SPLINE, $(PPF)), HAVE_OPT_SPLINE)
        OBJ += cubic_spline.o
        GSLOBJ += spline_interpolate.o
      endif
    endif
  # MKL objects
    ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
      MKLOBJ += EigensolverFEAST_MKL_2D.o
    endif
  # Armadillo objects (require C++)
    ifeq ($(findstring HAVE_ARMA_INSTALLED, $(PPF)), HAVE_ARMA_INSTALLED)
      ARMAOBJ += EigensolverArmadillo_2D.o
    endif


# Additional linked libraries, library paths and include directories
  # GNU Scientific Library
    ifeq ($(findstring HAVE_GSL_INSTALLED, $(PPF)), HAVE_GSL_INSTALLED)
    # additional libraries
      LIB += `pkg-config --cflags --libs gsl`
    endif

  # Intel MKL
    ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
      MKLPATH = /usr/local/intel
      MKLINCDIR += -I$(MKLPATH)/mkl/include
      MKLLIBDIR += -L$(MKLPATH)/mkl/lib/intel64
      MKLLIBDIR += -L$(MKLPATH)/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64_lin
     #MKLLIBDIR += -L/opt/intel/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin

    # additional libraries
      LIB += -lmkl_core
      LIB += -lmkl_intel_ilp64
      LIB += -lmkl_intel_thread
      LIB += -liomp5
      LIB += -lpthread
      LIB += -ldl
      LIB += -m64
      LIB += -DMKL_ILP64
    endif

  # Armadillo ARPACK
    ifeq ($(findstring HAVE_ARMA_INSTALLED, $(PPF)), HAVE_ARMA_INSTALLED)
      ARMAINCDIR += -I/home/mjs/libs/armadillo-6.400.3/include/
      ARMAINCDIR += -I/home/mjs/libs/xianyi-OpenBLAS-afedc8e/lapack/
      ARMAINCDIR += -I/home/mjs/libs/arpack-arma-master/include
      ARMALIBDIR += -L/home/mjs/libs/xianyi-OpenBLAS-afedc8e/
      ARMALIBDIR += -L/usr/lib64

    # additional libraries
      LIB += -lopenblas
      LIB += -lstdc++
    endif


all: $(EXE)
# Build object files out of C-source files
%.o : %.c Makefile
	$(CC) $(CFLAGS) $(PPF) -c $<

# Build Intel MKL objects
ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
$(MKLOBJ): $(subst .o,.c,$(MKLOBJ))
	$(CC) $(CFLAGS) $(PPF) $(MKLINCDIR) -c $<
endif

# Build Armadillo ARPACK objects (require C++)
ifeq ($(findstring HAVE_ARMA_INSTALLED, $(PPF)), HAVE_ARMA_INSTALLED)
$(ARMAOBJ): $(subst .o,.cpp,$(ARMAOBJ))
	$(CppC) $(CFLAGS) $(PPF) $(ARMAINCDIR) -c $<
endif


# link all objects to create the executable
$(EXE): $(OBJ) $(GSLOBJ) $(MKLOBJ) $(ARMAOBJ)
	$(CC) $(CFLAGS) $(LIB) $(MKLLIBDIR) $(ARMALIBDIR) $(OBJ) $(GSLOBJ) $(MKLOBJ) $(ARMAOBJ) -o $@


# allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# remove all generated binary files
clean:
	rm -f $(OBJ) $(GSLOBJ) $(MKLOBJ) $(ARMAOBJ) $(EXE)
