# Compiler
  CC   = gcc
  CppC = g++ -Wno-unused-parameter
# List of compiler flags
  CFLAGS += -g#                     # Enable debug symbols
  CFLAGS += -Og#                    # Set optimisation level, should be g if debug symbols are enabled
  CFLAGS += -march=native#          # Tune for current chipset, don't bother about backwards compatibility
 #CFLAGS += -mtune=native#          # Tune for current chipset, remain backwards compatible

  CFLAGS += -Wall#                  # Enable base set of warnings
  CFLAGS += -Wextra#                # Enable additional warnings
  CFLAGS += -Werror#                # Treat warnings as errors
 #CFLAGS = -g -w#                   # Disable all warnings

# Preprocessor flags (compile time flags)
  PPF += -D HAVE_GSL_INSTALLED#     # gnu scientific library support
  PPF += -D HAVE_MKL_INSTALLED#     # intel math kernel library support
  PPF += -D HAVE_ARMA_INSTALLED#    # Armadillo C++ linear algebra library


# Resulting executable
  EXEPATH = $(if ${MyLocalPath}, ${MyLocalPath}, bin)
  EXE = $(EXEPATH)/TheNumerov


# List of resulting object files
# Standard objects
    OBJ += main.o
    OBJ += GetSettingsControlFile.o
    OBJ += ParseControlFile.o
    OBJ += GetSettingsGetopt.o
    OBJ += CheckCoordinateSpacing.o
    OBJ += InputFunction.o
    OBJ += InputCoriolisCoefficients.o
    OBJ += TextOut.o
    OBJ += OutputSettings.o
    OBJ += MetaGetStencil.o
    OBJ += FillNumerovStencils.o
    OBJ += FillDerivativeStencils.o
    OBJ += MetaEigensolver.o
    OBJ += Help.o
    OBJ += MetaInterpolation.o
    OBJ += nx1dInterpolation.o
    OBJ += Integrators.o
  # GSL objects
    ifeq ($(findstring HAVE_GSL_INSTALLED, $(PPF)), HAVE_GSL_INSTALLED)
      GSLOBJ += GenFiniteDifferenceStencils.o
    endif
  # MKL objects
    ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
      MKLOBJ += Fill_MKL.o
      MKLOBJ += SolverFEAST_MKL.o
    endif
  # Armadillo objects (require C++)
    ifeq ($(findstring HAVE_ARMA_INSTALLED, $(PPF)), HAVE_ARMA_INSTALLED)
      ARMAOBJ += Fill_Armadillo.o
      ARMAOBJ += SolverARPACK_Armadillo.o
    endif


# List of linked libraries
  LIB += -lm

# Additional linked libraries, library paths and include directories
  # GNU scientific library
    ifeq ($(findstring HAVE_GSL_INSTALLED, $(PPF)), HAVE_GSL_INSTALLED)
      GSLINC = `pkg-config --cflags gsl`

    # additional libraries
      LIB += `pkg-config --libs gsl`
    endif

  # Armadillo ARPACK
    ifeq ($(findstring HAVE_ARMA_INSTALLED, $(PPF)), HAVE_ARMA_INSTALLED)
      ARMAINC = `pkg-config --cflags armadillo`

    # additional libraries
      LIB += `pkg-config --libs armadillo`
      LIB += -lstdc++
    endif

  # Intel MKL
    ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
      MKLINC += `pkg-config --cflags mkl`

    # additional libraries
      LIB += `pkg-config --libs mkl`
    endif


all: $(EXE)
gitversion.h:
	echo "static const char *gitversion = \"$(shell git describe --tags --always)\";" > $@

# Build object files out of C-source files
%.o : %.c Makefile gitversion.h
	$(CC) $(CFLAGS) $(PPF) -c $<

# Build GNU GSL objects
ifeq ($(findstring HAVE_GSL_INSTALLED, $(PPF)), HAVE_GSL_INSTALLED)
$(GSLOBJ): $(GSLOBJ:.o=.c)
	$(CC) $(CFLAGS) $(PPF) $(GSLINC) -c $?
endif

# Build Intel MKL objects
ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
$(MKLOBJ): $(MKLOBJ:.o=.c)
	$(CC) $(CFLAGS) $(PPF) $(MKLINC) -c $?
endif

# Build Armadillo ARPACK objects (require C++)
ifeq ($(findstring HAVE_ARMA_INSTALLED, $(PPF)), HAVE_ARMA_INSTALLED)
$(ARMAOBJ): $(ARMAOBJ:.o=.cpp)
	$(CppC) $(CFLAGS) $(PPF) $(ARMAINC) -c $?
endif


# link all objects to create the executable
$(EXE): $(OBJ) $(GSLOBJ) $(MKLOBJ) $(ARMAOBJ)
	$(CC) $(CFLAGS) $(LIB) $(OBJ) $(GSLOBJ) $(MKLOBJ) $(ARMAOBJ) -o $@


# allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# remove all generated binary files
clean:
	rm -f $(OBJ) $(GSLOBJ) $(MKLOBJ) $(ARMAOBJ) $(EXE) gitversion.h
