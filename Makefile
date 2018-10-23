# Compilers
  CC   = gcc
  CppC = g++ -Wno-unused-parameter

# List of compiler flags
  OPT  = -g -Og -march=native  # Set optimisation level to debug
  WARN = -Wall -Wextra -Werror # Enable base set and additional warnings and treat them as errors
  INC  = -I.

# Preprocessor flags (compile time flags)
  PPF += -D HAVE_GSL_INSTALLED#     # gnu scientific library support
  PPF += -D HAVE_MKL_INSTALLED#     # intel math kernel library support
  PPF += -D HAVE_ARMA_INSTALLED#    # Armadillo C++ linear algebra library

# Resulting executable
  EXEPATH = $(if ${MyLocalPath}, ${MyLocalPath}, bin)
  EXE = $(EXEPATH)/TheNumerov

# List of source files
# Standard source files
  SRC += main.c
  SRC += GetDefaultSettings.c
  SRC += GetSettingsControlFile.c
  SRC += ParseControlFile.c
  SRC += GetSettingsGetopt.c
  SRC += CheckCoordinateSpacing.c
  SRC += InputFunction.c
  SRC += InputCoriolisCoefficients.c
  SRC += TextOut.c
  SRC += OutputSettings.c
  SRC += MetaGetStencil.c
  SRC += FillNumerovStencils.c
  SRC += FillDerivativeStencils.c
  SRC += MetaEigensolver.c
  SRC += Help.c
  SRC += MetaInterpolation.c
  SRC += nx1dInterpolation.c
  SRC += Integrators.c
  # GSL sources
  ifeq ($(findstring HAVE_GSL_INSTALLED, $(PPF)), HAVE_GSL_INSTALLED)
    GSLSRC += GSL/GenFiniteDifferenceStencils.c
  endif
  # MKL sources
  ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
    MKLSRC += MKL/Fill_MKL_1D.c
    MKLSRC += MKL/Fill_MKL_2D.c
    MKLSRC += MKL/SolverFEAST_MKL.c
  endif
  # Armadillo sources (require C++)
  ifeq ($(findstring HAVE_ARMA_INSTALLED, $(PPF)), HAVE_ARMA_INSTALLED)
    ARMASRC += Armadillo/Fill_Armadillo_1D.cpp
    ARMASRC += Armadillo/Fill_Armadillo_2D.cpp
    ARMASRC += Armadillo/SolverARPACK_Armadillo.cpp
  endif

# Resulting objects
  OBJ     = $(notdir     $(SRC:.c=.o))
  GSLOBJ  = $(notdir  $(GSLSRC:.c=.o))
  MKLOBJ  = $(notdir  $(MKLSRC:.c=.o))
  ARMAOBJ = $(notdir $(ARMASRC:.cpp=.o))

# List of linked libraries
  LIB += -lm
  # GNU scientific library
  ifeq ($(findstring HAVE_GSL_INSTALLED, $(PPF)), HAVE_GSL_INSTALLED)
    GSLINC   = `pkg-config --cflags gsl`
    LIB     += `pkg-config --libs gsl`
  endif
  # Intel MKL
  ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
    MKLINC   = `pkg-config --cflags mkl`
    LIB     += `pkg-config --libs mkl`
  endif
  # Armadillo ARPACK
  ifeq ($(findstring HAVE_ARMA_INSTALLED, $(PPF)), HAVE_ARMA_INSTALLED)
    ARMAINC  = `pkg-config --cflags armadillo`
    LIB     += `pkg-config --libs armadillo` -lstdc++
  endif


all: Makefile gitversion.h $(EXE)
# Build gitversion header
gitversion.h:
	echo "static const char *gitversion = \"$(shell git describe --tags --always)\";" > $@

# Build object files out of C-source files
$(OBJ): $(SRC)
	$(CC) $(OPT) $(WARN) $(INC) $(PPF) -c $?

# Build GNU GSL objects
$(GSLOBJ): $(GSLSRC)
	$(CC) $(OPT) $(WARN) $(INC) $(PPF) $(GSLINC) -c $?

# Build Intel MKL objects
$(MKLOBJ): $(MKLSRC)
	$(CC) $(OPT) $(WARN) $(INC) $(PPF) $(MKLINC) -c $?

# Build Armadillo ARPACK objects (require C++)
$(ARMAOBJ): $(ARMASRC)
	$(CppC) $(OPT) $(WARN) $(INC) $(PPF) $(ARMAINC) -c $?


# link all objects to create the executable
$(EXE): $(OBJ) $(GSLOBJ) $(MKLOBJ) $(ARMAOBJ)
	$(CC) $(OPT) $(WARN) $(INC) $(LIB) $(OBJ) $(GSLOBJ) $(MKLOBJ) $(ARMAOBJ) -o $@


# allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# remove all generated binary files
clean:
	rm -f $(OBJ) $(GSLOBJ) $(MKLOBJ) $(ARMAOBJ) $(EXE) gitversion.h
