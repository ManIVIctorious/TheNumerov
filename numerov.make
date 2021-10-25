# vim: set filetype=make :

# Executable name
# PROGRAM_VERSION = test
# EXEDIR  = bin
  EXENAME = TheNumerov

# Additional settings, including warnings, optimisation
#   level, libraries, includes and packages (pkg-config files)
  CC       = gcc
  WARN     = -Wall -Wextra -Werror -Wno-unused-parameter
  LIB      = -lm
  INC      = -Iinclude/
  PACKAGES = gsl mkl armadillo

# Include preprocessor flags based on used packages
  ifeq ($(findstring gsl, $(PACKAGES)), gsl)
    PPF += -D HAVE_GSL_INSTALLED  # gnu scientific library support
  endif
  ifeq ($(findstring mkl, $(PACKAGES)), mkl)
    PPF += -D HAVE_MKL_INSTALLED  # intel math kernel library support
  endif
  ifeq ($(findstring armadillo, $(PACKAGES)), armadillo)
    PPF += -D HAVE_ARMA_INSTALLED # Armadillo C++ linear algebra library
  endif

# List of source files
  include make.src
