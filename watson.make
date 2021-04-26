# vim: set filetype=make :

# Executable name
# PROGRAM_VERSION = test
# EXEDIR  = bin
  EXENAME = TheWatson

# Additional libraries, includes and packages
  CC       = gcc
  WARN     = -Wall -Wextra -Werror
  LIB      = -lm
  INC      = -Iinclude/
  PACKAGES = gsl

# Preprocessor flags
# PPF += -Ddebug_normalmode
# PPF += -Ddebug_coords
# PPF += -Ddebug_moment_of_inertia
# PPF += -Ddebug_mu

# List of source files
  include make.src
