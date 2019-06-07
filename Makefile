
include make.def

# Compiler and compiler flags
    CC = gcc
  ifndef CppC
    CppC = g++ -Wno-unused-parameter
  endif
  ifndef OPT
    OPT  = -O2 -march=native
  endif
  ifndef WARN
    WARN  = -Wall -Wextra -Werror
    WARN += -Wno-stringop-truncation -Wno-stringop-overflow
  endif

# Linked libraries
  ifdef PACKAGES
    LIB += `pkg-config --libs $(PACKAGES)`
  endif

# includes
  # GNU scientific library
  ifeq ($(findstring HAVE_GSL_INSTALLED, $(PPF)), HAVE_GSL_INSTALLED)
    GSLINC   = `pkg-config --cflags gsl`
  endif
  # Intel MKL
  ifeq ($(findstring HAVE_MKL_INSTALLED, $(PPF)), HAVE_MKL_INSTALLED)
    MKLINC   = `pkg-config --cflags mkl`
  endif
  # Armadillo ARPACK
  ifeq ($(findstring HAVE_ARMA_INSTALLED, $(PPF)), HAVE_ARMA_INSTALLED)
    ARMAINC  = `pkg-config --cflags armadillo`
    LIB += -lstdc++
  endif

# Exexutable Directory
  ifndef EXEDIR
    EXEDIR = $(if ${MyLocalPath}, ${MyLocalPath}, bin)
  endif
  ifndef EXE
    EXE = $(EXEDIR)/$(EXENAME)
  endif

# Resulting objects
  OBJ     = $(notdir     $(SRC:.c=.o))
  GSLOBJ  = $(notdir  $(GSLSRC:.c=.o))
  MKLOBJ  = $(notdir  $(MKLSRC:.c=.o))
  ARMAOBJ = $(notdir $(ARMASRC:.cpp=.o))

  TOTSRC  = $(SRC) $(GSLSRC) $(MKLSRC) $(ARMASRC)
  TOTOBJ  = $(OBJ) $(GSLOBJ) $(MKLOBJ) $(ARMAOBJ)

# set path list searched for targets and prerequisites
  VPATH = $(sort $(dir $(TOTSRC)))

# version files
# GITHEAD changes with every branch switch
# GITREF  changes with every new commit
  GITHEAD = .git/HEAD
  GITREF  = $(shell sed "s+^\s*ref:\s*+.git/+" .git/HEAD)
  ifndef PROGRAM_VERSION
    PROGRAM_VERSION = $(shell git describe --tags --always)
  endif

.Phony: all
all: Makefile make.def make.src $(EXE)
# Create gitversion source and object file:
gitversion.c: $(GITHEAD) $(GITREF)
	echo "const char *gitversion = \"$(PROGRAM_VERSION)\";" > $@
gitversion.o: gitversion.c
	$(CC) -w $(PPF) -c $?

# Build object files out of C-source files
$(OBJ): %.o: %.c
	$(CC) $(OPT) $(WARN) $(INC) $(PPF) -c $?

# Build GNU GSL objects
$(GSLOBJ): %.o: %.c
	$(CC) $(OPT) $(WARN) $(INC) $(PPF) $(GSLINC) -c $?

# Build Intel MKL objects
$(MKLOBJ): %.o: %.c
	$(CC) $(OPT) $(WARN) $(INC) $(PPF) $(MKLINC) -c $?

# Build Armadillo ARPACK objects (require C++)
$(ARMAOBJ): %.o: %.cpp
	$(CppC) $(OPT) $(WARN) $(INC) $(PPF) $(ARMAINC) -c $?

# Link all objects to create the executable
$(EXE): $(TOTOBJ) $(EXEDIR) gitversion.o
	$(CC) $(OPT) $(WARN) $(LIB) $(TOTOBJ) gitversion.o -o $@

# Create executable directory
$(EXEDIR):
	mkdir -p $(EXEDIR)

# Allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# Remove all generated binary files
clean:
	rm -f $(TOTOBJ) gitversion.o gitversion.c $(EXE)
	rmdir -p $(EXEDIR) 2>/dev/null
