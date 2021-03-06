
include ../watson.make

# Compiler and compiler flags
  ifndef CC
    CC = gcc
  endif
  ifndef OPT
    OPT  = -O2 -march=native -flto
  endif
  ifndef WARN
    WARN = -Wall -Wextra -Werror
  endif

# Linked libraries
  ifdef PACKAGES
    LIB += `pkg-config --libs $(PACKAGES)`
  endif

# includes
  # GNU scientific library
  ifeq ($(findstring gsl, $(PACKAGES)), gsl)
    GSLINC   = `pkg-config --cflags gsl`
  endif

# Exexutable Directory
  ifndef EXEDIR
    EXEDIR = $(if ${MyLocalPath}, ${MyLocalPath}, ../bin)
  endif
  ifndef EXE
    EXE = $(EXEDIR)/$(EXENAME)
  endif

# Resulting objects
  OBJ     = $(notdir     $(SRC:.c=.o))
  GSLOBJ  = $(notdir  $(GSLSRC:.c=.o))

  TOTSRC  = $(SRC) $(GSLSRC)
  TOTOBJ  = $(OBJ) $(GSLOBJ) gitversion.o

# set path list searched for targets and prerequisites
  VPATH = $(sort $(dir $(TOTSRC)))

# set program version
  ifndef PROGRAM_VERSION
    PROGRAM_VERSION = $(shell git describe --tags --always --broken)
  endif

.PHONY: all
all: Makefile make.src $(EXE)

.PHONY: Watson
Watson: all

# Create gitversion object file
#   * define a version string
#   * do nothing if this matches the contents of the gitversion.c file
#   * else recreate gitversion.c
#   * then build gitversion.o
VERSIONSTRING = const char *gitversion = "$(PROGRAM_VERSION)";
makeVersionSrc = $(shell echo '$(VERSIONSTRING)' | cmp -s -- gitversion.c || echo '$(VERSIONSTRING)' > gitversion.c)
gitversion.o: gitversion.c | $(makeVersionSrc)
	$(CC) -w $(PPF) -c $?

# Build object files out of C-source files
$(OBJ): %.o: %.c
	$(CC) $(OPT) $(WARN) $(PPF) $(INC) -c $?

# Build GNU GSL objects
$(GSLOBJ): %.o: %.c
	$(CC) $(OPT) $(WARN) $(PPF) $(INC) $(GSLINC) -c $?

# Link all objects to create the executable
$(EXE): $(TOTOBJ) | $(EXEDIR)
	$(CC) $(OPT) $(WARN) $(LIB) $(TOTOBJ) -o $@

# Create executable directory
$(EXEDIR):
	mkdir -p $(EXEDIR)

# Create "do nothing" recipe to forcefully build objects
# similar to .PHONY, but also works if target is a file and exists
.PHONY: FORCE
FORCE:


# Allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# Call ctags with corresponding exclude pattern
.PHONY: tags
tags: FORCE
	ctags --recurse --exclude="test/*" --exclude="stash/*"

# Remove all generated binary files
.PHONY: clean
clean:
	rm -f $(TOTOBJ) gitversion.[oc] $(EXE)
	rmdir -p --ignore-fail-on-non-empty $(EXEDIR)
