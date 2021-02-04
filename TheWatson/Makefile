
include make.def

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

# Linked libraries and includes
  ifdef PACKAGES
    INC += `pkg-config --cflags $(PACKAGES)`
    LIB += `pkg-config --libs   $(PACKAGES)`
  endif

# Exexutable Directory
  ifndef EXEDIR
    EXEDIR = $(if ${MyLocalPath}, ${MyLocalPath}, ../bin)
  endif
  ifndef EXE
    EXE = $(EXEDIR)/$(EXENAME)
  endif

# Resulting objects
  OBJ = $(notdir $(SRC:.c=.o))

# set path list searched for targets and prerequisites
  VPATH = $(sort $(dir $(SRC)))

# set program version
  ifndef PROGRAM_VERSION
    PROGRAM_VERSION = $(shell git describe --tags --always)
  endif

.Phony: all
all: Makefile make.def $(EXE)

# Create gitversion source and object files
# First define a version string, if this string matches the contents
# of the gitversion.c file do nothing, else recreate this file
VERSIONSTRING = const char *gitversion = "$(PROGRAM_VERSION)";
gitversion.c: FORCE
	@echo '$(VERSIONSTRING)' | cmp -s -- $@ || echo '$(VERSIONSTRING)' > $@
gitversion.o: gitversion.c
	$(CC) -w $(PPF) -c $?

# Build object files out of C-source files
$(OBJ): %.o: %.c
	$(CC) $(OPT) $(WARN) $(INC) $(PPF) -c $?

# Link all objects to create the executable
$(EXE): $(OBJ) $(EXEDIR) gitversion.o
	$(CC) $(OPT) $(WARN) $(LIB) $(OBJ) gitversion.o -o $@

# Create executable directory
$(EXEDIR):
	mkdir -p $(EXEDIR)

# Create "do nothing" recipe to forcefully build objects
# similar to .Phony, but also works if target is a file and exists
.Phony: FORCE
FORCE:


# Allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# Call ctags with corresponding exclude pattern
.Phony: tags
tags: FORCE
	ctags --recurse --exclude="test/*" --exclude="stash/*"

# Remove all generated binary files
.Phony: clean
clean:
	rm -f $(OBJ) gitversion.[oc] $(EXE)
	rmdir -p --ignore-fail-on-non-empty $(EXEDIR)
