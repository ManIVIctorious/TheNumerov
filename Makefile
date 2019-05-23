
include make.def

# Compiler and compiler flags
  ifndef CC
    CC = gcc
  endif
  ifndef OPT
    OPT  = -O2 -march=native
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
    EXEDIR = $(if ${MyLocalPath}, ${MyLocalPath}, bin)
  endif
  ifndef EXE
    EXE = $(EXEDIR)/$(EXENAME)
  endif

# Resulting objects
  OBJ = $(notdir $(SRC:.c=.o))


.Phony: all
all: $(EXE) Makefile make.def
# Build object files out of C-source files
$(OBJ): %.o : %.c
	$(CC) $(OPT) $(WARN) $(INC) $(PPF) -c $?

# Link all objects to create the executable
$(EXE): $(OBJ) $(EXEDIR)
	$(CC) $(OPT) $(WARN) $(INC) $(LIB) $(OBJ) -o $@

# Create executable directory
$(EXEDIR):
	mkdir -p $(EXEDIR)

# Allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# Remove all generated binary files
clean:
	rm -f $(OBJ) $(EXE)
	rmdir -p $(EXEDIR)
