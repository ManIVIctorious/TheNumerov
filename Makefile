
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
  ifndef EXE1
    EXE1 = $(EXEDIR)/$(EXENAME1)
  endif
  ifndef EXE2
    EXE2 = $(EXEDIR)/$(EXENAME2)
  endif

# Resulting objects
  OBJ1 = $(notdir $(SRC1:.c=.o))
  OBJ2 = $(notdir $(SRC2:.c=.o))

  OBJ = $(OBJ1) $(OBJ2)

.Phony: all
all: $(EXE1) $(EXE2) Makefile make.def
# Build object files out of C-source files
$(OBJ): %.o : %.c
	$(CC) $(OPT) $(WARN) $(INC) $(PPF) -c $?

# Link all objects to create the executable
$(EXE1): $(OBJ1) $(EXEDIR)
	$(CC) $(OPT) $(WARN) $(INC) $(LIB) $(OBJ1) -o $@
$(EXE2): $(OBJ2) $(EXEDIR)
	$(CC) $(OPT) $(WARN) $(INC) $(LIB) $(OBJ2) -o $@

# Create executable directory
$(EXEDIR):
	mkdir -p $(EXEDIR)

# Allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# Remove all generated binary files
clean:
	rm -f $(OBJ) $(EXE1) $(EXE2)
	rmdir -p $(EXEDIR)
