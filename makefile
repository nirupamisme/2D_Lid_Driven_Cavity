CC = g++

# Folders
BDIR = bin
IDIR = inc
LDIR = lib
ODIR = obj
SDIR = src

# Flags
CFLAGS = -I$(IDIR) -std=c++11 -O3

# For clang compiler
# CFLAGS = -I$(IDIR) -std=c++11 -O3

# Macros
_DEPS = LBM.h IOobject.h	# Add header files here
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))
_SRC = LBM.cpp IOobject.cpp main.cpp	# Add source files here
SRC = $(patsubst %, $(SDIR)/%, $(_SRC))
_OBJ = LBM.o IOobject.o main.o	# Add object files here
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -I.

$(BDIR)/lidDrivenCavity: $(OBJ)
	$(CC) -o $@ $^ $(CFALGS)

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

.PHONY: run
run:
	./lidDrivenCavity
