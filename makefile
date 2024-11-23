# Compiler
MPICXX = mpicxx

# Source file
SRC = spod_program.cpp
#SRC = read_data.cpp

# Output executable
TARGET = spod_program
#TARGET = read_data

# Environment variables
PETSC_DIR = /home1/03119/srmurth2/lib/petsc-3.21.2
SLEPC_DIR = /home1/03119/srmurth2/lib/slepc-3.21.1
PETSC_ARCH = arch-linux-c-debug

# Include directories
INCLUDES = -I$(PETSC_DIR)/include \
           -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
           -I$(SLEPC_DIR)/include \
           -I$(SLEPC_DIR)/$(PETSC_ARCH)/include \
           -I/opt/apps/intel24/impi21/fftw3/3.3.10/include \
           -I/opt/apps/intel24/impi21/phdf5/1.14.3/x86_64/include

# Library directories
LIB_DIRS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
           -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib \
           -L/opt/apps/intel24/impi21/fftw3/3.3.10/lib \
           -L/opt/apps/intel24/impi21/phdf5/1.14.3/x86_64/lib

# Libraries
LIBS = -lpetsc -lslepc -lfftw3_mpi -lfftw3 -lhdf5

# Compiler flags
CXXFLAGS = $(INCLUDES)

# Linker flags
LDFLAGS = $(LIB_DIRS) $(LIBS)

# Build rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(MPICXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET)

