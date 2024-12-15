# Compiler
FC = gfortran

# Detect the operating system
UNAME_S := $(shell uname -s)

# Use the correct lapack and blas implementation
# based on the OS

# For linux
ifeq ($(UNAME_S), Linux)
	ACCEL = -llapack -lblas
# For MacOS
else ifeq ($(UNAME_S), Darwin)
	ACCEL = -framework Accelerate
else
	ACCEL = 
endif

# Flags for the compiler
FFLAGS = -Wall -fbacktrace -fcheck=all

# Name of the executable
EXEC = opt

# Object files
OBJS = energy.o constants.o gradient.o functions.o main.o

# Default target: Compile and link the program
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS) $(ACCEL)

# Compile the module file
energy.o: energy.f90
	$(FC) $(FFLAGS) -c energy.f90

# Compile the module file
constants.o: constants.f90
	$(FC) $(FFLAGS) -c constants.f90

# Compile the module file
gradient.o: gradient.f90
	$(FC) $(FFLAGS) -c gradient.f90 $(ACCEL)

# Compile the module file
functions.o: functions.f90
	$(FC) $(FFLAGS) -c functions.f90

# Compile the main program file
main.o: main.f90
	$(FC) $(FFLAGS) -c main.f90  

# Clean the directory (remove object files, module files, and executable)
clean:
	rm -f $(OBJS) *.mod $(EXEC)