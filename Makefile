# Compiler
FC = gfortran

# Flags for the compiler
FFLAGS = -Wall -fbacktrace -fcheck=all

LASTF = -llapack -lblas

# Name of the executable
EXEC = opt

# Object files
OBJS = energy.o constants.o gradient.o functions.o main.o

# Default target: Compile and link the program
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS) $(LASTF)

# Compile the module file
energy.o: energy.f90
	$(FC) $(FFLAGS) -c energy.f90

# Compile the module file
constants.o: constants.f90
	$(FC) $(FFLAGS) -c constants.f90

# Compile the module file
gradient.o: gradient.f90
	$(FC) $(FFLAGS) -c gradient.f90 $(LASTF)

# Compile the module file
functions.o: functions.f90
	$(FC) $(FFLAGS) -c functions.f90

# Compile the main program file
main.o: main.f90
	$(FC) $(FFLAGS) -c main.f90  

# Clean the directory (remove object files, module files, and executable)
clean:
	rm -f $(OBJS) *.mod $(EXEC)