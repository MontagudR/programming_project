# Compiler
FC = gfortran

# Flags for the compiler
FFLAGS = -Wall -fbacktrace -fcheck=all

LASTF = -llapack -lblas

# Name of the executable
EXEC = opt

# Object files
OBJS = parse_coordinates.o constants.o cartesian_gradient.o internal_gradient.o main.o

# Default target: Compile and link the program
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS) $(LASTF)

# Compile the module file
parse_coordinates.o: parse_coordinates.f90
	$(FC) $(FFLAGS) -c parse_coordinates.f90

# Compile the module file
constants.o: constants.f90
	$(FC) $(FFLAGS) -c constants.f90

# Compile the module file
cartesian_gradient.o: cartesian_gradient.f90
	$(FC) $(FFLAGS) -c cartesian_gradient.f90

# Compile the module file
internal_gradient.o: internal_gradient.f90
	$(FC) $(FFLAGS) -c internal_gradient.f90 $(LASTF)

# Compile the main program file
main.o: main.f90
	$(FC) $(FFLAGS) -c main.f90  

# Clean the directory (remove object files, module files, and executable)
clean:
	rm -f $(OBJS) *.mod $(EXEC)