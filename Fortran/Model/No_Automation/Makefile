# This is how you comment in makefile
FC = ifort -parallel
FFLAGS = -c

# You can change the name of the executable
EXECUTABLE = No_Automation
F90_FILES := $(wildcard *.f90)
OBJECTS := $(patsubst %.f90, %.o, $(F90_FILES))
$(EXECUTABLE) : $(OBJECTS)
	$(FC) -o $@ $^
%.o: %.f90
	$(FC) $(FFLAGS) $< -o $@

# Below you need to list all the dependencies
main.o: Global_Variable.o
main.o: Simulation.o
Simulation.o: Global_Variable.o
Simulation.o: My_Function.o
clean :
	rm -f *.o *.mod $(EXECUTABLE)
