# This is how you comment in makefile
FC = ifort -parallel
FFLAGS = -I /N/soft/rhel7/nlopt/2.4.2/include -c
FFLAGS2 = -I /N/soft/rhel7/nlopt/2.4.2/include
FFLAGS3 = -L /N/soft/rhel7/nlopt/2.4.2/lib -lnlopt -lm 

# You can change the name of the executable
EXECUTABLE = Two_Group
F90_FILES := $(wildcard *.f90)
OBJECTS := $(patsubst %.f90, %.o, $(F90_FILES))
$(EXECUTABLE) : $(OBJECTS)
	$(FC) $(FFLAGS2) -o $@ $^ $(FFLAGS3)
%.o: %.f90
	$(FC) $(FFLAGS) $< -o $@

# Below you need to list all the dependencies
main.o: Global_Variable.o
main.o: Simulation.o
Simulation.o: Global_Variable.o
Simulation.o: My_Function.o
Simulation.o: Burkardt_fsolve.o
clean :
	rm -f *.o *.mod $(EXECUTABLE)
