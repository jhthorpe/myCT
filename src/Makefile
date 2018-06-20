myutil=/home/james.thorpe/bin/myUtils/lib
FC=gfortran
FFLAGS= -fcheck=bounds -g

objects = input.o levels_HO.o energy_HO.o energy_VPT2.o enumerate_HO.o

all: $(objects)
	$(FC) $(FFLAGS) ct.f90 -o myCT $(objects) $(myutil)/myUtils.a -I$(myutil)

clean:
	rm -f *.o myCT

%.o : %.f90
	$(FC) $(FFLAGS) -c $<
