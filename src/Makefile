myutil=/home/james.thorpe/bin/myUtils/lib
FC=gfortran
#FC=ifort
FFLAGS= -fcheck=bounds -O3

objects = input.o energy_HO.o energy_VPT2.o hash_1Dint4.o enumerate_HO.o levels_HO.o init_HO.o

all: $(objects) $(needutil)
	$(FC) $(FFLAGS) ct.f90 -o myCT $(objects) $(myutil)/myUtils.a -I$(myutil)

clean:
	rm -f *.o myCT

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

levels_HO.o: levels_HO.f90 init_HO.o enumerate_HO.o 
	$(FC) $(FFLAGS) -c levels_HO.f90 init_HO.o enumerate_HO.o $(myutil)/myUtils.a -I$(myutil)

#enumerate_HO.o: enumerate_HO.f90 hash_1Dint4.o
#	$(FC) $(FFLAGS) -c enumerate_HO.f90 hash_1Dint4.o $(myutil)/myUtils.a -I$(myutil)

#enumerate_HO.o: enumerate_HO.f90
#	$(FC) $(FFLAGS) -c enumerate_HO.f90 $(myutil)/myUtils.a -I$(myutil)

enumerate_HO.o: enumerate_HO.f90 energy_HO.o hash_1Dint4.o
	$(FC) $(FFLAGS) -c enumerate_HO.f90 energy_HO.o hash_1Dint4.o $(myutil)/myUtils.a -I$(myutil)
