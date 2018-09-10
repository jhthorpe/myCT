myutil=/Users/jamesthorpe/local-hbar/myUtils/lib
FC=gfortran
#FC=ifort
FFLAGS= -fcheck=bounds -O3

objects = input.o energy_HO.o energy_VPT2.o hash_1Dint4.o recall_HO.o levels_HO.o init_HO.o stat_model.o

all: $(objects) $(myutil)
	$(FC) $(FFLAGS) ct.f90 -o myCT $(objects) $(myutil)/myUtils.a -I$(myutil)

clean:
	rm -f *.o myCT

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

levels_HO.o: levels_HO.f90 init_HO.o recall_HO.o 
	$(FC) $(FFLAGS) -c levels_HO.f90 init_HO.o recall_HO.o $(myutil)/myUtils.a -I$(myutil)

#recall_HO.o: recall_HO.f90 hash_1Dint4.o
#	$(FC) $(FFLAGS) -c recall_HO.f90 hash_1Dint4.o $(myutil)/myUtils.a -I$(myutil)

#recall_HO.o: recall_HO.f90
#	$(FC) $(FFLAGS) -c recall_HO.f90 $(myutil)/myUtils.a -I$(myutil)

recall_HO.o: recall_HO.f90 energy_HO.o hash_1Dint4.o
	$(FC) $(FFLAGS) -c recall_HO.f90 energy_HO.o hash_1Dint4.o $(myutil)/myUtils.a -I$(myutil)

stat_model.o: stat_model.f90
	$(FC) $(FFLAGS) -c stat_model.f90
