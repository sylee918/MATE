#FC = mpif90 -fc=ifort
FC = ifort

flags = -O2 -mcmodel=medium
#flags = -mcmodel=medium -g -check all -traceback

rccx: rccx.x

rccx.x: Module_CIMI.o CXrate_from_CIMI.o sub.o geopack_2008.o TS04c.o 
	$(FC) -o rccx.x Module_CIMI.o CXrate_from_CIMI.o sub.o geopack_2008.o TS04c.o 


Module_CIMI.o: Module_CIMI.f90
	$(FC) -c $(flags) Module_CIMI.f90

CXrate_from_CIMI.o: CXrate_from_CIMI.f90
	$(FC) -c $(flags) CXrate_from_CIMI.f90

sub.o: sub.f90
	$(FC) -c $(flags) sub.f90

geopack_2008.o: geopack_2008.f
	$(FC) -c $(flags) geopack_2008.f

TS04c.o: TS04c.f
	$(FC) -c $(flags) TS04c.f


clean:
	rm *.o *.x *.mod