#FC = ifort
FC = mpif90 -fc=ifort

#flags = -mcmodel=medium -g -check all -warn external -warn declarations,interfaces -traceback
#flags = -fopenmp -mcmodel=medium -g -check all -traceback
#flags = -O2 -qopenmp -mcmodel=medium
flags = -O2 -mcmodel=medium
#flags = -O2 -axCORE-AVX512,CORE-AVX2 -xAVX
#flags = -mcmodel=medium -check all -traceback -g
#flags = -mcmodel=medium -check all -warn all,nodec,interfaces -gen_interfaces -traceback -fpe0 -ftrapuv # -fpstkchk
#flags = -mcmodel=medium -xCOMMON-AVX512
#flags = -O2 -xAVX

all: MATE.x

MATE.x: main1_trace.o init.o trace.o PSD.o extra_tools.o IO_utils.o
	$(FC) -o MATE.x main1_trace.o init.o trace.o PSD.o extra_tools.o IO_utils.o

main1_trace.o: main1_trace.f90
	$(FC) -c $(flags) main1_trace.f90
trace.o: trace.f90
	$(FC) -c $(flags) trace.f90
init.o: init.f90
	$(FC) -c $(flags) init.f90
extra_tools.o: extra_tools.f90
	$(FC) -c $(flags) extra_tools.f90
IO_utils.o: IO_utils.f90
	$(FC) -c $(flags) IO_utils.f90
PSD.o: PSD.f90
	$(FC) -c $(flags) PSD.f90


clean:
	rm *.o *.x
clean2:
	rm *.o *.x *.data *.dat
clean3:
	rm *.o *.x *.data *.dat *__genmod*

