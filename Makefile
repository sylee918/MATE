#FC = ifort
FC = mpif90 -fc=ifort

flags = -O2 -mcmodel=medium
#flags = -mcmodel=medium -g -check all -warn external -warn declarations,interfaces -traceback
#flags = -O2 -qopenmp -mcmodel=medium
#flags = -fopenmp -mcmodel=medium -g -check all -traceback
#flags = -O2 -axCORE-AVX512,CORE-AVX2 -xAVX
#flags = -mcmodel=medium -g -check all -traceback
#flags = -mcmodel=medium -check all -warn all,nodec,interfaces -gen_interfaces -traceback -fpe0 -ftrapuv # -fpstkchk
#flags = -mcmodel=medium -xCOMMON-AVX512
#flags = -O2 -xAVX

all: MATE.x

esc: MATE_esc.x

MATE.x: Module_Setting.o Module_MATE.o Module_BC.o Module_Lyman.o Module_Vol.o Module_ChargeExchange.o main.o init.o trace.o PSD.o IO_utils.o
	$(FC) -o MATE.x Module_Setting.o Module_MATE.o Module_BC.o Module_Lyman.o Module_Vol.o Module_ChargeExchange.o main.o init.o trace.o PSD.o IO_utils.o

#MATE_esc.x: Module_MATE.o Module_Setting.o Module_Lyman.o Module_BC.o main_esc.o init.o trace.o PSD.o IO_utils.o
#	$(FC) -o MATE_esc.x Module_MATE.o Module_Setting.o Module_Lyman.o Module_BC.o main_esc.o init.o trace.o PSD.o IO_utils.o


Module_Setting.o: Module_Setting.f90
	$(FC) -c $(flags) Module_Setting.f90

Module_BC.o: Module_BC.f90
	$(FC) -c $(flags) Module_BC.f90

Module_Lyman.o: Module_Lyman.f90
	$(FC) -c $(flags) Module_Lyman.f90

Module_MATE.o: Module_MATE.f90
	$(FC) -c $(flags) Module_MATE.f90

Module_Vol.o: Module_Vol.f90
	$(FC) -c $(flags) Module_Vol.f90

Module_ChargeExchange.o: Module_ChargeExchange.f90
	$(FC) -c $(flags) Module_ChargeExchange.f90

#main_esc.o: main_esc.f90
#	$(FC) -c $(flags) main_esc.f90

main.o: main.f90
	$(FC) -c $(flags) main.f90
trace.o: trace.f90
	$(FC) -c $(flags) trace.f90
init.o: init.f90
	$(FC) -c $(flags) init.f90
IO_utils.o: IO_utils.f90
	$(FC) -c $(flags) IO_utils.f90
PSD.o: PSD.f90
	$(FC) -c $(flags) PSD.f90
	
#extra_tools.o: extra_tools.f90
#	$(FC) -c $(flags) extra_tools.f90


clean:
	rm *.o *.x *.mod
clean2:
	rm *.o *.x *.data *.dat
clean3:
	rm *.o *.x *.data *.dat *__genmod*

