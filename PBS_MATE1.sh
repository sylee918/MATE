#PBS -S /bin/csh
#PBS -N slee122
#PBS -l select=47:ncpus=27:mpiprocs=27:model=bro+1:ncpus=28:mpiprocs=28:model=bro
##PBS -l select=5:ncpus=24:mpiprocs=24:model=bro+1:ncpus=25:mpiprocs=25:model=bro
##PBS -l select=11:ncpus=12:mpiprocs=12:model=bro+1:ncpus=13:mpiprocs=13:model=bro
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -W group_list=s3015
#PBS -m e
#PBS -q devel
##PBS -q normal
##PBS -q long

module load comp-intel mpi-hpe
module load pkgsrc

cd $PBS_O_WORKDIR

#mpirun -np 145 ./MATE.x
mpirun -np 1297 ./MATE.x
