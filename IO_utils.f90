   Subroutine write_density_4D(density_4D,iday)

      USE PHYSICS_TAG
      USE SETTING
      IMPLICIT NONE
      
      real*8 density_4D(nRadial,nLong,nLat_NS,ntperday)
      real, dimension(:,:,:,:), allocatable :: real_density_4D
      integer iday, nlen
      character*10 dayst
      character*100 filename

      allocate(real_density_4D(nRadial,nLong,nLat_NS,ntperday))
      real_density_4D = real(density_4D)

      write(dayst, '(I7.7)') iday
      filename = trim(outdir) // 'MATE_nH_' // trim(tag_phys) // '_' // trim(tag0) // '_' // trim(dayst) // '.data'
      inquire(iolength=nlen) real_density_4D
      open(file=filename,unit=45,form='unformatted',access='direct',recl=nlen,status='replace')
      write(45,rec=1) real_density_4D
      close(45)

      deallocate(real_density_4D)

      return
   End



   Subroutine Make_Parameters_OutFile

      USE PHYSICS_TAG
      USE SETTING
      USE MPI_MATE, only: nprocs
      IMPLICIT NONE

      character*100 filename

      filename = 'MATE_Parameters_' // trim(Runname_in_10char) // '.in'
      open(file=filename,unit=123,status='replace')
      write(123,*) nvel, nRadial, nEnergy
      write(123,*) nRadial, nLon, nLat_NS, ntperday

      write(123,*) "Above paramters are ..."
      write(123,*) "    [nvel, nRadial, nEnergy]"
      write(123,*) "    [nRadial, nLon, nLat_NS, ntperday]"
      write(123,*) "Start_Time_in_YYYYDOY = ", Start_Time_in_YYYYDOY
      write(123,*) "End_Time_in_YYYYDOY   = ", End_Time_in_YYYYDOY
      write(123,*) "Output interval       = ", Output_Time_Interval_in_Minute, " [minutes]"
      write(123,*) "Exobase BC:             ", ExobaseBC_Model_Name
      write(123,*) "Physics:                ", tag_phys
      write(123,*) "MPI:                    ", nprocs

      close(123)


   End