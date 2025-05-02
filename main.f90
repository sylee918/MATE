   Program main

      USE SETTING
      USE INTEGRATED_INITIALIZATION
      USE MPI_MATE
      IMPLICIT NONE

      include "mpif.h"
      external Trace_particle, Calculate_Density
      external write_density_4D, Make_Parameters_OutFile

      real*8, dimension(nRadial,nLon,nLat_NS,ntperday) :: number_density_4D, number_density_4D_MPI
      real*8, allocatable, dimension(:,:,:) :: ptl
      integer, allocatable, dimension(:,:) :: flags
      real*8 :: number_density_0D

      integer doy, iday, ihour, iminute, it, year, hour, nLon0
      real*8 current_time
      integer:: N_REDUCE

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
!      call Calculate_Local_nRadial
!      nR_loc = nRadial
      nR_loc = 1

      call Initialize_Setting

      if (rank .eq. 0) call Make_Parameters_OutFile()  ! It's not module, just making .in file

      allocate(ptl(nvel,nEnergy,7), flags(nvel,nEnergy))

      do iday=start_ydoy, end_ydoy
         number_density_4D_MPI=0.d0; number_density_4D=0.d0
         do it=1,ntperday  ! hour loop
            current_time = iday*1.d0 + it*(time_resolution/86400.d0)
            ihour = it*(time_resolution/3600.d0)
            iminute = it*(time_resolution/60.d0)-ihour*60
            print*, 'Current time:', iday, ihour, iminute

            do ilat=nLat,nLat_NS
               lat = latitudeNS_range(ilat)
               if (ilat .eq. 1 .or. ilat .eq. nLat_NS) then; nLon0=1; else; nLon0=nLong; endif  ! North & South poles
               do ilon=1,nLon0
                  lon = longitude_range(ilon)
                  i1 = (ilon-1 + (ilat-nLat)*nLong)
                  do irad=1,nRadial
                     rad = radial_distance_range(irad)
                     i2 = (i1-1)*nRadial + irad
                     if (rank .eq. i2) then
                        print '(a, f5.2, i3, i3)', "(RAD, LON, LAT) = ", rad/Re, int(lon*180/pi), int(lat*180/pi)

                        call Init_Particles(ptl)
                        call Trace_particle(ptl, flags, current_time)
                        call Calculate_Density(ptl, flags, current_time, number_density_0D)
                        number_density_4D_MPI(irad,ilon,ilat,it) = number_density_0D

                        if (lat .gt. 0) then    ! N/S symmetry
                           ptl(:,:,4) = -ptl(:,:,4)
                           ptl(:,:,7) = -ptl(:,:,7)
                           call Calculate_Density(ptl, flags, current_time, number_density_0D)
                           number_density_4D_MPI(irad,ilon,nLat_NS+1-ilat,it) = number_density_0D
                        endif
                     endif
                  enddo ! irad
               enddo ! ilon
            enddo ! ilat
         enddo ! ihour

         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         N_REDUCE = nR_loc * nLon * nLat_NS * ntperday
         call MPI_REDUCE(number_density_4D_MPI, number_density_4D, N_REDUCE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         
         if (rank .eq. 0) then
            do it=1,ntperday
               do ilon=2,nLong
                  number_density_4D(:,ilon,1,it)       = number_density_4D(:,1,1,it)         ! South pole
                  number_density_4D(:,ilon,nLat_NS,it) = number_density_4D(:,1,nLat_NS,it)   ! North pole
               enddo
            enddo ! it

            call write_density_4D(number_density_4D, iday)
         endif

      enddo ! iday

      deallocate(ptl,flags)

      call MPI_FINALIZE(ierr)

   End Program
