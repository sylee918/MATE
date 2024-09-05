      Program main

      include "mpif.h"
      include "constants.inc"
 
      external Init_Parameter, Init_Particles, Trace_particle, Calculate_Density
      external Get_exobaseBC, read_Lya_Bph, write_density_4D
      external MPI_INIT, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_FINALIZE, MPI_BARRIER, MPI_REDUCE

      real*8, dimension(N_vel_directions,nRadial,nEnergy,7) :: ptl
      integer :: flags(N_vel_directions,nRadial,nEnergy)

      real*8 radial_distance_range(nRadial), energy_range(nEnergy) 
      real*8 longitude_range(nLong), latitude_range(nLat), latitudeNS_range(nLat_NS), radial_boundary(2), tmax
      real*8 lon,lat
      integer ilon, ilat, nLon0

      real*8 number_density_1D(nRadial)
      real*8, dimension(nRadial,nLon,nLat_NS,ntperday) :: number_density_4D, number_density_4D_MPI
      real*8, dimension(nbx,nby,nbtperday,start_ydoy-nt_bwd_bc:end_ydoy) :: nH_BC, TH_BC
      character*70 BC_dir, Lya_dir, outdir
      character*100 filename_BC
      character*30 tag
      integer rank, nprocs, ierr, il, N_REDUCE

      integer doy, iday, ihour, it, year, hour
      real*8, dimension(start_ydoy_index:end_ydoy_index) :: Lya, bph
      real*8 current_time
      character*10 yearst, dayst

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
!      rank=0
!   print*, N_vel_directions

      call Init_Parameter(radial_distance_range, energy_range, longitude_range, latitude_range, latitudeNS_range, radial_boundary, tmax)

      ! load BC
      year = int(start_ydoy/1000)
      write(yearst, '(I4.4)') year
!      BC_dir = '/nobackup/slee122/MATE/MSIS/BC/' // trim(yearst) // '/'
      BC_dir = '/nobackup/slee122/MATE/MSIS/BC/'
      call Get_exobaseBC(BC_dir, nH_BC, TH_BC, rank)
!      print*, 'main', minval(nH_BC), minval(TH_BC)

      ! load Lya
      Lya_dir = '/nobackup/slee122/MATE/OMNI/indices_for_MSIS_1963-2023_2.txt'
      call read_Lya_Bph(Lya_dir, Lya, bph)
      
      do iday=start_ydoy, end_ydoy
         number_density_4D_MPI=0.d0; number_density_4D=0.d0
         do ihour=0,23
!if (rank .eq. int(ihour/2)) then
!   hour = int(ihour/2)*2 + mod(ihour,2)
            it=ihour+1
            current_time = iday*1.d0 + ihour/24.d0
            print*, 'Current time:', iday, ihour

            do ilat=7,nLat_NS
               lat = latitudeNS_range(ilat)
               if (ilat .eq. 1 .or. ilat .eq. nLat_NS) then;  nLon0=1  ! North & South poles
               else;  nLon0=nLong;  endif
               do ilon=1,nLon0
                  lon = longitude_range(ilon)
                  il = ilon-1 + (ilat-1-6)*nLong
                  if (rank .eq. il) then
                     print*, 'LON & LAT = ', int(lon*180/pi), int(lat*180/pi), '[deg]'

                     call Init_Particles(ptl, radial_distance_range, energy_range, lon,lat)
                     call Trace_particle(ptl, flags, radial_boundary, tmax, Lya, current_time)
                     call Calculate_Density(ptl, flags, current_time, nH_BC, TH_BC, number_density_1D, bph, rank)
                     number_density_4D_MPI(:,ilon,ilat,it) = number_density_1D

                     if (lat .gt. 0) then    ! N/S symmetry
                        ptl(:,:,:,4) = -ptl(:,:,:,4)
                        ptl(:,:,:,7) = -ptl(:,:,:,7)
                        call Calculate_Density(ptl, flags, current_time, nH_BC, TH_BC, number_density_1D, bph, rank)
                        number_density_4D_MPI(:,ilon,nLat_NS+1-ilat,it) = number_density_1D
                     endif

                  endif
               enddo ! ilon
            enddo ! ilat
!endif
         enddo ! ihour

         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         N_REDUCE = nRadial * nLon * nLat_NS * ntperday
!         print*, "N_REDUCE: ", N_REDUCE
         call MPI_REDUCE(number_density_4D_MPI, number_density_4D, N_REDUCE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         
         if (rank .eq. 0) then
            do it=1,ntperday
               do ilon=2,nLong
                  number_density_4D(:,ilon,1,it)       = number_density_4D(:,1,1,it)         ! South pole
                  number_density_4D(:,ilon,nLat_NS,it) = number_density_4D(:,1,nLat_NS,it)   ! North pole
               enddo
            enddo ! it

            outdir = "/nobackup/slee122/MATE/0826/GRC/"
            write(dayst, '(I7.7)') iday
            tag = '_GRC_' // trim(dayst)
            call write_density_4D(number_density_4D, outdir, tag)
         endif

      enddo ! iday

      call MPI_FINALIZE(ierr)

      End Program
