      Program main

      use Module_for_NVelocityDirection
      use Module_Physics_tag

      include "mpif.h"
      include "Setting.inc"

      external Init_Parameter, Init_Particles, Trace_particle
      external Get_exobaseBC, read_Lya_Bph, write_density_4D
      external MPI_INIT, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_FINALIZE, MPI_BARRIER, MPI_REDUCE
      external Forward_Tracing_particle, Calculate_Escaping_Flux_constBC, Write_ESC_FLUX_2D
      external Make_Parameters_OutFile

      real*8, allocatable, dimension(:,:,:,:) :: b_ptl, f_ptl
      integer, allocatable, dimension(:,:,:) :: b_flags, f_flags

      real*8 radial_distance_range(nRadial), energy_range(nEnergy) 
      real*8 longitude_range(nLong), latitude_range(nLat), latitudeNS_range(nLat_NS), radial_boundary(2), tmax
      real*8 lon,lat
      integer ilon, ilat, nLon0

      real*8, dimension(nbx,nby) :: esc_flux, esc_flux_MPI
      real*8, dimension(nbx,nby,nbtperday,start_ydoy-nt_bwd_bc:end_ydoy) :: nH_BC, TH_BC
      character*30 tag
      integer rank, nprocs, ierr, il, N_REDUCE

      integer doy, iday, ihour, iminute, it, year, hour
      real*8, dimension(start_ydoy_index:end_ydoy_index) :: Lya, bph
      real*8 current_time
      character*10 yearst, dayst

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      call Init_Parameter(radial_distance_range, energy_range, longitude_range, latitude_range, latitudeNS_range, radial_boundary, tmax)
      ! load BC
      year = int(start_ydoy/1000);  write(yearst, '(I4.4)') year
      Lya = 3.d0 ; bph=0.d0
!      call Get_exobaseBC(nH_BC, TH_BC, rank)
!      call read_Lya_Bph(Lya, bph);  if (i_Photoionization .eq. 0) bph = 0.d0
      ! call modules
      call Physics_tag(); call gen_points_for_NV(); print*, "N_vel_directions = ", N_vel_directions
      call Make_Parameters_OutFile()
      ! End Initialization

      allocate(b_ptl(N_vel_directions,nRadial,nEnergy,7)) ; allocate(b_flags(N_vel_directions,nRadial,nEnergy))
      allocate(f_ptl(N_vel_directions,nRadial,nEnergy,7)) ; allocate(f_flags(N_vel_directions,nRadial,nEnergy))

      esc_flux=0.d0; esc_flux_MPI=0.d0
!      do iday=start_ydoy, end_ydoy
!         number_density_2D_MPI=0.d0; number_density_2D=0.d0
!         do it=1,ntperday
!            current_time = iday*1.d0 + it*(time_resolution/86400.d0)
!            ihour = it*(time_resolution/3600.d0)
!            iminute = it*(time_resolution/60.d0)-ihour*60
!            print*, 'Current time:', iday, ihour, iminute

            do ilat=nLat,nLat_NS
               lat = latitudeNS_range(ilat)
               if (ilat .eq. 1 .or. ilat .eq. nLat_NS) then; nLon0=1; else; nLon0=nLong; endif  ! North & South poles
               do ilon=1,nLon0
                  lon = longitude_range(ilon)
                  il = ilon-1 + (ilat-nLat)*nLong
                  if (rank .eq. il) then
                     print*, '  LON & LAT = ', int(lon*180/pi), int(lat*180/pi), '[deg]'

                     call Init_Particles(b_ptl, radial_distance_range, energy_range, lon,lat)
                     f_ptl = b_ptl
                     call Trace_particle(b_ptl, b_flags, radial_boundary, tmax, Lya, current_time)
                     call Forward_Tracing_particle(f_ptl,f_flags, b_flags, radial_boundary, tmax, Lya, current_time)
                     call Calculate_Escaping_Flux_constBC(b_ptl, f_ptl, f_flags, esc_flux_MPI, rank)

                     if (lat .gt. 0) then    ! N/S symmetry
                        b_ptl(:,:,:,4) = -b_ptl(:,:,:,4)
                        b_ptl(:,:,:,7) = -b_ptl(:,:,:,7)
                        call Calculate_Escaping_Flux_constBC(b_ptl, f_ptl, f_flags, esc_flux_MPI, rank)
                     endif

                  endif
               enddo ! ilon
            enddo ! ilat
!         enddo ! ihour

         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         N_REDUCE = nRadial * nLon * nLat_NS * ntperday
         call MPI_REDUCE(esc_flux_MPI, esc_flux, N_REDUCE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         
         if (rank .eq. 0) then
!            do it=1,ntperday
!               do ilon=2,nLong
!                  esc_flux(ilon,1)       = esc_flux(1,1)         ! South pole
!                  esc_flux(ilon,nLat_NS) = esc_flux(1,nLat_NS)   ! North pole
!               enddo
!            enddo ! it

            write(dayst, '(I7.7)') iday
            tag = '_' // trim(tag0)
            call Write_ESC_FLUX_2D(esc_flux, tag)
         endif

!      enddo ! iday

      deallocate(b_ptl,b_flags,f_ptl,f_flags)

      call MPI_FINALIZE(ierr)

      End Program
