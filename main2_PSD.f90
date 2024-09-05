   Program main2_PSD

!      use omp_lib
      include "mpif.h"
      include "constants.inc"
      external Init_Parameter, Initialize, Calculate_Density_Of_Each_Particle, Generate_tag
      external read_exobaseBC, write_density_4D
      external MPI_INIT, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_FINALIZE, MPI_REDUCE, MPI_BARRIER

      real*8, dimension(N_vel_directions,nRadial,nEnergy,7) :: init, fin
      integer, dimension(N_vel_directions,nRadial,nEnergy) :: flags
!      real*8, dimension(:,:,:,:), allocatable :: init, fin
!      integer, dimension(:,:,:), allocatable :: flags
      character*70 input_dir, BC_dir, outdir
      character*100 filename_BC
      character*30 tag
      real*8 radial_distance_range(nRadial), energy_range(nEnergy)
      real*8 longitude_range(nLong), latitude_range(nLat), latitudeNS_range(nLat_NS)
      real*8 radial_boundary(2), tmax
      real*8 number_density_1D(nRadial)
!      real*8, dimension(nRadial,nLon,nLat_NS) :: number_density_3D, number_density_3D_MPI
      real*8, dimension(nRadial,nLon,nLat_NS,ntperday) :: number_density_4D, number_density_4D_MPI
      real*8, dimension(nbx,nby,nbt) :: nH_BC, TH_BC
      real*8 lon, lat, current_time
      integer ilon, ilat, nLon0, it, iday
      integer rank, nprocs, ierr, N_REDUCE
      character*10 dayst

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

!      call OMP_set_num_threads(2)
      
!      input_dir = "/nobackup/slee122/MATE/0212/f03/GR/"
      input_dir = "/nobackup/slee122/MATE/0212/G/"
!      input_dir = "/u/slee122/MATE/f03/GRC/"
      call Init_Parameter(radial_distance_range, energy_range, longitude_range, latitude_range, latitudeNS_range, radial_boundary,tmax)

      filename_BC = '/nobackup/slee122/MATE/MSIS/MSIS_20200101-20200110_5deg_5min.bc'
      call read_exobaseBC(filename_BC, nH_BC,TH_BC, rank+12)

      do iday=1,10

      number_density_4D_MPI=0.d0
      do ilat=1,nLat_NS
         if (rank .eq. ilat-1) then 

            lat = latitudeNS_range(ilat)
            if (ilat .eq. 1 .or. ilat .eq. nLat_NS) then; nLon0=1  ! North & South poles
            else;  nLon0=nLong;  endif

            do ilon=1,nLon0
               lon = longitude_range(ilon)
               print*, 'rank:', rank, ', LON & LAT = ', int(lon*180/pi), int(lat*180/pi), '[deg]'

               call Generate_tag(lon,abs(lat), tag)
               print*, tag

               init=0.d0; fin=0.d0; flags=0
               call Initialize(input_dir, init, fin, flags, tag, rank)
               if (lat .lt. 0) then
                  print*, 'rank:', rank, ', LAT = ', int(lat*180/pi), tag
                  fin(:,:,:,4) = -fin(:,:,:,4)
                  fin(:,:,:,7) = -fin(:,:,:,7)
               endif

               do it=1,ntperday
                  current_time = it*tres + (iday-1)*86400
                  print*, 'rank:', rank, ', Time = ', current_time, '[s]'

                  call Calculate_Density_Of_Each_Particle(fin, flags, current_time, nH_BC, TH_BC, number_density_1D, rank)
                  number_density_4D_MPI(:,ilon,ilat,it) = number_density_1D
               enddo ! it

            enddo ! ilon
         endif ! rank
      enddo ! ilat

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      N_REDUCE = nRadial * nLon * nLat_NS * ntperday
      print*, "N_REDUCE: ", N_REDUCE
      call MPI_REDUCE(number_density_4D_MPI, number_density_4D, N_REDUCE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      
      if (rank .eq. 0) then
         do it=1,ntperday
            do ilon=2,nLong
               number_density_4D(:,ilon,1,it)       = number_density_4D(:,1,1,it)         ! South pole
               number_density_4D(:,ilon,nLat_NS,it) = number_density_4D(:,1,nLat_NS,it)   ! North pole
            enddo
         enddo ! it

         outdir = "/nobackup/slee122/MATE/0418/"
         write(dayst, '(I3.3)') iday
         tag = '_G_' // trim(dayst)
         call write_density_4D(number_density_4D, outdir, tag)
      endif

      enddo ! iday

      call MPI_FINALIZE(ierr)

   End Program
