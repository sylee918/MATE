   Program main

      USE SETTING
      USE MPI_MATE
      USE SET_VELOCITY_DIRECTION
      USE GRID_PARAMETERS
      USE EXOBASE_BC
      USE SOLAR_LYMAN_ALPHA
      USE CHARGEEXCHANGE
      USE PHYSICS_TAG
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
      integer:: N_REDUCE, dnLon0
      real*8, dimension(7) :: one
      real*8 :: nps1, nH1
      real*8 :: trace_t0, trace_t1, trace_time_total
      real*8 :: calc_t0, calc_t1, calc_time_total
      integer :: global_task_idx, total_valid_tasks
      integer :: tasks_per_proc, remainder
      integer :: my_start_idx, my_end_idx


      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      if (i_Full_3D .eq. 1) call Calculate_Local_nRadial

      if (nprocs .ne. nRadial * (nLon * (nLat-1) + 1) .and. rank .eq. 0) then
         print*, "nprocs", nprocs, "nRadial", nRadial, "nLon", nLon, "nLat", nLat
      endif

      call gen_points_for_NV
      call Init_Parameter
      call Get_exobaseBC
      call Physical_tag
      call read_Lya_Bph  
         if (ExobaseBC_Model_Name .eq. "CONST") then; Lya = 4.d0; endif
         if (i_Photoionization .eq. 0) then; bph = 0.d0; endif
         if (i_Photoionization .eq. 1 .and. ExobaseBC_Model_Name .eq. "CONST") then; bph = 1.5d-7; endif

      if (i_ChargeExchange .eq. 1 .or. i_ChargeExchange .eq. 3) then
         call Read_Plasmasphere
         call Read_Exosphere
      endif
      if (i_ChargeExchange .eq. 2 .or. i_ChargeExchange .eq. 3) then
         call Get_Beta_RCCX() 
      endif
      if (i_ChargeExchange .eq. 0) then
         beta_RCCX = 0.d0
         nps = 0.d0
      endif
      

      if (rank .eq. 0) call Make_Parameters_OutFile()  ! It's not module, just making .in file

      allocate(ptl(nvel,nEnergy,7), flags(nvel,nEnergy))

      total_valid_tasks = 0
      do ilat = 1, nLat_NS
         ! --- 위도별 경도 간격 설정 로직 (기존 코드 유지) ---
         if (ilat .eq. 1 .or. ilat .eq. nLat_NS) then
            nLon0 = 1
         else
            nLon0 = nLong
         endif
    
         dnLon0 = 1 ! 기본값 (Full mode or Equator)
         if (i_Three_Slices .eq. 1) then
!            if (ilat .gt. nLat .and. ilat .lt. nLat_NS) then
            if (ilat .ne. nLat .and. ilat .ne. 1 .and. ilat .ne. nLat_NS) then
               dnLon0 = nLon0 / 4 ! 90도 간격
            endif
         endif
    
         do ilon = 1, nLon0, dnLon0
            do irad = 1, nRadial
               total_valid_tasks = total_valid_tasks + 1
            end do
         end do
      end do

      tasks_per_proc = total_valid_tasks / nProcs
      remainder = mod(total_valid_tasks, nProcs)

      ! 랭크별 시작/끝 인덱스 계산 (나머지 처리 포함)
      if (rank < remainder) then
         tasks_per_proc = tasks_per_proc + 1
         my_start_idx = rank * tasks_per_proc + 1
      else
         my_start_idx = rank * tasks_per_proc + remainder + 1
      endif
      my_end_idx = my_start_idx + tasks_per_proc - 1

      global_task_idx = 0 ! 카운터 초기화

      do iday=start_ydoy, end_ydoy
         number_density_4D_MPI=0.d0; number_density_4D=0.d0
         do it=1,ntperday  ! hour loop

            current_time = iday*1.d0 + it*(time_resolution/86400.d0)
            ihour = it*(time_resolution/3600.d0)
            iminute = it*(time_resolution/60.d0)-ihour*60
            if (rank .eq. 0) print*, 'Current time:', iday, ihour, iminute
            
            global_task_idx = 0

            do ilat=1,nLat_NS
               lat = latitudeNS_range(ilat)
               if (ilat .eq. 1 .or. ilat .eq. nLat_NS) then; nLon0=1; else; nLon0=nLong; endif  ! North & South poles
               dnLon0 = 1
               if (i_Three_Slices .eq. 1) then
!                  if (ilat .gt. nLat .and. ilat .lt. nLat_NS) then
                  if (ilat .ne. nLat .and. ilat .ne. 1 .and. ilat .ne. nLat_NS) then
                     dnLon0=nLon0/4
                  endif
               endif

               do ilon=1,nLon0,dnLon0
                  lon = longitude_range(ilon)
!                  i1 = ilon-1 + (ilat-nLat)*nLong      ! starts from 0
                  do irad=1,nRadial
                     rad = radial_distance_range(irad)
                     global_task_idx = global_task_idx + 1  

                     if (i_Full_3D .eq. 1) then
                        grid_point_idx = i1*nRadial + irad-1
                        if (grid_point_idx >= start_grid .and. grid_point_idx <= end_grid) then
                           print '(a, f5.2, i4, i4)', "(RAD, LON, LAT) = ", rad/Re, int(lon*180/pi), int(lat*180/pi)

                           call Calculate_Density(current_time, number_density_0D)                      
                           number_density_4D_MPI(irad,ilon,ilat,it) = number_density_0D
                        endif
                     endif

                     if (i_Three_Slices .eq. 1) then
!                       ** 핵심: 현재 인덱스가 내 담당 구간인지 확인 **
                        if (global_task_idx >= my_start_idx .and. global_task_idx <= my_end_idx) then
                           ! === [여기에 실제 물리 계산 코드 삽입] ===
                            print *, "Rank", rank, "computing idx:", global_task_idx                    
                           print '(a, i4, i4)', "(LON, LAT) = ", int(lon*180/pi), int(lat*180/pi)

                           call Calculate_Density(current_time, number_density_0D)                      
                           number_density_4D_MPI(irad,ilon,ilat,it) = number_density_0D
                        endif 
                     endif

                  enddo ! irad

               enddo ! ilon
            enddo ! ilat
         enddo ! ihour

         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         N_REDUCE = nRadial * nLon * nLat_NS * ntperday
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

      print*, "maxnH", rank, maxval(number_density_4D), maxval(number_density_4D_MPI)

      deallocate(ptl,flags)

      call MPI_FINALIZE(ierr)

   End Program
