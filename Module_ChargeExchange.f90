Module ChargeExchange

   USE SETTING
   IMPLICIT NONE

   integer, parameter :: nx=201, ny=201, nz=201, nh=101, nMLT=24, nphi=24, nrho=101
   real*8, dimension(nh,nMLT,nz) :: nps, Tps
   real*8 :: rho(nh), MLT(nMLT), zps(nz)
   real*8, dimension(nh) :: rho_ps
!   integer, allocatable, dimension(:,:) :: nstep
   real*8, dimension(nRadial,nLon,nLat_NS,ntperday) :: nH0
   real*8, dimension(nRadial_CX,nLon_CX,nLat_CX,ntperday_CX,start_ydoy-nt_bwd_CX:end_ydoy) :: beta_RCCX

contains

   Subroutine Calculate_ChargeExchange(iE,iv,ptl0,flag,current_time, ICX, PSD_CX)

      USE SETTING
      USE MPI_MATE, only: rank
      IMPLICIT NONE

      integer :: iE, iv
      real*8, dimension(7) :: ptl0
      integer :: flag
      real*8, intent(in) :: current_time
      real*8 :: ICX, PSD_CX
      real*8 :: vrel, sigma, vsig_1eV, fac, cexo2, fac2, T_PS_eV, T_PS_K, ICX_i
!      real*8, allocatable :: beta_dt(:), nH_traj(:), vel2(:)
      real*8, dimension(nstep) :: beta_dt, nH_traj, vel2
      integer :: i, istep

      beta_dt = 0.d0; nH_traj = 0.d0; vel2 = 0.d0
      call Trace_Again(iE,iv, ptl0, flag, current_time, beta_dt, nH_traj, vel2, istep)

!      allocate(beta_dt(nstep(iv,iE)), nH_traj(nstep(iv,iE)), vel2(nstep(iv,iE)))

      T_PS_eV = 1.d0 ! eV
      T_PS_K = T_PS_eV * 11604.525 ! K

      vrel = sqrt(2.d0*T_PS_eV*e/mH)*100.d0 ! cm/s
      sigma = 5.d-15 ! cm^2
      vsig_1eV = vrel * sigma

      fac = 2.d0*kb/mH
      cexo2 = fac*T_PS_K
      fac2 = 1.d0/(pi*cexo2)**1.5

      ! PSD_CX is the PSD of CX-created nH. Below is not necessary for RCCX.
!      PSD_CX = 0.d0
!      do i=1,istep
!         ICX_i = sum(beta_dt(1:i)) * vsig_1eV
!         PSD_CX = PSD_CX + abs(beta_dt(i)) * nH_traj(i) * vsig_1eV * exp(-vel2(i)/cexo2) * fac2 * exp(ICX_i)
!      enddo
!      ICX = sum(beta_dt) * vsig_1eV

      ICX = sum(beta_dt)
      ICX = abs(ICX)*(-1.d0)  ! Make sure to be negative.

!      deallocate(beta_dt, nH_traj, vel2)


   End Subroutine Calculate_ChargeExchange



   Subroutine Get_Beta_RCCX()

      USE SETTING
      IMPLICIT NONE

      real*8, dimension(nRadial_CX,nLon_CX,nLat_CX,ntperday_CX) :: beta_ring_current
      character(len=100) :: filename_RC
      character(len=7) :: ydoy_str, yearst
      integer :: iday, nlen, IO_unit
      logical :: iexist

      if (start_ydoy/1000 .eq. end_ydoy/1000) then
!         write(yearst, '(I4.4)') start_ydoy/1000

         do iday=start_ydoy-nt_bwd_CX,end_ydoy
            write(ydoy_str,'(I7.7)') iday
            filename_RC = trim(RCCX_dir) // "RCCX_p_" // trim(ydoy_str) //  ".data"
            call Read_beta_Ring_Current(filename_RC, beta_ring_current)
            beta_RCCX(:,:,:,:,iday) = beta_ring_current
         enddo
      endif

   End Subroutine Get_Beta_RCCX



   Subroutine Read_beta_Ring_Current(filename, beta_ring_current)

      USE SETTING
      IMPLICIT NONE

      real*8, dimension(nRadial_CX,nLon_CX,nLat_CX,ntperday_CX) :: beta_ring_current
      real, allocatable :: beta_ring_current_real(:,:,:,:)
      character(len=100) :: filename
      integer :: i, j, k, nlen, IO_unit
      logical :: iexist

      allocate(beta_ring_current_real(nRadial_CX,nLon_CX,nLat_CX,ntperday_CX))

      IO_unit = 111
      inquire(iolength=nlen) beta_ring_current_real
      open(file=filename,unit=IO_unit,form='unformatted',access='direct',recl=nlen,status='old')
      read(IO_unit,rec=1) beta_ring_current_real
      close(IO_unit)

      beta_ring_current = beta_ring_current_real*1.d0

      deallocate(beta_ring_current_real)

   End Subroutine Read_beta_Ring_Current  



   Subroutine Read_Plasmasphere()

      USE SETTING
      IMPLICIT NONE

      real, allocatable :: nps_real(:,:,:)
      character(len=200) :: filename
      integer :: i, j, k, nlen, IO_unit
      logical :: iexist

      do i=1,nz
         zps(i) = -10.d0 + (i-1)*0.1
      enddo

      do i=1,nh
         rho(i) = (i-1)*0.1
         rho_ps(i) = (i-1)*0.1  ! rho_ps 배열 초기화
      enddo

      do i=1,nMLT
         MLT(i) = (i-1)*15.d0
      enddo

      allocate(nps_real(nh,nMLT,nz))

      IO_unit = 110
      filename = "GCPM_example_cylindrical_kp0.dat"
      inquire(file=filename, exist=iexist)
      if (.not. iexist) then
         print*, "File is not exist: ", filename
      else
         inquire(iolength=nlen) nps_real
         open(file=filename,unit=IO_unit,form='unformatted', &
            access='direct',action='read',recl=nlen,status='old')
         read(IO_unit,rec=1) nps_real
         close(IO_unit)
      endif

      nps = nps_real*1.d0
      deallocate(nps_real)

      Tps = 0.d0 !! FIX ME!!


   End Subroutine Read_Plasmasphere


   Subroutine interpolate_plasmasphere(one, nps1)

      USE GRID_PARAMETERS
      IMPLICIT NONE
      
      real*8, dimension(7) :: one
      real*8 :: nps1
      
      real*8 :: x, y, z, rho, phi, z_coord
      real*8 :: rho_min, rho_max, phi_min, phi_max, z_min, z_max
      real*8 :: phi_grid, phi_grid_next
      integer :: i_rho, i_phi, i_z
      integer :: i_rho1, i_rho2, i_phi1, i_phi2, i_z1, i_z2
      real*8 :: w_rho1, w_rho2, w_phi1, w_phi2, w_z1, w_z2
      real*8 :: nps_interp
      real*8 :: d_rho, d_phi, d_z
      
      ! 입자의 x, y, z 좌표 추출
      x = one(2)/Re
      y = one(3)/Re
      z = one(4)/Re
      
      ! Cylindrical coordinates로 변환
      rho = sqrt(x**2 + y**2)  ! Radial distance from z-axis
      phi = atan2(y, x)        ! Azimuthal angle
      z_coord = z              ! Height

!      print*, "rho, phi, z_coord", rho, phi, z_coord
      
      ! phi를 0-2π 범위로 정규화
      if (phi < 0.d0) phi = phi + 2.d0*pi
      
      ! Grid boundaries 확인
      rho_min = minval(rho_ps)
      rho_max = maxval(rho_ps)
      phi_min = 0.d0
      phi_max = 2.d0*pi
      z_min = minval(zps)
      z_max = maxval(zps)
      
      ! Boundary check - 만약 입자가 plasmasphere grid 밖에 있으면 0 반환
      if (rho < rho_min .or. rho > rho_max .or. &
          z_coord < z_min .or. z_coord > z_max) then
         nps1 = 0.d0
         return
      endif
      
      ! rho grid index 찾기
      i_rho1 = 1
      i_rho2 = nrho
      do i_rho = 1, nrho-1
         if (rho >= rho_ps(i_rho) .and. rho <= rho_ps(i_rho+1)) then
            i_rho1 = i_rho
            i_rho2 = i_rho + 1
            exit
         endif
      enddo
      
      ! phi grid index 찾기 (MLT를 phi로 변환)
      i_phi1 = 1
      i_phi2 = nphi
      do i_phi = 1, nphi-1
         phi_grid = 2.d0*pi * (i_phi-1) / (nphi-1)  ! Grid phi values
         phi_grid_next = 2.d0*pi * i_phi / (nphi-1)
         if (phi >= phi_grid .and. phi <= phi_grid_next) then
            i_phi1 = i_phi
            i_phi2 = i_phi + 1
            exit
         endif
      enddo
      
      ! z grid index 찾기
      i_z1 = 1
      i_z2 = nz
      do i_z = 1, nz-1
         if (z_coord >= zps(i_z) .and. z_coord <= zps(i_z+1)) then
            i_z1 = i_z
            i_z2 = i_z + 1
            exit
         endif
      enddo
      
      ! Interpolation weights 계산
      if (i_rho2 > i_rho1) then
         d_rho = rho_ps(i_rho2) - rho_ps(i_rho1)
         w_rho1 = (rho_ps(i_rho2) - rho) / d_rho
         w_rho2 = (rho - rho_ps(i_rho1)) / d_rho
      else
         w_rho1 = 1.d0
         w_rho2 = 0.d0
      endif
      
      if (i_phi2 > i_phi1) then
         d_phi = 2.d0*pi / (nphi-1)
         phi_grid = 2.d0*pi * (i_phi1-1) / (nphi-1)
         w_phi1 = (phi_grid + d_phi - phi) / d_phi
         w_phi2 = (phi - phi_grid) / d_phi
      else
         w_phi1 = 1.d0
         w_phi2 = 0.d0
      endif
      
      if (i_z2 > i_z1) then
         d_z = zps(i_z2) - zps(i_z1)
         w_z1 = (zps(i_z2) - z_coord) / d_z
         w_z2 = (z_coord - zps(i_z1)) / d_z
      else
         w_z1 = 1.d0
         w_z2 = 0.d0
      endif
      
      ! 3D Trilinear interpolation 수행
      nps_interp = 0.d0
      
      ! 8개 corner points에 대한 interpolation
      nps_interp = nps_interp + &
                   w_rho1 * w_phi1 * w_z1 * nps(i_rho1, i_phi1, i_z1) + &
                   w_rho2 * w_phi1 * w_z1 * nps(i_rho2, i_phi1, i_z1) + &
                   w_rho1 * w_phi2 * w_z1 * nps(i_rho1, i_phi2, i_z1) + &
                   w_rho2 * w_phi2 * w_z1 * nps(i_rho2, i_phi2, i_z1) + &
                   w_rho1 * w_phi1 * w_z2 * nps(i_rho1, i_phi1, i_z2) + &
                   w_rho2 * w_phi1 * w_z2 * nps(i_rho2, i_phi1, i_z2) + &
                   w_rho1 * w_phi2 * w_z2 * nps(i_rho1, i_phi2, i_z2) + &
                   w_rho2 * w_phi2 * w_z2 * nps(i_rho2, i_phi2, i_z2)
      
      nps1 = nps_interp
      
   End Subroutine interpolate_plasmasphere


   Subroutine nearest_grid_plasmasphere(one, nps1)

      USE SETTING
      USE GRID_PARAMETERS
      IMPLICIT NONE
      
      real*8, dimension(7) :: one
      real*8 :: nps1
      
      real*8 :: x, y, z, rho, phi, z_coord
      real*8 :: rho_min, rho_max, phi_min, phi_max, z_min, z_max
      integer :: i_rho_nearest, i_phi_nearest, i_z_nearest
      real*8 :: drho, dphi, dz
      
      ! 입자의 x, y, z 좌표 추출
      x = one(2)/Re
      y = one(3)/Re
      z = one(4)/Re
      
      ! Cylindrical coordinates로 변환
      rho = sqrt(x**2 + y**2)  ! Radial distance from z-axis
      phi = atan2(y, x)        ! Azimuthal angle
      z_coord = z              ! Height

      rho_min = minval(rho_ps)
      rho_max = maxval(rho_ps)
      phi_min = 0.d0
      phi_max = 2.d0*pi
      z_min = minval(zps)
      z_max = maxval(zps)

      ! phi를 0-2π 범위로 정규화
      if (phi < 0.d0) phi = phi + 2.d0*pi
      
      drho=0.1
      i_rho_nearest = nint((rho-rho_min)/drho) + 1
      if (i_rho_nearest .le. 0) i_rho_nearest = 1
      if (i_rho_nearest .gt. nrho) then
         nps1 = 0.d0
         return
      endif
      
      dphi=2.d0*pi/nphi
      i_phi_nearest = nint((phi-phi_min)/dphi) + 1
      if (i_phi_nearest .le. 0) i_phi_nearest = i_phi_nearest + nphi
      if (i_phi_nearest .gt. nphi) i_phi_nearest = i_phi_nearest - nphi

      dz=0.1
      i_z_nearest = nint((z_coord-z_min)/dz) + 1
      if (i_z_nearest .le. 0 .or. i_z_nearest .gt. nz) then
         nps1 = 0.d0
         return
      endif

!      print*, "i_rho_nearest", i_rho_nearest, i_phi_nearest, i_z_nearest
   
      nps1 = nps(i_rho_nearest, i_phi_nearest, i_z_nearest)

      return
      
   End Subroutine nearest_grid_plasmasphere


   Subroutine Read_Exosphere()

      USE SETTING
      IMPLICIT NONE

      character(len=200) :: filename
      integer :: i, j, k, nlen, IO_unit
      logical :: iexist
      real*4, allocatable :: nH_temp(:,:,:,:)

      allocate(nH_temp(nRadial,nLon,nLat_NS,ntperday))

      IO_unit = 120
!      filename = trim(outdir)//"MATE_nH_GRCX_CXtest1_1000004.data"
      filename = trim(outdir)//"MATE_nH_GRC_00_1000008.data"
      inquire(file=filename, exist=iexist)
      if (.not. iexist) then
         print*, "File is not exist: ", filename
      else
         inquire(iolength=nlen) nH_temp
         open(file=filename,unit=IO_unit,form='unformatted', &
            access='direct',action='read',recl=nlen,status='old')
         read(IO_unit,rec=1) nH_temp
         close(IO_unit)
      endif
      
      nH0 = nH_temp*1.d0
      deallocate(nH_temp)

      return

   End Subroutine Read_Exosphere


   Subroutine interpolate_exosphere(one, nH1)

      USE SETTING
      USE GRID_PARAMETERS
      IMPLICIT NONE
      
      real*8, dimension(7) :: one
      real*8 :: nH1
      
      real*8 :: x, y, z, r, longitude, latitude
      real*8 :: r_min, r_max, lon_min, lon_max, lat_min, lat_max
      integer :: i_r, i_lon, i_lat, i_time
      integer :: i_r1, i_r2, i_lon1, i_lon2, i_lat1, i_lat2
      real*8 :: w_r1, w_r2, w_lon1, w_lon2, w_lat1, w_lat2
      real*8 :: nH_interp
      real*8 :: d_r, d_lon, d_lat
      
      ! Read exosphere data
      !call Read_Exosphere(nH0)
!      print*, "maxval(nH0)", maxval(nH0)
      
      ! 입자의 x, y, z 좌표 추출
      x = one(2)/Re
      y = one(3)/Re
      z = one(4)/Re
      
      ! Spherical coordinates로 변환
      r = sqrt(x**2 + y**2 + z**2)           ! Radial distance from Earth center
      longitude = atan2(y, x)                 ! Longitude (0 to 2π)
      latitude = asin(z/r)                    ! Latitude (-π/2 to π/2)

!      print*, "r, longitude, latitude", r, longitude*180/pi, latitude*180/pi
      
      ! longitude를 0-2π 범위로 정규화
      if (longitude < 0.d0) longitude = longitude + 2.d0*pi
      
      ! Grid boundaries 확인
      r_min = minval(radial_distance_range)/Re
      r_max = maxval(radial_distance_range)/Re
      lon_min = 0.d0
      lon_max = 2.d0*pi
      lat_min = minval(latitudeNS_range)
      lat_max = maxval(latitudeNS_range)
      
      ! Boundary check - 만약 입자가 exosphere grid 밖에 있으면 0 반환
      if (r < r_min .or. r > r_max .or. &
          longitude < lon_min .or. longitude > lon_max .or. &
          latitude < lat_min .or. latitude > lat_max) then
         nH1 = 0.d0
         return
      endif
      
      ! r grid index 찾기
      i_r1 = 1
      i_r2 = nRadial
      do i_r = 1, nRadial-1
         if (r >= radial_distance_range(i_r)/Re .and. r <= radial_distance_range(i_r+1)/Re) then
            i_r1 = i_r
            i_r2 = i_r + 1
            exit
         endif
      enddo
      
      ! longitude grid index 찾기
      i_lon1 = 1
      i_lon2 = nLon
      do i_lon = 1, nLon-1
         if (longitude >= longitude_range(i_lon) .and. longitude <= longitude_range(i_lon+1)) then
            i_lon1 = i_lon
            i_lon2 = i_lon + 1
            exit
         endif
      enddo
      
      ! latitude grid index 찾기 (latitude 인덱스가 작으면 남극, 최대는 북극)
      i_lat1 = 1
      i_lat2 = nLat_NS
      do i_lat = 1, nLat_NS-1
         if (latitude >= latitudeNS_range(i_lat) .and. latitude <= latitudeNS_range(i_lat+1)) then
            i_lat1 = i_lat
            i_lat2 = i_lat + 1
            exit
         endif
      enddo
      
      ! Interpolation weights 계산
      if (i_r2 > i_r1) then
         d_r = (radial_distance_range(i_r2) - radial_distance_range(i_r1))/Re
         w_r1 = ((radial_distance_range(i_r2)/Re) - r) / d_r
         w_r2 = (r - (radial_distance_range(i_r1)/Re)) / d_r
      else
         w_r1 = 1.d0
         w_r2 = 0.d0
      endif
      
      if (i_lon2 > i_lon1) then
         d_lon = longitude_range(i_lon2) - longitude_range(i_lon1)
         w_lon1 = (longitude_range(i_lon2) - longitude) / d_lon
         w_lon2 = (longitude - longitude_range(i_lon1)) / d_lon
      else
         w_lon1 = 1.d0
         w_lon2 = 0.d0
      endif
      
      if (i_lat2 > i_lat1) then
         d_lat = latitudeNS_range(i_lat2) - latitudeNS_range(i_lat1)
         w_lat1 = (latitudeNS_range(i_lat2) - latitude) / d_lat
         w_lat2 = (latitude - latitudeNS_range(i_lat1)) / d_lat
      else
         w_lat1 = 1.d0
         w_lat2 = 0.d0
      endif
      
      ! 3D Trilinear interpolation 수행 (시간은 첫 번째 시간 스텝 사용)
      i_time = 1
      nH_interp = 0.d0
      
      ! 8개 corner points에 대한 interpolation
      nH_interp = nH_interp + &
                   w_r1 * w_lon1 * w_lat1 * nH0(i_r1, i_lon1, i_lat1, i_time) + &
                   w_r2 * w_lon1 * w_lat1 * nH0(i_r2, i_lon1, i_lat1, i_time) + &
                   w_r1 * w_lon2 * w_lat1 * nH0(i_r1, i_lon2, i_lat1, i_time) + &
                   w_r2 * w_lon2 * w_lat1 * nH0(i_r2, i_lon2, i_lat1, i_time) + &
                   w_r1 * w_lon1 * w_lat2 * nH0(i_r1, i_lon1, i_lat2, i_time) + &
                   w_r2 * w_lon1 * w_lat2 * nH0(i_r2, i_lon1, i_lat2, i_time) + &
                   w_r1 * w_lon2 * w_lat2 * nH0(i_r1, i_lon2, i_lat2, i_time) + &
                   w_r2 * w_lon2 * w_lat2 * nH0(i_r2, i_lon2, i_lat2, i_time)
      
      nH1 = nH_interp
      
   End Subroutine interpolate_exosphere


   Subroutine nearest_grid_exosphere(current_time, one, nH1, beta_RCCX1)

      USE SETTING
      USE GRID_PARAMETERS
      IMPLICIT NONE
      
      real*8, intent(in) :: one(7), current_time 
      real*8, intent(out) :: nH1, beta_RCCX1
      
      real*8 :: x, y, z, r, longitude, latitude
      real*8 :: r_min, r_max, lon_min, lon_max, lat_min, lat_max
      integer :: i_r_nearest, i_lon_nearest, i_lat_nearest
      real*8 :: dr1, dlon1, dlat1

      integer :: it_nearest, iday, iyear
      real*8 :: year_doy_frac, frac_day, hour_val
      integer :: days_in_year
      logical :: is_leap_year

      ! Read exosphere data
!      call Read_Exosphere(nH0)
      
      ! 입자의 x, y, z 좌표 추출
      x = one(2)/Re
      y = one(3)/Re
      z = one(4)/Re
      
      ! Spherical coordinates로 변환
      r = sqrt(x**2 + y**2 + z**2)           ! Radial distance from Earth center
      longitude = atan2(y, x)                 ! Longitude (0 to 2π)
      latitude = asin(z/r)                    ! Latitude (-π/2 to π/2)

      r_min = minval(radial_distance_range)/Re
      r_max = maxval(radial_distance_range)/Re
      lon_min = 0.d0
      lon_max = 2.d0*pi
      lat_min = minval(latitudeNS_range)
      lat_max = maxval(latitudeNS_range)

      ! longitude를 0-2π 범위로 정규화
      if (longitude < 0.d0) longitude = longitude + 2.d0*pi
      
     
      ! r grid index 찾기 (nearest grid point)
      dr1=0.5
      i_r_nearest = nint((r-r_min)/dr1) + 1
      if (i_r_nearest .lt. 1) i_r_nearest = 1
      if (i_r_nearest .gt. nRadial) then
         nH1 = 0.d0
         return
      endif

      dlon1=2.d0*pi/nLon
      i_lon_nearest = nint((longitude-lon_min)/dlon1) + 1
      if (i_lon_nearest .le. 0) i_lon_nearest = i_lon_nearest + nLon
      if (i_lon_nearest .gt. nLon) i_lon_nearest = i_lon_nearest - nLon

      dlat1=pi/(nLat_NS-1)
      i_lat_nearest = nint((latitude-lat_min)/dlat1) + 1

      year_doy_frac = mod(current_time, 1000.d0)
      iday = int(year_doy_frac)
      
      frac_day = year_doy_frac - dble(iday)
      
      hour_val = frac_day * 24.d0
      
      it_nearest = nint(hour_val) + 1
      
      ! 경계 처리: 23:30 이상(hour_val > 23.5)이 되어 반올림으로 25가 될 경우
      ! 0시(index 1)로 순환
      if (it_nearest > 24) then
         it_nearest = 1
         iday=iday+1
      endif
      if (it_nearest < 1) then
         it_nearest = 24
         iday=iday-1
      endif
      if (iday < Beta_CX_Start_Time_in_YYYYDOY) then
         iday = Beta_CX_Start_Time_in_YYYYDOY
         it_nearest = 1
      endif

      ! -----------------------------------------------------------


!      nH1 = nH0(i_r_nearest, i_lon_nearest, i_lat_nearest, it_nearest)
!      print*, 'nearest', i_r_nearest, i_lon_nearest, i_lat_nearest, it_nearest, iday
      beta_RCCX1 = beta_RCCX(i_r_nearest, i_lon_nearest, i_lat_nearest, it_nearest, iday)

      return

   End Subroutine nearest_grid_exosphere



   Subroutine Trace_Again(iE,iv, ptl0,flag, current_time, beta_dt, nH_traj, vel2, istep)

      USE SETTING
      USE MPI_MATE, only: rank
      USE GRID_PARAMETERS, only: radial_boundary
      USE SOLAR_LYMAN_ALPHA, only: Lya
      IMPLICIT NONE
      external rk4, calculate_final_timestep

      integer :: iE, iv, istep
      real*8, dimension(7) :: ptl0
      real*8, dimension(7) :: one, old
      real*8, intent(in) :: current_time
!      real*8, dimension(nstep(iv,iE)) :: beta_dt, nH_traj, vel2
      real*8, dimension(nstep) :: beta_dt, nH_traj, vel2
      integer :: flag     ! 0: orbiting Earth t<tmax;   1: into exobase;  2: out of outer boundary;  3: orbiting but t>tmax
      real*8 :: radial_distance, radial_distance_old
      integer :: i
      real*8 :: dt, vt,vt_old,dv, ds
      real*8, parameter :: max_ds = 1.d6
      real*8 :: x0, f0, trace_time
      integer :: ydoy, ii
      real*8 ::  nps1, nH1, beta_RCCX1


      istep = 0; flag=0
      one = ptl0

      radial_distance_old = sqrt(one(2)**2+one(3)**2+one(4)**2)
      vt = sqrt(one(5)**2 + one(6)**2 + one(7)**2)
      trace_time = current_time - 1e-5

      do while (abs(one(1)) < tmax)

         istep = istep + 1
         radial_distance_old = sqrt(one(2)**2+one(3)**2+one(4)**2)
         vt_old = sqrt(one(5)**2 + one(6)**2 + one(7)**2)

         dt = -1.d0*max_ds / vt_old     ! -1e6 or 4e6 is a "factor" in python code. The maximum distance jump at single time step.
         if (mod(trace_time-1,1000.0) .gt. 500) then               ! eg. trace_time=2010000.98, then it should be 2009365.98, 
            ii=1000-mod(int(trace_time-1),1000)                  ! eg. trace_time-1 = 2009999.98, ii=1000-999=1
            if (mod(int((trace_time-1)/1000),4) .eq. 0) then
               trace_time = int((trace_time-1)/1000)*1000 + (367-ii) + mod(trace_time,1.0)      ! For leap years (400-year period is not applied)
            else
               trace_time = int((trace_time-1)/1000)*1000 + (366-ii) + mod(trace_time,1.0)      ! eg. trace_time = 2009000+365+0.98 = 2009365.98
            endif
         endif

         ydoy = int(trace_time)
         f0 = Lya(ydoy)
         if (int(trace_time + dt/86400) .ne. int(trace_time)) then
            dt = (int(trace_time)-trace_time)*86400 - 1e-5    ! trace_time always hits the time (00:00:00) for daily-varying Lya.
            if (abs(dt) .lt. 1e-6) then
               print*, "ERROR: dt is too small"
               stop
            endif
         endif

         old = one
   101 continue
         call rk4(one,dt,f0)
         trace_time = current_time + one(1)/86400  ! one(2) < 0
         ! FIX ME (if time cross year)

         radial_distance = sqrt(one(2)**2+one(3)**2+one(4)**2)
         vt = sqrt(one(5)**2 + one(6)**2 + one(7)**2)
         ds = sqrt((old(2)-one(2))**2+(old(3)-one(3))**2+(old(4)-one(4))**2)
         dv = vt-vt_old

         ! If the solution is diverging ...
         if (abs(ds/radial_distance_old) .gt. 1e-1 .or. abs(ds)/max_ds .gt. 1.2 .or. dv/vt_old .gt. 10) then
               dt=dt/2
               one=old
               goto 101
         endif

!         call nearest_grid_plasmasphere(one, nps1)
!         call nearest_grid_exosphere(one, nH1)
!         call nearest_grid_exosphere(current_time, one, nH1, beta_RCCX1)
         call nearest_grid_exosphere(trace_time, one, nH1, beta_RCCX1)
         !call interpolate_plasmasphere(one, nps1)
         !call interpolate_exosphere(one, nH1)
!         beta_dt(istep) = nps1*dt
         beta_dt(istep) = beta_RCCX1*dt
         nH_traj(istep) = nH1
         vel2(istep) = vt**2
!            print*, "beta_CX1, nps1, dt", beta_CX1, nps1, dt, rank
!            stop

         radial_distance = sqrt(one(2)**2+one(3)**2+one(4)**2)
         if (radial_distance .lt. radial_boundary(1)) then
            flag = 1
            call calculate_final_timestep(old,one,dt,f0)
               istep = istep + 1
!               call nearest_grid_plasmasphere(one, nps1)
               call nearest_grid_exosphere(trace_time, one, nH1, beta_RCCX1)
               !call interpolate_plasmasphere(one, nps1)
               !call interpolate_exosphere(one, nH1)
!               beta_dt(istep) = nps1*dt
               beta_dt(istep) = beta_RCCX1*dt
               nH_traj(istep) = nH1
               vel2(istep) = vt**2
         else if (radial_distance .gt. radial_boundary(2)) then
            flag = 2
         endif

         if (flag > 0) then
            exit
         endif

         if (istep .ge. nstep) then
            print*, "WARNING: istep >= nstep"
         endif

      enddo ! end while

   ptl0 = one

   return
End Subroutine Trace_Again

End Module ChargeExchange