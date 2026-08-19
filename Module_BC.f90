   Module EXOBASE_BC

      USE SETTING
      USE MPI_MATE, only: rank
      IMPLICIT NONE

      real*8, dimension(nbx,nby,nbtperday) :: nH_temp, TH_temp
      real*8, dimension(nbx,nby,nbtperday,start_ydoy-nt_bwd_bc:end_ydoy) :: nH_BC, TH_BC
      character*100 filename_BC

      contains

      Subroutine Get_exobaseBC

         integer iday, maxdoy
         character*7 ydoy_str, yearst

         if (ExobaseBC_Model_Name .eq. "CONST") then
            call Set_Lunar_Surface_BC()
            return
         endif


         if (start_ydoy/1000 .eq. end_ydoy/1000) then
            write(yearst, '(I4.4)') start_ydoy/1000

            do iday=start_ydoy-nt_bwd_bc,end_ydoy
!               nH_temp = 0.d0 ; TH_temp=0.d0
               write(ydoy_str,'(I7.7)') iday
               write(yearst, '(I4.4)') start_ydoy/1000
               filename_BC = trim(BC_dir) // trim(yearst) // "/" // trim(ExobaseBC_Model_Name) // "_" // trim(ydoy_str) //  ".bc"
               call read_exobaseBC
               nH_BC(:,:,:,iday) = nH_temp
               TH_BC(:,:,:,iday) = TH_temp
               if (minval(nH_temp) .lt. 1e-15 .or. minval(TH_temp) .lt. 1e-15) then
                  print*, 'extra_tools', minval(nH_temp), minval(TH_temp)
                  print*, "ERROR: BC has zero values."
                  stop
               endif
            enddo

         else


         maxdoy=365
         if (mod(start_ydoy/1000,4) .eq. 0) then
            maxdoy=366
         endif

         ! eg. 2012360 - 2012366
         do iday=start_ydoy-nt_bwd_bc, (start_ydoy/1000)*1000+maxdoy
!            nH_temp = 0.d0 ; TH_temp=0.d0
            write(ydoy_str,'(I7.7)') iday
            write(yearst, '(I4.4)') start_ydoy/1000
            filename_BC = trim(BC_dir) // trim(yearst) // "/" // trim(ExobaseBC_Model_Name) // "_" // trim(ydoy_str) //  ".bc"
            call read_exobaseBC()
            nH_BC(:,:,:,iday) = nH_temp
            TH_BC(:,:,:,iday) = TH_temp
            if (minval(nH_temp) .lt. 1e-15 .or. minval(TH_temp) .lt. 1e-15) then
               print*, 'extra_tools', minval(nH_temp), minval(TH_temp)
               print*, "ERROR: BC has zero values."
               stop
            endif
         enddo

         ! eg. 2013001 - 2013012
         do iday=(end_ydoy/1000)*1000+1, end_ydoy
!            nH_temp = 0.d0 ; TH_temp=0.d0
            write(ydoy_str,'(I7.7)') iday
            write(yearst, '(I4.4)') end_ydoy/1000
            filename_BC = trim(BC_dir) // trim(yearst) // "/" // trim(ExobaseBC_Model_Name) // "_" // trim(ydoy_str) //  ".bc"
            call read_exobaseBC()
            nH_BC(:,:,:,iday) = nH_temp
            TH_BC(:,:,:,iday) = TH_temp
            if (minval(nH_temp) .lt. 1e-15 .or. minval(TH_temp) .lt. 1e-15) then
               print*, 'extra_tools', minval(nH_temp), minval(TH_temp)
               print*, "ERROR: BC has zero values."
               stop
            endif
         enddo

         endif

         return
      End


      Subroutine read_exobaseBC

         IMPLICIT NONE

         real, dimension(:,:,:), allocatable :: nH_real, TH_real
         integer nlen, thread_num, IO_unit
         logical iexist

         allocate(nH_real(nbx,nby,nbtperday),TH_real(nbx,nby,nbtperday))

         if (rank .eq. 0) print*, "Read exobase BC file: ", filename_BC

         inquire(file=filename_BC, exist=iexist)
         if (.not. iexist) then
            print*, "File is not exist: ", filename_BC
         else
            inquire(iolength=nlen) nH_real
            nlen=nlen*2

            open(file=filename_BC,newunit=IO_unit,form='unformatted', &
               access='direct',action='read',recl=nlen,status='old')
            read(IO_unit,rec=1) nH_real, TH_real
            close(IO_unit)
         endif

         nH_temp = nH_real*1.d0
         TH_temp = TH_real*1.d0
         deallocate(nH_real,TH_real)

         return
      End


      Subroutine Set_Lunar_Surface_BC()

         IMPLICIT NONE

         integer :: ilon, ilat, iday
         real*8  :: lon_deg, lat_deg, lon_rad, lat_rad
         real*8  :: cos_SZA, cos_term, d_lon
         real*8  :: T_dayside, T_term, T_val, n_val
         real*8, parameter :: T_ss = 390.0d0        ! Subsolar temperature (K)
         real*8, parameter :: T_min = 100.0d0       ! Nightside minimum temperature (K)
         real*8, parameter :: decay_lon = 45.0d0    ! Nightside cooling decay scale (deg)
         real*8, parameter :: T0_ref = 120.0d0      ! Reference temperature (K)
         real*8, parameter :: n0_ref = 20.0d0       ! Reference density at T0 (cm^-3)
         real*8, parameter :: rad_factor = pi / 180.0d0

         if (rank .eq. 0) then
            print*, ">>> Initializing Lunar Surface BC (Hurley et al. 2015 + Noncondensable Gas Law)"
            print*, "    T_ss =", T_ss, "K, T_min =", T_min, "K, n(120K) =", n0_ref, "cm^-3"
         endif

         do ilat = 1, nby
            lat_deg = -90.0d0 + (ilat - 0.5d0) * bc_res
            lat_rad = lat_deg * rad_factor

            ! Terminator temperature at dusk (lon = 89.9 deg) for this latitude
            cos_term = cos(lat_rad) * cos(89.9d0 * rad_factor)
            if (cos_term > 0.0d0) then
               T_term = max(T_ss * (cos_term)**0.25d0, T_min)
            else
               T_term = T_min
            endif

            do ilon = 1, nbx
               lon_deg = -180.0d0 + (ilon - 0.5d0) * bc_res
               lon_rad = lon_deg * rad_factor

               ! 1) Dayside (-90 deg <= lon <= +90 deg)
               if (lon_deg >= -90.0d0 .and. lon_deg <= 90.0d0) then
                  cos_SZA = cos(lat_rad) * cos(lon_rad)
                  if (cos_SZA > 0.0d0) then
                     T_dayside = T_ss * (cos_SZA)**0.25d0
                  else
                     T_dayside = 0.0d0
                  endif
                  T_val = max(T_dayside, T_min)

               ! 2) Nightside (lon > 90 deg or lon < -90 deg)
               else
                  if (lon_deg > 90.0d0) then
                     d_lon = lon_deg - 90.0d0
                  else
                     d_lon = lon_deg + 270.0d0
                  endif
                  T_val = (T_term - T_min) * exp(-d_lon / decay_lon) + T_min
               endif

               ! Noncondensable gas law: n * T^2.5 = const => n = n0 * (T0 / T)^2.5
               n_val = n0_ref * (T0_ref / T_val)**2.5d0

               nH_temp(ilon, ilat, :) = n_val
               TH_temp(ilon, ilat, :) = T_val
            enddo
         enddo

         do iday = start_ydoy - nt_bwd_bc, end_ydoy
            nH_BC(:,:,:,iday) = nH_temp
            TH_BC(:,:,:,iday) = TH_temp
         enddo

         if (rank .eq. 0) then
            print*, "    Lunar BC set successfully. Min/Max T (K):", minval(TH_temp), maxval(TH_temp)
            print*, "    Min/Max Density (cm^-3):", minval(nH_temp), maxval(nH_temp)
         endif

         return
      End Subroutine Set_Lunar_Surface_BC

   END MODULE EXOBASE_BC
