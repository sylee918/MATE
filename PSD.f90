   Subroutine Calculate_Density(current_time, number_density_0D)
      ! "cdensity" in python code
!         use omp_lib
      USE SETTING
      USE TIME_UTILS, only: ydoy_add_days, ydoy_diff_days, ydoy_add_days_int
      USE EXOBASE_BC, only: nH_BC, TH_BC, total_bc_days
      USE SOLAR_LYMAN_ALPHA, only: bph
      USE MPI_MATE, only: nR_loc
      USE VOLUME_ELEMENT
      USE ChargeExchange
      use, intrinsic :: ieee_arithmetic
      external GSE2SPH

      real*8, dimension(nvel,nEnergy,7) :: ptl
      integer, dimension(nvel,nEnergy) :: flags
      real*8, intent(in) :: current_time
      real*8, intent(out) :: number_density_0D

      real*8, dimension(7) :: ptl0
      real*8, dimension(:,:), allocatable :: each_n
      real*8 pos(3), vel(3), vel2
      real*8 temp_BC, n_BC, vel_BC(3), fac, PhaseSpaceDensity
      integer iE,iv, i
      real*8 finlon, finlat, cexo2
       real*8 t0, t1, Iph, ICX, PSD_CX, PSD_exobase
!       real*8 :: cx_t0, cx_t1, cx_time_total
!       integer :: cx_calls
      integer iflon, iflat, it, quotient
      integer idoy, iday, day_idx, cur_day


      vel_BC = 0.d0;
      allocate(dV2(nEnergy,nvel), solid_angle(nvel))
      call calculate_Velocity_Volume_Element(dV2)

      allocate(each_n(nvel,nEnergy))
      each_n = 0.d0
      number_density_0D = 0.d0
!      cx_time_total = 0.d0
!      cx_calls = 0

      fac = 2.d0*kb/mH
      call Init_Particles(ptl)

      do iE=1,nEnergy
         do iv=1,nvel

            ! Trace here
            ptl0 = ptl(iv,iE,:)
            call Calculate_ChargeExchange(iE,iv, ptl0, flags(iv,iE), current_time, ICX, PSD_CX)
            ptl(iv,iE,:) = ptl0

            t0 = ydoy_add_days(current_time, ptl(iv,iE,1)/86400.d0)    ! unit day
            idoy = int(t0)                               ! yyyy+doy
            t1 = (t0 - idoy)*86400.d0                    ! hms in seconds
            it = floor(t1/tb_res)+1
            if (it < 1) it = 1
            if (it > nbtperday) it = nbtperday
            if (idoy .lt. BC_Start_Time_in_YYYYDOY) then
               idoy = BC_Start_Time_in_YYYYDOY
               it = 1
            endif



            if (flags(iv,iE) .eq. 1) then ! Exobase-origin particle
               do i=1,3
                  pos(i) = ptl(iv,iE,i+1)
                  vel(i) = ptl(iv,iE,i+4)
               enddo

               call GSE2SPH(pos,finlon,finlat)
               iflon=floor(finlon/bc_res)+(180/bc_res)+1              
               iflat=floor(finlat/bc_res)+(90/bc_res)+1     ! -90 < lat < 90
               if (iflat .eq. 180/bc_res+1) then
                  iflon = iflon + 180/bc_res
                  iflat = 180/bc_res
               endif
               if (iflon .ge. 360/bc_res+1) then
                  quotient = int(iflon/(360/bc_res))
                  iflon = iflon - (360/bc_res)*quotient
               endif

               day_idx = ydoy_diff_days(idoy, BC_Start_Time_in_YYYYDOY) + 1
               if (day_idx < 1) day_idx = 1
               if (day_idx > total_bc_days) day_idx = total_bc_days

               n_BC    = nH_BC(iflon,iflat,it,day_idx)
               temp_BC = TH_BC(iflon,iflat,it,day_idx)

               !! ** FIX ME (above): Trilinear interpolation is desired for more accurate calculation.
               !!                    Current code is just the 0th-order interpolation.
            else ! Plasmasphere-origin particle
               n_BC = 0.d0
               temp_BC = 1.d0  !! Set any finite value to avoid NaN
            endif

               if (i_Photoionization .eq. 1) then
                  if (idoy .eq. int(current_time)) then
                     Iph = bph(idoy) * abs(ptl(iv,iE,1))
                  else
                     Iph = bph(idoy) * (86400.d0 - t1)
                     cur_day = ydoy_add_days_int(idoy, 1)
                     do while (cur_day < int(current_time))
                        Iph = Iph + bph(cur_day)*86400.d0
                        cur_day = ydoy_add_days_int(cur_day, 1)
                     enddo
                     Iph = Iph + bph(int(current_time)) * (current_time-int(current_time))*86400.d0
                  endif
               else
                  Iph = 0.d0
               endif

               cexo2 = fac*temp_BC
               vel2 = sum(vel*vel)
               if (i_ChargeExchange .eq. 0) then
                  ICX = 0.d0     ! loss rate by CX
                  PSD_CX = 0.d0  ! PSD of CX-created H
               endif
               if (i_ChargeExchange .eq. 2) then  ! RCCX only
                  PSD_CX = 0.d0
               endif

               if (ICX .gt. 0.d0) then
                  print*, "ICX is positive", ICX
                  stop
               endif

               !vel = (vel - vel_BC)
               PSD_exobase = n_BC * exp(-vel2/cexo2) / (pi*cexo2)**1.5 * exp(-Iph + ICX) ! dt is negative, so ICX is already negative.
               each_n(iv,iE) = (PSD_exobase + PSD_CX) * dV2(iE,iv)

         enddo
      enddo
       number_density_0D = sum(each_n(:,:))
!       if (cx_calls > 0) then
!          print *, 'ChargeExchange total time (s) =', cx_time_total, ' average per call (s) =', cx_time_total / cx_calls
!       else
!          print *, 'ChargeExchange was not called.'
!       endif
!       stop

      deallocate(each_n)
      deallocate(dV2,solid_angle)
      return
   End



   Subroutine GSE2SPH(pos,finlon,finlat)
   !  Just transform GSE to GEO without considering Earth's rotation. FIX IT when considering the temporal effect of Earth's rotation.

      USE SETTING
      real*8 pos(3)
      real*8 finlon,finlat

      if (pos(2) .gt. 0) then
         if (pos(1) .gt. 0) then
            finlon = atan(pos(2)/pos(1))
         else
            finlon = atan(pos(2)/pos(1)) + pi
         endif
      else
         if (pos(1) .lt. 0) then
            finlon = atan(pos(2)/pos(1)) + pi
         else
            finlon = atan(pos(2)/pos(1)) + 2*pi
         endif
      endif
      finlat = atan(pos(3)/sqrt(pos(1)*pos(1)+pos(2)*pos(2)))

      finlon = finlon * 180.d0/pi
      finlat = finlat * 180.d0/pi

      return
   End

