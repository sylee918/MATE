   Subroutine Calculate_Density(fin, flags, current_time, number_density_0D)
      ! "cdensity" in python code
!         use omp_lib
      USE SETTING
      USE EXOBASE_BC, only: nH_BC, TH_BC
      USE SOLAR_LYMAN_ALPHA, only: bph
      USE MPI_MATE, only: nR_loc
      USE VOLUME_ELEMENT
      USE ChargeExchange
      use, intrinsic :: ieee_arithmetic
      external GSE2SPH

      real*8, dimension(nvel,nEnergy,7), intent(in) :: fin
      integer, dimension(nvel,nEnergy), intent(in) :: flags
      real*8, intent(in) :: current_time
      real*8, intent(out) :: number_density_0D

      real*8, dimension(:,:), allocatable :: each_n
      real*8 pos(3), vel(3), vel2
      real*8 temp_BC, n_BC, vel_BC(3), fac, number_density
      integer iE,iv, i
      real*8 finlon, finlat, cexo2
      real*8 t0, t1, Iph, ICX
      integer iflon, iflat, it, quotient
      integer idoy, iday


      vel_BC = 0.d0;
      allocate(dV2(nEnergy,nvel), solid_angle(nvel))
      call calculate_Velocity_Volume_Element(dV2)

      allocate(each_n(nvel,nEnergy))
      each_n = 0.d0
      number_density_0D = 0.d0

      fac = 2.d0*kb/mH

      do iE=1,nEnergy
         do iv=1,nvel
            t0 = current_time + fin(iv,iE,1)/86400.    ! unit day
            idoy = int(t0)                               ! yyyy+doy
            t1 = (t0 - idoy)*86400.                       ! hms in seconds
            it = floor(t1/tb_res)+1
            if (idoy .lt. start_ydoy-nt_bwd_bc) then ; idoy=start_ydoy-nt_bwd_bc ; it=1 ; endif
            if (flags(iv,iE) .eq. 1) then
               do i=1,3
                  pos(i) = fin(iv,iE,i+1)
                  vel(i) = fin(iv,iE,i+4)
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

               n_BC    = nH_BC(iflon,iflat,it,idoy)
               temp_BC = TH_BC(iflon,iflat,it,idoy)

               !! ** FIX ME (above): Trilinear interpolation is desired for more accurate calculation.
               !!                    Current code is just the 0th-order interpolation.

               if (i_Photoionization .eq. 1) then
                  if (idoy .eq. int(current_time)) then
                     Iph = bph(idoy) * abs(fin(iv,iE,1))
                  else
                     Iph = bph(idoy) * (86400.-t1)
                     do iday=idoy+1,int(current_time)-1
                        Iph = Iph + bph(iday)*86400.
                     enddo
                     Iph = Iph + bph(iday) * (current_time-int(current_time))*86400.
                  endif
               endif

               call Calculate_ChargeExchange(iE,iv, current_time, ICX)

               !vel = (vel - vel_BC)
               cexo2 = fac*temp_BC
               vel2 = sum(vel*vel)
               number_density = n_BC * exp(-vel2/cexo2) / (pi*cexo2)**1.5 * exp(-Iph) * exp(-ICX)
               each_n(iv,iE) = number_density * dV2(iE,iv)

            endif
         enddo
      enddo
      number_density_0D = sum(each_n(:,:))

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

