MODULE SOLAR_LYMAN_ALPHA

   USE SETTING
   IMPLICIT NONE

   real*8, dimension(start_ydoy_index:end_ydoy_index) :: Lya, bph

contains

   Subroutine read_Lya_Bph

      IMPLICIT NONE

      character(len=80) :: line
      integer :: year, doy, i, yyyydoy, IO_unit
      real :: f10_7, f107a, ap, lyman_alpha, beta_ph, factor

      open(newunit=IO_unit, file=trim(Lya_dir), status='old', action='read')
      read(IO_unit, '(A)', iostat=i) line   ! Skip the header line
      do while (.true.)
         read(IO_unit, '(A)', iostat=i) line
         if (i /= 0) exit
         read(line, *, iostat=i) year, doy, f10_7, f107a, ap, lyman_alpha, beta_ph
         yyyydoy = year * 1000 + doy
         if (yyyydoy >= start_ydoy_index .and. yyyydoy <= end_ydoy_index) then
!               print *, "Year:", year, "DOY:", doy, "Lyman-alpha:", lyman_alpha
            Lya(yyyydoy) = lyman_alpha
            bph(yyyydoy) = beta_ph
         end if
      end do
      close(IO_unit)

      ! Convert line-integrated Lya to line-centered Lya [Emerich et al., 2005]
      factor = (h*c/121.6d-9)*1e11*1e4       !  121.6e-9 m for the wavelength of Lyman-alpha, 1e4 for m2->cm2, and 1e12 from Emmerich et al. (2005)
      Lya = 0.64*(Lya/factor)**1.21        

      return
   End
   
END MODULE SOLAR_LYMAN_ALPHA



