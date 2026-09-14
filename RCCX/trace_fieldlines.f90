Program trace_fieldlines

   use constants
   use cread1
   use cread2
   use cgrid
   use cfield

!   implicit none

   integer, parameter :: np = 10000
   integer :: m, npf1, npf2, npf3, npf4, iout1, iout2, iout3, iout4
   real :: ut, phi_mid, phi_dusk
   real :: xi1, yi1, zi1, xf1, yf1, zf1
   real :: xi2, yi2, zi2, xf2, yf2, zf2
   real :: xi3, yi3, zi3, xf3, yf3, zf3
   real :: xi4, yi4, zi4, xf4, yf4, zf4
   real, dimension(np) :: xa1, ya1, za1, ra1, ba1
   real, dimension(np) :: xa2, ya2, za2, ra2, ba2
   real, dimension(np) :: xa3, ya3, za3, ra3, ba3
   real, dimension(np) :: xa4, ya4, za4, ra4, ba4
   character(len=120) :: CIMI_flux_file

   COMMON /GEOPACK1/ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI, &
          SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33, &
          E11,E21,E31,E12,E22,E32,E13,E23,E33

   is = 1
   CIMI_flux_file = "/home/sylee/CIMI/output/June2008/2008164_MATE_h.fls"

   intB = 1 ; imod = 2 ; ires = 1
   itype = 1 ; tstart = 0. ; dt = 1. 
   tmax = 3600 ; tint = 3600.
   rlim = 120.0 ; rmn = 0.999
   err = 0.0001 ; dsmax = 0.05

   call readInputData()

   ! Advance to it = 70
   do it2 = 1, 70
      t = real(it2 - 1) * 3600.0
      call calculate_Parmod
      call read_CIMI_flux(CIMI_flux_file, it2, is)
   end do

   write(*,*) 'Field Line Tracing for it = 70, iday =', iday2, 'ihour =', ihour
   write(*,*) 'Dipole Tilt PSI (rad) =', psi, 'parmod(1:4) =', parmod(1:4)

   ! Recalculate Geopack coordinate transformations
   call recalc_08(2008, iday2, ihour, imin, isec, -400.0, 0.0, 0.0)

   phi_mid  = 175.0 * pi / 180.0  ! Lon = 175 deg
   phi_dusk = 90.0 * pi / 180.0   ! 18 LT (Dusk, 90 deg)

   ! Line 1: Lon = 175 deg, MLAT = 60 deg
   xi1 = 1.001 * cos(60.0 * pi / 180.0) * cos(phi_mid)
   yi1 = 1.001 * cos(60.0 * pi / 180.0) * sin(phi_mid)
   zi1 = 1.001 * sin(60.0 * pi / 180.0)
   call traceF(imod, intB, xi1, yi1, zi1, 1.0, dsmax, err, rlim, rmn, parmod, psi, &
               np, xf1, yf1, zf1, xa1, ya1, za1, ra1, ba1, npf1, iout1)
   write(*,*) 'Line 1 (Lon 175, MLAT 60): Traced', npf1, 'pts, R_max =', maxval(ra1(1:npf1))

   ! Line 2: Lon = 175 deg, MLAT = 65 deg
   xi2 = 1.001 * cos(65.0 * pi / 180.0) * cos(phi_mid)
   yi2 = 1.001 * cos(65.0 * pi / 180.0) * sin(phi_mid)
   zi2 = 1.001 * sin(65.0 * pi / 180.0)
   call traceF(imod, intB, xi2, yi2, zi2, 1.0, dsmax, err, rlim, rmn, parmod, psi, &
               np, xf2, yf2, zf2, xa2, ya2, za2, ra2, ba2, npf2, iout2)
   write(*,*) 'Line 2 (Lon 175, MLAT 65): Traced', npf2, 'pts, R_max =', maxval(ra2(1:npf2))

   ! Line 3: Dusk (90 deg, 18 LT), MLAT = 60 deg
   xi3 = 1.001 * cos(60.0 * pi / 180.0) * cos(phi_dusk)
   yi3 = 1.001 * cos(60.0 * pi / 180.0) * sin(phi_dusk)
   zi3 = 1.001 * sin(60.0 * pi / 180.0)
   call traceF(imod, intB, xi3, yi3, zi3, 1.0, dsmax, err, rlim, rmn, parmod, psi, &
               np, xf3, yf3, zf3, xa3, ya3, za3, ra3, ba3, npf3, iout3)
   write(*,*) 'Line 3 (Dusk 18LT, MLAT 60): Traced', npf3, 'pts, R_max =', maxval(ra3(1:npf3))

   ! Line 4: Dusk (90 deg, 18 LT), MLAT = 65 deg
   xi4 = 1.001 * cos(65.0 * pi / 180.0) * cos(phi_dusk)
   yi4 = 1.001 * cos(65.0 * pi / 180.0) * sin(phi_dusk)
   zi4 = 1.001 * sin(65.0 * pi / 180.0)
   call traceF(imod, intB, xi4, yi4, zi4, 1.0, dsmax, err, rlim, rmn, parmod, psi, &
               np, xf4, yf4, zf4, xa4, ya4, za4, ra4, ba4, npf4, iout4)
   write(*,*) 'Line 4 (Dusk 18LT, MLAT 65): Traced', npf4, 'pts, R_max =', maxval(ra4(1:npf4))

   ! Save all 4 traced lines to ASCII file
   open(unit=50, file='fieldlines_SM_it070.dat', status='replace')
   write(50, '(4I8)') npf1, npf2, npf3, npf4
   write(50, '(A)') '# Line 1: Lon 175 deg, MLAT 60 deg'
   do m = 1, npf1
      write(50, '(5E16.7)') xa1(m), ya1(m), za1(m), ra1(m), ba1(m)
   end do
   write(50, '(A)') '# Line 2: Lon 175 deg, MLAT 65 deg'
   do m = 1, npf2
      write(50, '(5E16.7)') xa2(m), ya2(m), za2(m), ra2(m), ba2(m)
   end do
   write(50, '(A)') '# Line 3: Dusk 18 LT (90 deg), MLAT 60 deg'
   do m = 1, npf3
      write(50, '(5E16.7)') xa3(m), ya3(m), za3(m), ra3(m), ba3(m)
   end do
   write(50, '(A)') '# Line 4: Dusk 18 LT (90 deg), MLAT 65 deg'
   do m = 1, npf4
      write(50, '(5E16.7)') xa4(m), ya4(m), za4(m), ra4(m), ba4(m)
   end do
   close(50)

   write(*,*) 'Successfully saved fieldlines_SM_it070.dat for Lon 175 deg & Dusk (18LT)!'

End Program trace_fieldlines
