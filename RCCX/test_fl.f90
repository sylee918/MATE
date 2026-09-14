Program test_fl
   use constants
   use cread1
   use cread2
   use cgrid
   use cfield
   use MATEgrid, only: nRadial, nLon, nLat, drMATE, dangleMATE, dtMATE, rMATE, lonMATE, latMATE
!   implicit none
   integer, parameter :: np = 1000
   real, dimension(np) :: xa1, ya1, za1, rs1, bs1
   real, dimension(np) :: xa2, ya2, za2, rs2, bs2
   real, dimension(je) :: dE
   real, dimension(ig) :: dPA, domega
   real :: dlgE, E1, E2, PA1, PA2
   character(len=80) :: CIMI_flux_file
   integer :: is, it2, i, j, k, iE, iPA, npf1, npf2, iout
   real :: xi, yi, zi, phi1, xlat1, xf, yf, zf, bm_mirror, max_totN

   COMMON /GEOPACK1/ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI, &
          SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33, &
          E11,E21,E31,E12,E22,E32,E13,E23,E33

   is = 1
   CIMI_flux_file = '/home/sylee/CIMI/output/June2008/2008164_MATE_h.fls'

   call readInputData()
   do it2 = 1, 70
      t = real(it2 - 1) * 3600.0
      call calculate_Parmod
      call read_CIMI_flux(CIMI_flux_file, it2, is)
   enddo

   dlgE = log10(gride(is,2)) - log10(gride(is,1))
   do iE = 1, je
      E2 = 10.0**(log10(gride(is,iE)) + dlgE/2.0)
      E1 = 10.0**(log10(gride(is,iE)) - dlgE/2.0)
      dE(iE) = E2 - E1
   enddo
   do iPA = 1, ig
      if (iPA == 1) then; PA1 = 0.0; else; PA1 = asin(gridy(iPA-1)); endif
      if (iPA == ig) then; PA2 = pi/2.0; else; PA2 = asin(gridy(iPA+1)); endif
      dPA(iPA) = (PA2 - PA1) / 2.0
      domega(iPA) = 4.0 * pi * gridy(iPA) * dPA(iPA)
   enddo

   print*, 'Tracing for j=24 (Dusk)...'
   j = 24
   max_totN = 0.0
   do i = 1, iba(j)
      xlat1 = 0.0
      phi1 = xmlto(i,j) * pi / 12.0 + pi
      xi = ro(i,j) * cos(xlat1) * cos(phi1)
      yi = ro(i,j) * cos(xlat1) * sin(phi1)
      zi = 0.0
      call traceF(imod, intB, xi, yi, zi, 1.0, 0.2, 0.0001, 12.0, 0.8, parmod, psi, &
                  np, xf, yf, zf, xa1, ya1, za1, rs1, bs1, npf1, iout)
      call traceF(imod, intB, xi, yi, zi, -1.0, 0.2, 0.0001, 12.0, 0.8, parmod, psi, &
                  np, xf, yf, zf, xa2, ya2, za2, rs2, bs2, npf2, iout)
      
      ! Calculate flux at equator (k=1)
      do iPA = 1, ig
         bm_mirror = (bo(i,j)*1.0e9)/(gridy(iPA)**2)
         if (bs1(1) < bm_mirror) then
            do iE = 1, je
               max_totN = max(max_totN, fl(i,j,iE,iPA))
            enddo
         endif
      enddo
   enddo
   print*, 'Max equatorial differential flux passed mirror check:', max_totN
   print*, 'Field line point counts npf1, npf2:', npf1, npf2
   print*, 'bs1(1) vs bo*1e9:', bs1(1), bo(1, j)*1.0e9
End Program test_fl
