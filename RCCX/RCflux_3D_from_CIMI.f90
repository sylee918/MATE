Program RCflux_3D_from_CIMI

   use constants
   use cread1
   use cread2
   use cgrid
   use cfield

!   implicit none

   ! High-Resolution Grid Specifications:
   ! Radial: 1.0 to 10.0 Re with dr = 0.1 Re (91 points)
   ! Longitude: 0 to 359 deg with dlon = 1.0 deg (360 points)
   ! Latitude: -90 to +90 deg with dlat = 1.0 deg (181 points)
   integer, parameter :: nRadial = 91
   integer, parameter :: nLon    = 360
   integer, parameter :: nLat    = 181
   real, parameter    :: drGrid = 0.1, dangleGrid = 1.0

   real, dimension(nRadial) :: rGrid
   real, dimension(nLon)    :: lonGrid
   real, dimension(nLat)    :: latGrid

   integer, parameter :: np = 2500
   integer :: indx(np), iba0(ip), mstop(4)
   real :: xk3(np), bm1(np), rm(np), dss(np), dssi(np), yint(np), &
           yint1(np), yinth(np), h3(np), bba(np), zabs(np), &
           rs(np), xs(np), ys(np), zs(np), bs(np), tya3(np)
   real, dimension(np) :: xa1, ya1, za1, rs1, bs1
   real, dimension(np) :: xa2, ya2, za2, rs2, bs2
   real, dimension(np) :: xa, ya, za, bs_all
   character(len=120)  :: CIMI_flux_file
   character(len=120)  :: out_file

   real, dimension(je) :: dE
   real, dimension(ig) :: dPA, domega
   real, dimension(je, ig) :: diff_flux_eq

   ! 3D and 4D Grid Accumulators (SM Coordinates)
   real, allocatable :: tot_Nflux_3D(:,:,:), tot_Eflux_3D(:,:,:), weight_3D(:,:,:)
   real, allocatable :: omni_flux_4D(:,:,:,:)

   real, allocatable :: totN_1D(:), totE_1D(:)
   real, allocatable :: omni_1D(:,:), fl_1D(:,:,:)
   real :: bm_mirror, Bfactor
   real :: r_sm, lat_sm, lon_sm, dr, dlon, dlat, weight, diff_deg
   integer :: ir_m, ilon_m, ilat_m, i_off, j_off, k_off
   real :: wr_w(-1:1), wlon_w(-1:1), wlat_w(-1:1)
   integer :: ii, jj, kk, npf, npf1, npf2, iout
   real :: d, d2, dlgE, E1, E2, PA1, PA2
   real :: xi, yi, zi, phi1, xlat1, xf, yf, zf
   integer :: i, j, k, iE, iPA, it2, it_target

   COMMON /GEOPACK1/ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI, &
          SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33, &
          E11,E21,E31,E12,E22,E32,E13,E23,E33

   is = 1  ! H+ (Proton)
   it_target = 70

   CIMI_flux_file = "/home/sylee/CIMI/output/June2008/2008164_MATE_h.fls"
   out_file = "RC_3D_flux_SM_it070_highres.dat"

   intB = 1 ; imod = 2 ; ires = 1
   itype = 1 ; tstart = 0. ; dt = 1. 
   tmax = 3600 ; tint = 3600.

   rb = 10.0
   hlosscone = 100.
   rc = (re_m + hlosscone*1000.) / re_m

   dir = 1.0 ; dir1 = -1.0
   rlim = 1.2 * rb ; rmin = 0.8
   err = 0.0001 ; dsmax = 0.05

   print*, '==================================================================='
   print*, '  CIMI High-Resolution 3D Ring Current Calculator (Geopack TS04)'
   print*, '  Target Time: it = 70 | Coordinate System: SM'
   print*, '  Grid: R(1.0-10.0 Re, dr=0.1) x Lon(0-360, dlon=1 deg) x Lat(-90..90, dlat=1 deg)'
   print*, '  Grid Dimensions:', nRadial, 'x', nLon, 'x', nLat, '=', nRadial*nLon*nLat, 'cells'
   print*, '==================================================================='

   ! Initialize High-Resolution Grid Coordinates
   do i = 1, nRadial
      rGrid(i) = 1.0 + real(i - 1) * drGrid
   enddo
   do j = 1, nLon
      lonGrid(j) = real(j - 1) * dangleGrid
   enddo
   do k = 1, nLat
      latGrid(k) = -90.0 + real(k - 1) * dangleGrid
   enddo

   ! Allocate High-Resolution Arrays
   allocate(tot_Nflux_3D(nRadial, nLon, nLat))
   allocate(tot_Eflux_3D(nRadial, nLon, nLat))
   allocate(weight_3D(nRadial, nLon, nLat))
   allocate(omni_flux_4D(nRadial, nLon, nLat, je))

   tot_Nflux_3D = 0.0
   tot_Eflux_3D = 0.0
   omni_flux_4D = 0.0
   weight_3D    = 0.0

   call readInputData()

   ! Read CIMI flux up to it_target = 70
   print*, 'Loading CIMI flux time steps 1 to 70...'
   do it2 = 1, it_target
      t = real(it2 - 1) * 3600.0
      call calculate_Parmod
      call read_CIMI_flux(CIMI_flux_file, it2, is)
   enddo

   print*, '==================================================================='
   print*, '--> Target Time Step it2 = 70 Successfully Loaded!'
   print*, '    Equatorial Max Differential Flux:', maxval(fl)
   print*, '    Equatorial Max B-field (nT):', maxval(bo)*1.0e9
   print*, '==================================================================='

   ! Compute Energy Integration Weights (dE in keV)
   dlgE = log10(gride(is,2)) - log10(gride(is,1))
   do iE = 1, je
      E2 = 10.0**(log10(gride(is,iE)) + dlgE/2.0)
      E1 = 10.0**(log10(gride(is,iE)) - dlgE/2.0)
      dE(iE) = E2 - E1
   enddo

   ! Compute Pitch Angle Integration Weights (dPA in radians, solid angle domega)
   do iPA = 1, ig
      if (iPA == 1) then
         PA1 = 0.0
      else
         PA1 = asin(gridy(iPA-1))
      endif
      if (iPA == ig) then
         PA2 = pi / 2.0
      else
         PA2 = asin(gridy(iPA+1))
      endif
      dPA(iPA) = (PA2 - PA1) / 2.0
      domega(iPA) = 4.0 * pi * gridy(iPA) * dPA(iPA)
   enddo

   ! Start 3D Magnetic Field Line Tracing in SM Coordinates
   print*, 'Starting High-Resolution Geopack TS04 Field Line Tracing & TSC Deposition...'
   jloop: do j = 1, ip
      iloop: do i = 1, iba(j)
         
         xlat1 = 0.0
         phi1 = xmlto(i,j) * pi / 12.0 + pi  ! +X corresponds to Noon in SM
         xi = ro(i,j) * cos(xlat1) * cos(phi1)
         yi = ro(i,j) * cos(xlat1) * sin(phi1)
         zi = ro(i,j) * sin(xlat1)

         ! Northward Field Line Tracing (fine step dsmax = 0.05 Re)
         xa1 = 0.0; ya1 = 0.0; za1 = 0.0; rs1 = 0.0; bs1 = 0.0; npf1 = 0; iout = 0
         call traceF(imod, intB, xi, yi, zi, dir, dsmax, err, rlim, rmin, parmod, psi, &
                     np, xf, yf, zf, xa1, ya1, za1, rs1, bs1, npf1, iout)

         ! Southward Field Line Tracing (fine step dsmax = 0.05 Re)
         xa2 = 0.0; ya2 = 0.0; za2 = 0.0; rs2 = 0.0; bs2 = 0.0; npf2 = 0; iout = 0
         call traceF(imod, intB, xi, yi, zi, dir1, dsmax, err, rlim, rmin, parmod, psi, &
                     np, xf, yf, zf, xa2, ya2, za2, rs2, bs2, npf2, iout)

         npf = npf1 + npf2
         if (npf < 1) cycle

         ! Combine Northward and Southward field line points
         do k = 1, npf1
            xa(k) = xa1(k); ya(k) = ya1(k); za(k) = za1(k); bs_all(k) = bs1(k)
         enddo
         do k = 1, npf2
            xa(npf1 + k) = xa2(k); ya(npf1 + k) = ya2(k); za(npf1 + k) = za2(k); bs_all(npf1 + k) = bs2(k)
         enddo

         ! Allocate 1D field line flux arrays
         allocate(totN_1D(npf), totE_1D(npf), omni_1D(npf, je), fl_1D(npf, je, ig))
         totN_1D = 0.0
         totE_1D = 0.0
         omni_1D = 0.0
         fl_1D   = 0.0

         ! Map Equatorial Flux along the Field Line using Liouville & Mirror Condition
         do iPA = 1, ig
            bm_mirror = (bo(i,j) * 1.0e9) / (gridy(iPA)**2)
            do k = 1, npf
               if (bs_all(k) < bm_mirror) then
                  do iE = 1, je
                     fl_1D(k, iE, iPA) = fl(i, j, iE, iPA)
                     omni_1D(k, iE) = omni_1D(k, iE) + fl(i, j, iE, iPA) * domega(iPA)
                  enddo
               endif
            enddo
         enddo

         ! Compute Total Number Flux and Total Energy Flux along Field Line
         do k = 1, npf
            do iE = 1, je
               totN_1D(k) = totN_1D(k) + omni_1D(k, iE) * dE(iE)
               totE_1D(k) = totE_1D(k) + omni_1D(k, iE) * (dE(iE) * gride(is, iE))
            enddo
         enddo

         ! ===============================================================
         ! [TSC Interpolation] Map Field Line Points to High-Res 3D SM Grid
         ! ===============================================================
         do k = 1, npf
            ! 1. Spherical SM Coordinates (Direct SM, no GSE conversion)
            r_sm   = sqrt(xa(k)**2 + ya(k)**2 + za(k)**2)
            lat_sm = atan2(za(k), sqrt(xa(k)**2 + ya(k)**2)) * 180.0 / pi
            lon_sm = atan2(ya(k), xa(k)) * 180.0 / pi
            if (lon_sm < 0.0) lon_sm = lon_sm + 360.0

            ! 2. High-Res Grid Index
            ir_m   = nint((r_sm - 1.0) / drGrid) + 1
            ilon_m = nint(lon_sm / dangleGrid) + 1
            if (ilon_m < 1) ilon_m = ilon_m + nLon
            if (ilon_m > nLon) ilon_m = ilon_m - nLon
            ilat_m = nint((lat_sm + 90.0) / dangleGrid) + 1

            ! 3. Bounds Check & Distance Weights
            if (ir_m >= 1 .and. ir_m <= nRadial .and. &
                ilon_m >= 1 .and. ilon_m <= nLon .and. &
                ilat_m >= 1 .and. ilat_m <= nLat) then

               dr = (r_sm - rGrid(ir_m)) / drGrid

               diff_deg = lon_sm - lonGrid(ilon_m)
               if (diff_deg > 180.0)  diff_deg = diff_deg - 360.0
               if (diff_deg < -180.0) diff_deg = diff_deg + 360.0
               dlon = diff_deg / dangleGrid

               dlat = (lat_sm - latGrid(ilat_m)) / dangleGrid

               ! TSC Weight Calculation (Smooth 3x3x3 cloud)
               d = dr
               d2 = d * d
               wr_w(0)  = 0.75 - d2
               wr_w(-1) = 0.5 * (0.5 - d)**2
               wr_w(1)  = 0.5 * (0.5 + d)**2

               d = dlon
               d2 = d * d
               wlon_w(0)  = 0.75 - d2
               wlon_w(-1) = 0.5 * (0.5 - d)**2
               wlon_w(1)  = 0.5 * (0.5 + d)**2

               d = dlat
               d2 = d * d
               wlat_w(0)  = 0.75 - d2
               wlat_w(-1) = 0.5 * (0.5 - d)**2
               wlat_w(1)  = 0.5 * (0.5 + d)**2

               ! 4. Accumulate into 3D and 4D High-Res Grids
               do i_off = -1, 1
                  ii = ir_m + i_off
                  if (ii < 1 .or. ii > nRadial) cycle

                  do j_off = -1, 1
                     jj = ilon_m + j_off
                     if (jj < 1)    jj = jj + nLon
                     if (jj > nLon) jj = jj - nLon

                     do k_off = -1, 1
                        kk = ilat_m + k_off
                        if (kk < 1 .or. kk > nLat) cycle

                        weight = wr_w(i_off) * wlon_w(j_off) * wlat_w(k_off)

                        tot_Nflux_3D(ii, jj, kk) = tot_Nflux_3D(ii, jj, kk) + totN_1D(k) * weight
                        tot_Eflux_3D(ii, jj, kk) = tot_Eflux_3D(ii, jj, kk) + totE_1D(k) * weight
                        weight_3D(ii, jj, kk)    = weight_3D(ii, jj, kk) + weight

                        do iE = 1, je
                           omni_flux_4D(ii, jj, kk, iE) = omni_flux_4D(ii, jj, kk, iE) + omni_1D(k, iE) * weight
                        enddo
                     enddo
                  enddo
               enddo
            endif
         enddo ! k loop

         deallocate(totN_1D, totE_1D, omni_1D, fl_1D)
      enddo iloop
      if (mod(j, 6) == 0) then
         write(*, '(A, I3, A, I3)') '  Completed MLT sector j = ', j, ' / ', ip
      endif
   enddo jloop

   ! Normalize by TSC Weights
   do kk = 1, nLat
      do jj = 1, nLon
         do ii = 1, nRadial
            if (weight_3D(ii, jj, kk) > 0.0) then
               tot_Nflux_3D(ii, jj, kk) = tot_Nflux_3D(ii, jj, kk) / weight_3D(ii, jj, kk)
               tot_Eflux_3D(ii, jj, kk) = tot_Eflux_3D(ii, jj, kk) / weight_3D(ii, jj, kk)
               do iE = 1, je
                  omni_flux_4D(ii, jj, kk, iE) = omni_flux_4D(ii, jj, kk, iE) / weight_3D(ii, jj, kk)
               enddo
            endif
         enddo
      enddo
   enddo

   print*, '==================================================================='
   print*, '  High-Resolution 3D Field Line Tracing & TSC Deposition Complete!'
   print*, '  Max 3D Total Number Flux: ', maxval(tot_Nflux_3D), ' cm^-2 s^-1'
   print*, '  Max 3D Total Energy Flux: ', maxval(tot_Eflux_3D), ' keV cm^-2 s^-1'
   print*, '==================================================================='

   ! Write High-Resolution Output Binary Data File
   open(unit=50, file=trim(out_file), form='unformatted', status='replace', access='stream')
   
   ! Header & Grid Meta
   write(50) nRadial, nLon, nLat, je, ig, it_target
   write(50) rGrid
   write(50) lonGrid
   write(50) latGrid
   write(50) (gride(is, iE), iE=1,je)
   write(50) (gridy(iPA), iPA=1,ig)
   
   ! 3D and 4D Datasets
   write(50) tot_Nflux_3D
   write(50) tot_Eflux_3D
   write(50) omni_flux_4D
   
   close(50)
   print*, '--> High-Res Output saved to: ', trim(out_file)

   deallocate(tot_Nflux_3D, tot_Eflux_3D, weight_3D, omni_flux_4D)

End Program RCflux_3D_from_CIMI
