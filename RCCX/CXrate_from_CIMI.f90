Program CXrate_from_CIMI

   use constants
   use cread1
   use cread2
   use cgrid
   use cfield
   use MATEgrid

!   IMPLICIT NONE
   integer, parameter :: np=1000
   integer indx(np),iba0(ip),mstop(4)
   real xk3(np),bm1(np),rm(np),dss(np),dssi(np),yint(np),&
       yint1(np),yinth(np),h3(np),bba(np),zabs(np), &
       rs(np),xs(np),ys(np),zs(np),bs(np),tya3(np)
   real, dimension(np) :: xa1, ya1, za1, rs1, bs1
   real, dimension(np) :: xa2, ya2, za2, rs2, bs2
   character(len=80) :: CIMI_flux_file
   real, dimension(ns,je) :: CXsigma
   real, dimension(ir,ip,ig) :: beta2D
   real, allocatable :: beta1D(:), np1D(:), temp(:,:,:,:)
   real, allocatable :: xGSE(:), yGSE(:), zGSE(:), xGSW(:), yGSW(:), zGSW(:)
   real :: bm_mirror
   real :: r_gse, lat_gse, lon_gse, dr, dlon, dlat, weight, wr, wlon, wlat
   integer :: ir_m, ilon_m, ilat_m, i_off, j_off, k_off, it_m
   real :: wr_w(-1:1), wlon_w(-1:1), wlat_w(-1:1)  ! 3개 포인트 (-1, 0, 1)
   integer :: ii, jj, kk
   real :: d, d2, d3
   real, dimension(ir,ip) :: np2D
   COMMON /GEOPACK1/ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI, &
          SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33, &
          E11,E21,E31,E12,E22,E32,E13,E23,E33

   is=1  ! H+
!   is=2  ! O+

   iRC=0
   iPS=1
      

   if (is==1) then
      CIMI_flux_file="/home/sylee/CIMI/output/June2008/2008164_MATE_h.fls"
   elseif (is==2) then
      CIMI_flux_file="/home/sylee/CIMI/output/June2008/2008164_MATE_o.fls"
   endif

   
   nday=11

   intB=1 ; imod=2 ; ires=1
   itype=1 ; tstart=0. ; dt=1. 
   tmax = 3600 ! 259200
   tint = 3600.

   rb=10
   hlosscone=100.                       ! alt in km
   rc=(re_m+hlosscone*1000.)/re_m       ! losscone in Re 


   dir=1. ; dir1=-1.
   rlim=1.2*rb ; rmin=0.8
!   rlim=3*rb ; rmin=0.1
   err=0.0001 ; dsmax=0.2  

   call readInputData()
   call calculate_Parmod     ! get parmod from CIMI output

print*, intB, imod, ires, itype, tstart, dt, tmax, tint, rb, hlosscone, rc
!stop

   it2=0
!do iday1=2008164,2008164+nday-1
!do iday1=2008164,2008165
do iday1=2008164,2008164+nday-1
   betaMATE = 0.
   nPS_MATE = 0.
   do it=1,nt
   !do it=1,1
      it2=it2+1

      call read_CIMI_flux(CIMI_flux_file,it2,is)
!      if (iday1<2008168) cycle
!      if (it<4) cycle

      print*, 'maxval(fl)', maxval(fl)
      print*, 'maxval(density)', maxval(density)
      call ChargeExchangeCrossSection(CXsigma)
      phi=xmlt*pi/12

      ! Start field line tracing.  Field line tracing from north to south.
      Bfactor=1.
      weight1 = 0.
      jloop: do j=1,ip
   !   jloop: do j=20,24
   !      print*, 'j, iba(j)', j, iba(j)
         iloop: do i=1,iba(j)
   !      iloop: do i=1,ir

   !         xlat1=xlati(i,j)*pi/180.
   !555      phi1=phi(j)+pi                       ! +x corresponing to noon 
   !         xi=rc*cos(xlat1)*cos(phi1)
   !         yi=rc*cos(xlat1)*sin(phi1)
   !         zi=rc*sin(xlat1)

            xlat1 = 0.
   555     phi1 = xmlto(i,j)*pi/12+pi
            xi=ro(i,j)*cos(xlat1)*cos(phi1)
            yi=ro(i,j)*cos(xlat1)*sin(phi1)
            zi=ro(i,j)*sin(xlat1)

            ! Northward
            npf=0
            xa1=0.; ya1=0.; za1=0. ; rs1=0.; bs1=0.; npf1=0.; iout=0
            call traceF(imod,intB,xi,yi,zi,dir,dsmax,err,rlim,rmin,parmod,psi, &
               np,xf,yf,zf,xa1,ya1,za1,rs1,bs1,npf1,iout)

            ! Southward
            xa2=0.; ya2=0.; za2=0. ; rs2=0.; bs2=0.; npf2=0.; iout=0
            call traceF(imod,intB,xi,yi,zi,dir1,dsmax,err,rlim,rmin,parmod,psi, &
               np,xf,yf,zf,xa2,ya2,za2,rs2,bs2,npf2,iout)

            npf=npf1+npf2
   !         npf=npf1+npf2-1


            beta2D=0. 
            dlgE=log10(gride(is,2))-log10(gride(is,1))
            do iE=1, je
               sig=CXsigma(is,iE)
               E2=10**(log10(gride(is,iE))+dlgE/2)
               E1=10**(log10(gride(is,iE))-dlgE/2)
               dE1=E2-E1
               do iPA=1, ig
                  if (fl(i,j,iE,iPA) < 1e-20) cycle
                  if (iPA.eq.1) then; PA1=0; else; PA1=asin(gridy(iPA-1)); endif
                  if (iPA.eq.ig) then; PA2=pi/2; else; PA2=asin(gridy(iPA+1)); endif
                  dPA1=(PA2-PA1)/2.
                  beta2D(i,j,iPA) = beta2D(i,j,iPA) + fl(i,j,iE,iPA)*sig*(2*pi*gridy(iPA)*dE1*dPA1)  ! unit s^-1
   !               beta2D(i,j,iPA) = beta2D(i,j,iPA) + fl(i,j,iE,iPA)*gride(is,iE)*gridy(iPA)*(2*pi*dE1*dPA1)  ! Eflux
               enddo 
            enddo 

            allocate(beta1D(npf), np1D(npf))
            beta1D=0.; np1D=0.
            do iPA=1,ig
               bm_mirror=bo(i,j)*1e9/gridy(iPA)**2
               do k=1,npf1
                  if (bs1(k).lt.bm_mirror) then
                     beta1D(k) = beta1D(k) + beta2D(i,j,iPA)
                     np1D(k) = np1D(k) + density(i,j)
                  endif
               enddo

               do k=1,npf2
                  if (bs2(k).lt.bm_mirror) then
                     beta1D(npf1+k) = beta1D(npf1+k) + beta2D(i,j,iPA)
                     np1D(npf1+k) = np1D(npf1+k) + density(i,j)
                  endif
               enddo
            enddo
!            print*, 'np1D_2', maxval(np1D)

            allocate(xGSE(npf), yGSE(npf), zGSE(npf), xGSW(npf), yGSW(npf), zGSW(npf))

            ! SM to GSE coordinate conversion
            do k=1,npf1
               call smgsw_08(xa1(k),ya1(k),za1(k),xGSW(k),yGSW(k),zGSW(k),1)
               call gswgse_08(xGSW(k),yGSW(k),zGSW(k),xGSE(k),yGSE(k),zGSE(k),1)

               !For test: SM=GSE
   !            xGSE(k) = xa1(k)
   !            yGSE(k) = ya1(k)
   !            zGSE(k) = za1(k)
            enddo

            do k=1,npf2
               call smgsw_08(xa2(k),ya2(k),za2(k),xGSW(k),yGSW(k),zGSW(k),1)
               call gswgse_08(xGSW(k),yGSW(k),zGSW(k),xGSE(k+npf1),yGSE(k+npf1),zGSE(k+npf1),1)

               !For test: SM=GSE
   !            xGSE(npf1+k) = xa2(k)
   !            yGSE(npf1+k) = ya2(k)
   !            zGSE(npf1+k) = za2(k)
            enddo

   ! ==============================================================================
            ! [TSC Interpolation] (Triangular Shaped Cloud)
            ! 특징: 모든 가중치가 양수(Positive)여서 '구멍'이 생기지 않음.
            ! 기준: 가장 가까운 격자점(Nearest Grid Point)을 중심으로 계산
            ! ==============================================================================

            do k=1,npf
               ! 1. 좌표 변환 (Cartesian -> Spherical)
               r_gse = sqrt(xGSE(k)**2 + yGSE(k)**2 + zGSE(k)**2)
               lat_gse = atan2(zGSE(k), sqrt(xGSE(k)**2 + yGSE(k)**2)) * 180.0/pi
               lon_gse = atan2(yGSE(k), xGSE(k)) * 180.0/pi
               if (lon_gse < 0.0) lon_gse = lon_gse + 360.0

               ! 2. 인덱스 찾기: [중요] int(버림) 대신 nint(반올림/가장 가까운 점) 사용
               !    TSC는 가장 가까운 점(Nearest)을 기준으로 -1, 0, 1 범위를 잡아야 함
               
               ! Radial
               ir_m = nint((r_gse - 2.0)/drMATE) + 1  ! 1-based index
               
               ! Longitude
               ilon_m = nint(lon_gse/dangleMATE) + 1
               if (ilon_m < 1) ilon_m = ilon_m + nLon
               if (ilon_m > nLon) ilon_m = ilon_m - nLon
               
               ! Latitude
               ilat_m = nint((lat_gse + 90.0)/dangleMATE) + 1

               ! 3. 거리 계산 (dr, dlon, dlat)
               !    가장 가까운 점(ir_m)으로 부터의 거리이므로 범위는 -0.5 ~ +0.5 가 됨
               !    (rMATE 등 배열 값은 해당 인덱스의 실제 좌표값이어야 함)
               
               ! Check bounds
               if (ir_m >= 1 .and. ir_m <= nRadial .and. &
                  ilon_m >= 1 .and. ilon_m <= nLon .and. &
                  ilat_m >= 1 .and. ilat_m <= nLat) then

                  ! -------------------------------------------------------
                  ! [수정됨] 거리(d) 계산 시 Periodicity 고려
                  ! -------------------------------------------------------
                  
                  ! Radial Distance (비주기적)
                  dr = (r_gse - rMATE(ir_m)) / drMATE
                  
                  ! Longitude Distance (주기적 보정 추가!)
                  ! 단순 차이 계산
                  diff_deg = lon_gse - lonMATE(ilon_m)
                  
                  ! 차이가 180도보다 크면 360도를 빼서 반대쪽 가까운 거리로 인식하게 함
                  ! 예: 359 - 0 = 359 -> -1도로 변환
                  if (diff_deg > 180.0)  diff_deg = diff_deg - 360.0
                  if (diff_deg < -180.0) diff_deg = diff_deg + 360.0
                  
                  dlon = diff_deg / dangleMATE
                  
                  ! Latitude Distance (비주기적)
                  dlat = (lat_gse - latMATE(ilat_m)) / dangleMATE

                  ! -------------------------------------------------------
                  ! TSC Weight Calculation (Always Positive)
                  ! -------------------------------------------------------
                  
                  ! --- Radial Weights ---
                  d = dr  ! -0.5 ~ 0.5
                  d2 = d * d
                  wr_w(0)  = 0.75 - d2               ! Center
                  wr_w(-1) = 0.5 * (0.5 - d)**2      ! Left
                  wr_w(1)  = 0.5 * (0.5 + d)**2      ! Right

                  ! --- Longitude Weights ---
                  d = dlon
                  d2 = d * d
                  wlon_w(0)  = 0.75 - d2
                  wlon_w(-1) = 0.5 * (0.5 - d)**2
                  wlon_w(1)  = 0.5 * (0.5 + d)**2

                  ! --- Latitude Weights ---
                  d = dlat
                  d2 = d * d
                  wlat_w(0)  = 0.75 - d2
                  wlat_w(-1) = 0.5 * (0.5 - d)**2
                  wlat_w(1)  = 0.5 * (0.5 + d)**2

                  ! 4. Accumulate Loop (3x3x3)
                  do i_off = -1, 1
                     ii = ir_m + i_off
                     if (ii < 1 .or. ii > nRadial) cycle

                     do j_off = -1, 1
                        jj = ilon_m + j_off
                        ! Periodic Boundary for Longitude
                        if (jj < 1)    jj = jj + nLon
                        if (jj > nLon) jj = jj - nLon

                        do k_off = -1, 1
                           kk = ilat_m + k_off
                           if (kk < 1 .or. kk > nLat) cycle

                           weight = wr_w(i_off) * wlon_w(j_off) * wlat_w(k_off)
                           
                           betaMATE(ii, jj, kk, it) = betaMATE(ii, jj, kk, it) + beta1D(k) * weight
                           nPS_MATE(ii, jj, kk, it) = nPS_MATE(ii, jj, kk, it) + np1D(k) * weight
                           weight1(ii, jj, kk)      = weight1(ii, jj, kk) + weight
                        enddo
                     enddo
                  enddo
               endif
            enddo ! k - npf

   !         print*, i, j, maxval(beta1D)

            deallocate(xGSE, yGSE, zGSE, xGSW, yGSW, zGSW)
            deallocate(beta1D, np1D)

         enddo iloop
         print*, 'it, j', it, j
      enddo jloop

      ! Normalize betaMATE by weight1
      where (weight1 > 0.0)
         betaMATE(:,:,:,it) = betaMATE(:,:,:,it) / weight1(:,:,:)
         nPS_MATE(:,:,:,it) = nPS_MATE(:,:,:,it) / weight1(:,:,:)
      end where

   !if (iday1==2008166 .and. it==19) exit
   !if (iday1==2008164 .and. it==2) exit
   !if (iday1==2008164 .and. it==4) exit

   enddo !it

   call write_4D(betaMATE,iday1,is,'RCCX')
   call write_4D(nPS_MATE,iday1,is,'nPS')
!   stop
!if (iday1==2008166 .and. it==19) stop
!if (iday1==2008164 .and. it==2) stop

enddo !iday

   print*, 'betaMATE normalization complete'
   print*, 'Sum of betaMATE:', sum(betaMATE), maxval(betaMATE)
   print*, 'Sum of weight1:', sum(weight1), maxval(weight1)


End Program CXrate_from_CIMI

