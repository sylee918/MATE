      Subroutine Initialize(input_dir, init, fin, flags, tag, thread_num)

         include "constants.inc"
         external Init_Parameter, read_ind_binary, read_fin_binary

         real*8, dimension(N_vel_directions,nRadial,nEnergy,7) :: init, fin
         integer, dimension(N_vel_directions,nRadial,nEnergy) :: flags

         real*8 radial_distance_range(nRadial), energy_range(nEnergy)
         real*8 longitude_range(nLong), latitude_range(nLat), latitudeNS_range(nLat_NS)
         real*8 radial_boundary(2), tmax
         character*70 input_dir
         character*30 tag
         integer thread_num

         call Init_Parameter(radial_distance_range, energy_range, longitude_range, latitude_range, latitudeNS_range, radial_boundary,tmax)
         init=0.d0
  
         print*, 'thread_num in Initialize', thread_num
         call read_ind_binary(flags, input_dir, tag, thread_num)
         call read_fin_binary(  fin, input_dir, tag, thread_num)

         return
      End


      Subroutine Calculate_Density(fin, flags, current_time, nH_BC, TH_BC, number_density_1D, bph, rank)
         ! "cdensity" in python code
!         use omp_lib
         use, intrinsic :: ieee_arithmetic
         include "constants.inc"
         external calculate_Velocity_Volume_Element
         external GSE2SPH

         real*8, dimension(N_vel_directions,nRadial,nEnergy,7) :: fin
         integer, dimension(N_vel_directions,nRadial,nEnergy) :: flags
         real*8, dimension(nEnergy,N_vel_directions) :: dV2
         real*8, dimension(start_ydoy_index:end_ydoy_index) :: bph
         real*8, dimension(:,:,:), allocatable :: each_n
         real*8 pos(3), vel(3), vel2
         real*8 temp_BC, n_BC, vel_BC(3), fac, number_density
         real*8 number_density_1D(nRadial), cexo2
         integer iR,iE,iv, i
         real*8, dimension(nbx,nby,nbtperday,start_ydoy-nt_bwd_bc:end_ydoy) :: nH_BC, TH_BC
         real*8 finlon, finlat
         real*8 current_time, t0, t1, Iph
         integer iflon, iflat, it, rank, quotient
         character*30 fn2D, fn3D
         integer idoy, iday

         vel_BC = 0.d0;
         call calculate_Velocity_Volume_Element(dV2)

         allocate(each_n(N_vel_directions,nRadial,nEnergy))
         each_n = 0.d0
         number_density_1D = 0.d0

         fac = 2.d0*kb/mH

         do iR=1,nRadial      ! Outermost iR-loop
            do iE=1,nEnergy
               do iv=1,N_vel_directions
                  t0 = current_time + fin(iv,iR,iE,1)/86400    ! unit day
                  idoy = int(t0)                               ! yyyy+doy
                  t1 = (t0 - idoy)*86400                       ! hms in seconds
                  it = floor(t1/tb_res)+1
                  if (idoy .lt. start_ydoy-nt_bwd_bc) then ; idoy=start_ydoy-nt_bwd_bc ; it=1 ; endif
                  if (flags(iv,iR,iE) .eq. 1) then
                     do i=1,3
                        pos(i) = fin(iv,iR,iE,i+1)
                        vel(i) = fin(iv,iR,iE,i+4)
                     enddo

                     call GSE2SPH(pos,finlon,finlat)
!                     iflon=floor(finlon/bc_res)+1                !   0 < lon < 360
!                    ** It is due to the longitude is defined from -180 to 180 in python, not 0 to 360.
!                    ** If it is defined from 0 to 360, then use the above one.
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

                     !! FIX ME (05.28.2024) : Temporarily fixed the 180 degree
                     !difference of MSIS longitude
!                     iflon = ilon + 180/bc_res
!                     if (iflon .gt. 360/bc_res) then
!                        iflon = iflon - 360/bc_res
!                     endif

                     n_BC    = nH_BC(iflon,iflat,it,idoy)
                     temp_BC = TH_BC(iflon,iflat,it,idoy)

                     !! ** FIX ME (above): Trilinear interpolation is desired for more accurate calculation.
                     !!                    Current code is just the 0th-order interpolation.

                     if (idoy .eq. int(current_time)) then
                        Iph = bph(idoy) * abs(fin(iv,iR,iE,1))
                     else
                        Iph = bph(idoy) * (86400-t1)
                        do iday=idoy+1,int(current_time)-1
                           Iph = Iph + bph(iday)*86400
                        enddo
                        Iph = Iph + bph(iday) * (current_time-int(current_time))
                     endif

!                     vel = (vel - vel_BC)
                     cexo2 = fac*temp_BC
                     vel2 = sum(vel*vel)
                     number_density = n_BC * exp(-vel2/cexo2) / (pi*cexo2)**1.5 * exp(-Iph)
                     each_n(iv,iR,iE) = number_density * dV2(iE,iv) !* dV1(iR)

                  endif
               enddo
            enddo
            number_density_1D(iR) = sum(each_n(:,iR,:))
         enddo
         deallocate(each_n)

         print*, number_density_1D

         return
      End



      Subroutine MSIS_averaged_over_exobase(MSIS_nH,MSIS_TH)

         include "constants.inc"
!         real*8, dimension(nRadial) :: MSIS_nH, MSIS_TH
         real*8, dimension(41) :: MSIS_nH, MSIS_TH

         MSIS_nH = [647815.75,644301.19,640786.75,637272.19,633479.62,624297.38,610958.62,593882.81,573313.12,549425.44,522383.88,492333.12,459380.84, &
          423633.91,385176.69,344093.22,300450.25,254312.52,207411.00,174808.81,143747.97,120508.17,100757.82,83995.781,71875.648,61286.113,52548.836, &
          46214.402,40904.266,36620.562,33277.266,30971.016,29303.338,28260.053,27836.275,27868.945,27901.617,27934.287,27966.957,27999.629,28032.299]
         MSIS_TH = [431.80115,434.83173,437.86230,440.89291,444.16327,452.08124,463.58334,478.30795,496.04541,516.64398,539.96216,565.87524,594.29028, &
          625.11517,658.27716,693.70386,731.33759,771.12256,812.71167,853.42889,895.14233,936.03601,976.58215,1016.9861,1055.1702,1092.6792,1129.0428, &
          1162.3274,1193.8517,1223.2327,1250.0793,1272.8732,1292.5026,1308.6052,1320.7113,1329.6329,1338.5547,1347.4763,1356.3981,1365.3197,1374.2415]

         return
      End



      Subroutine GSE2SPH(pos,finlon,finlat)
      !  Just transform GSE to GEO without considering Earth's rotation. FIX IT when considering the temporal effect of Earth's rotation.

         include "constants.inc"
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

