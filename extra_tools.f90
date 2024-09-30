      Subroutine gen_points(vel_dir)
         ! Generates points on a sphere at n values of theta.
         ! Points are roughly evenly spaced on the sphere \
         ! n : number of theta values for points, odd n are preferred for most even spacing \
         !  output : boolean; if 0, returns points as normal;
         !           if 1, returns array specifying how many points are at each theta value"

         use Module_for_NVelocityDirection
         include "Setting.inc"
         integer :: nsize, n2
         integer :: i, j, i0
         real*8 piset(nTheta), piset2(0:nTheta), piset2_i
         real*8 :: co1, co2, co3, co4
         real*8, dimension(N_vel_directions,3) :: vel_dir         
         real*8, dimension(:,:), allocatable :: coords
         real*8 :: t1,t2,t3

         do i=1,nTheta
            piset(i)=i*1.d0/(nTheta-1)*pi
            piset2(i)=i*2.d0/(nTheta-1)*pi
         enddo
         piset2(0)=0.d0;

         nsize=2
         do j=1,nTheta-2
            n2=nint((nTheta*2-2)*dsin(piset(j)))
            do i=0,n2-1
               nsize=nsize+1
            enddo
         enddo

         if (nsize .ne. N_vel_directions) then
            print*, "gen_points: <N_vel_directions> is not equal to <nsize>"
            print*, "nsize = ", nsize
            stop
         endif

         allocate(coords(nsize,3)); coords=0.d0

         i0=1+1 ! index for fortran
         do j=1,nTheta-2
            if (piset(j).eq.pi/2) then
               co1=1.d0; co2=0.d0
            else
               co1=dsin(piset(j)) ; co2=dcos(piset(j))
            endif

            n2=nint((nTheta*2-2)*co1)
            do i=0,n2-1
               piset2_i = 2.d0*i/n2*pi

               if (piset2_i.eq.pi/2 .or. piset2_i.eq.pi .or. piset2_i.eq.1.5*pi) then
                  co3=nint(dsin(piset2_i)); co4=nint(dcos(piset2_i))
               else
                  co3=dsin(piset2_i); co4=dcos(piset2_i)
               endif
               coords(i0,1) = co4*co1
               coords(i0,2) = co3*co1
               coords(i0,3) = co2

               i0=i0+1
         enddo
         enddo
         coords(1,3)=1.d0
         coords(nsize,3)=-1.d0

         t1=0.d0; t2=0.d0; t3=0.d0
         do i=1,nsize
            t1=t1+coords(i,1)
            t2=t2+coords(i,2)
            t3=t3+coords(i,3)
         enddo

         vel_dir = coords

         deallocate(coords)
  
         return
      End


      Subroutine gen_points_for_each_row(row)
         ! "gen_points" with output=1 in python code.

         include "Setting.inc"
         integer :: i, j, k, n2
         integer row(0:nTheta)
         real*8 piset(nTheta), piset2(0:nTheta)

         row = 0
         do i=1,nTheta
            piset(i)=i*1.d0/(nTheta-1)*pi
            piset2(i)=i*2.d0/(nTheta-1)*pi
         enddo
         piset2(0)=0.d0

         row(0) = 1
         row(1) = 1
         k=2
         do j=1,nTheta-2
            n2=nint((nTheta*2-2)*dsin(piset(j)))
            do i=0,n2-1
               row(k) = row(k) + 1
            enddo
            k=k+1
         enddo
         row(k)=1

         return
      End


      Subroutine Solid_Angle_For_Velocity_Volume_Element(solid_angle)
         ! 'solid_angle' = sin(theta).d(theta).d(phi)
         ! 'solanglist' in python code
         use Module_for_NVelocityDirection
         include "Setting.inc"
         external gen_points_for_each_row

         real*8 solid_angle(N_vel_directions)
         integer row(0:nTheta)
         real*8 dphi, lat0, latup, latdown
         integer thetasec, direc
         integer i, k

         call gen_points_for_each_row(row)

         thetasec = nTheta - 1

         do k=1,N_vel_directions
            direc = k-1;  i = 1
            do while (direc+1 .gt. row(i)) 
               direc = direc - row(i)
               i = i + 1
            enddo
            dphi = 2*pi/row(i)
            lat0 = pi/2 - pi*(i-1)/thetasec

            latup = lat0 + pi/thetasec/2
            latdown = lat0 - pi/thetasec/2
            if ( (i .eq. 1) ) then
               latup = pi/2
               latdown = pi/2 - pi/thetasec/2
            else if (i .eq. thetasec+1) then
               latup = -pi/2 + pi/thetasec/2
               latdown = -pi/2
            endif

            solid_angle(k) = (sin(latup)-sin(latdown))*dphi

         enddo

         return
      End


      Subroutine Radial_Component_For_Velocity_Volume_Element(v2dv)
         ! 'v2dv' = v^2 dv (v=vr for initial condition)
         ! 'vollist' in python code.
         include "Setting.inc"
         external Init_Parameter

         real*8, dimension(nEnergy) :: energy_to_speed, v2dv
         real*8 v_spacing(nEnergy+1), half_dv
         real*8 radial_distance_range(nRadial), energy_range(nEnergy), longitude_range(nLong), latitude_range(nLat), latitudeNS_range(nLat_NS)
         real*8 radial_boundary(2), tmax
         integer iE

         call Init_Parameter(radial_distance_range, energy_range, longitude_range, latitude_range, latitudeNS_range, radial_boundary, tmax)
         energy_to_speed = sqrt(energy_range*e*2.d0/mH)

         v_spacing(1)=energy_to_speed(1)/2
         do iE=2, nEnergy
            if (iE .lt. nEnergy) then
               half_dv = 0.5d0*(energy_to_speed(iE+1)-energy_to_speed(iE))
            endif
            !use old half_dv of iE=nEnergy for energy_to_speed(iE+1)
            v_spacing(iE) = energy_to_speed(iE)-half_dv
         enddo
         v_spacing(nEnergy+1) = energy_to_speed(nEnergy)+half_dv

         do iE=1, nEnergy
            v2dv(iE) = (v_spacing(iE+1)**3 - v_spacing(iE)**3)/3.d0
         enddo

         return
      End


      Subroutine calculate_Velocity_Volume_Element(dV2)

         use Module_for_NVelocityDirection
         include "Setting.inc"
         external Solid_Angle_For_Velocity_Volume_Element, Radial_Component_For_Velocity_Volume_Element

         real*8 solid_angle(N_vel_directions)
         real*8, dimension(nEnergy) :: v2dv
         real*8 dV2(nEnergy,N_vel_directions)
         integer iE, iv

         call Solid_Angle_For_Velocity_Volume_Element(solid_angle)
         call Radial_Component_For_Velocity_Volume_Element(v2dv)

         do iv=1,N_vel_directions
            do iE=1,nEnergy
               dV2(iE,iv) = v2dv(iE)*solid_angle(iv)
            enddo
         enddo

         return
      End


      Subroutine calculate_Configuration_Volume_Element(radial_distance_range, lat, dV1)
         ! For RadPres...
         include "Setting.inc"
         real*8 radial_distance_range(nRadial)
         real*8 r, lat, dr1, dV1(nRadial)
         real*8 dlat, dphi
         integer iR, ilat

         dr1 = 0.5d0
         dlat = 15.d0 *pi/180
         dphi = 15.d0 *pi/180
         ! Solid angle for configuration volume element
         do iR=1,nRadial
            r = radial_distance_range(iR)
!            lat = latitude_range(ilat)
!            dV1(iR,ilat) = ((r+dr1)**3 - r**3)/3.d0 * (cos(lat-dlat/2)-cos(lat+dlat/2))*dphi
            dV1(iR) = ((r+dr1)**3 - r**3)/3.d0 * (cos(lat-dlat/2)-cos(lat+dlat/2))*dphi
         enddo

         return
      End


      Subroutine Volume_Element(radial_distance_range,lat, dV)
         ! dV = dx^3 * dv^3
         use Module_for_NVelocityDirection
         include "Setting.inc"
         external calculate_Configuration_Volume_Element, calculate_Velocity_Volume_Element

         real*8 dV1(nRadial), dV2(nEnergy,N_vel_directions)
         real*8 radial_distance_range(nRadial)
         real*8, dimension(N_vel_directions,nRadial,nEnergy) :: dV
         real*8 lat
         integer iE,iR,iv

         call calculate_Configuration_Volume_Element(radial_distance_range, lat, dV1)
         call calculate_Velocity_Volume_Element(dV2)

         do iE=1,nEnergy
            do iR=1,nRadial
               do iv=1,N_vel_directions
                  dV(iv,iR,iE) = dV1(iR)*dV2(iE,iv)
               enddo
            enddo
         enddo

         return
      End


      Subroutine Generate_tag(lon,lat, tag)
         ! Generate "tag" in format "i3.3" considering negative latitudes.
         ! Example: tag = "_lon270_lat000"
         include "Setting.inc"
         real*8 lon, lat
         character*30 tag

         if (lat .ge. 0) then
            write(tag,'(A, i3.3, A, i3.3, A)') "_lon", nint(lon*180/pi), "_lat", nint(lat*180/pi)
         else
            write(tag,'(A, i3.3, A, i3.2, A)') "_lon", nint(lon*180/pi), "_lat", nint(lat*180/pi)
         endif

         return
      End



      Subroutine Get_exobaseBC(nH_BC, TH_BC, rank)

         include "Setting.inc"
         external read_exobaseBC

         real*8, dimension(nbx,nby,nbtperday,start_ydoy-nt_bwd_bc:end_ydoy) :: nH_BC, TH_BC
         real*8, dimension(nbx,nby,nbtperday) :: nH_temp, TH_temp
         integer iday, rank, maxdoy
         character*7 ydoy_str, yearst
         character*100 filename_BC

         if (start_ydoy/1000 .eq. end_ydoy/1000) then
         write(yearst, '(I4.4)') start_ydoy/1000

         do iday=start_ydoy-nt_bwd_bc,end_ydoy
            nH_temp = 0.d0 ; TH_temp=0.d0
            write(ydoy_str,'(I7.7)') iday
            write(yearst, '(I4.4)') start_ydoy/1000
            filename_BC = trim(BC_dir) // trim(yearst) // "/" // trim(ExobaseBC_Model_Name) // "_" // trim(ydoy_str) //  ".bc"
            call read_exobaseBC(filename_BC, nH_temp,TH_temp, rank+12)
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
            nH_temp = 0.d0 ; TH_temp=0.d0
            write(ydoy_str,'(I7.7)') iday
            write(yearst, '(I4.4)') start_ydoy/1000
            filename_BC = trim(BC_dir) // trim(yearst) // "/" // trim(ExobaseBC_Model_Name) // "_" // trim(ydoy_str) //  ".bc"
            call read_exobaseBC(filename_BC, nH_temp,TH_temp, rank+12)
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
            nH_temp = 0.d0 ; TH_temp=0.d0
            write(ydoy_str,'(I7.7)') iday
            write(yearst, '(I4.4)') end_ydoy/1000
            filename_BC = trim(BC_dir) // trim(yearst) // "/" // trim(ExobaseBC_Model_Name) // "_" // trim(ydoy_str) //  ".bc"
            call read_exobaseBC(filename_BC, nH_temp,TH_temp, rank+12)
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




