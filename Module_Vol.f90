MODULE VOLUME_ELEMENT

   USE SETTING
   USE SET_VELOCITY_DIRECTION, only: gen_points_for_each_row
   USE GRID_PARAMETERS
   IMPLICIT NONE

   real*8, allocatable, dimension(:,:) :: dV2
   real*8, allocatable, dimension(:) :: solid_angle


   contains

      Subroutine Solid_Angle_For_Velocity_Volume_Element(solid_angle)
      ! 'solid_angle' = sin(theta).d(theta).d(phi)
      ! 'solanglist' in python code
      real*8, dimension(nvel) :: solid_angle
      integer row(0:nTheta)
      real*8 dphi, lat0, latup, latdown
      integer thetasec, direc
      integer i, k

      call gen_points_for_each_row(row)

      thetasec = nTheta - 1

      do k=1,nvel
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

      real*8 v_spacing(nEnergy+1), half_dv, energy_to_speed(nEnergy), v2dv(nEnergy)
      integer iE

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

      integer iE, iv
      real*8, dimension(nEnergy,nvel) :: dV2
      real*8, dimension(nvel) :: solid_angle
      real*8, dimension(nEnergy) :: v2dv

      call Solid_Angle_For_Velocity_Volume_Element(solid_angle)
      call Radial_Component_For_Velocity_Volume_Element(v2dv)

      do iv=1,nvel
         do iE=1,nEnergy
            dV2(iE,iv) = v2dv(iE)*solid_angle(iv)
         enddo
      enddo

      return
   End


   Subroutine Radial_Component_Of_Velocity_Volume_Element_For_Flux(v3dv)
      ! 'v3dv' = v^3 dv (v=vr for initial condition)
      ! 'vollist' in python code.
      IMPLICIT NONE

      real*8, dimension(nEnergy) :: energy_to_speed, v3dv
      real*8 v_spacing(nEnergy+1), half_dv
      integer iE

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
         v3dv(iE) = (v_spacing(iE+1)**4 - v_spacing(iE)**4)/4.d0
      enddo

      return
   End


   Subroutine calculate_Velocity_Volume_Element_For_Flux(dV2)

      real*8 solid_angle(nvel)
      real*8, dimension(nEnergy) :: v3dv
      real*8 dV2(nEnergy,nvel)
      integer iE, iv

      call Solid_Angle_For_Velocity_Volume_Element(solid_angle)
      call Radial_Component_Of_Velocity_Volume_Element_For_Flux(v3dv)

      do iv=1,nvel
         do iE=1,nEnergy
            dV2(iE,iv) = v3dv(iE)*solid_angle(iv)
         enddo
      enddo

      return
   End


   Subroutine calculate_Configuration_Volume_Element(lat, dV1)
      ! For RadPres...

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


   Subroutine Volume_Element2(lat, dV)
      ! dV = dx^3 * dv^3

      real*8 dV1(nRadial)
      real*8, dimension(nvel,nRadial,nEnergy) :: dV
      real*8 lat
      integer iE,iR,iv

      call calculate_Configuration_Volume_Element(lat, dV1)
      call calculate_Velocity_Volume_Element(dV2)

      do iE=1,nEnergy
         do iR=1,nRadial
            do iv=1,nvel
               dV(iv,iR,iE) = dV1(iR)*dV2(iE,iv)
            enddo
         enddo
      enddo

      return
   End

END MODULE VOLUME_ELEMENT
