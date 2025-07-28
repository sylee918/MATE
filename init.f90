   Subroutine Init_Particles(ptl)

      USE SETTING
      USE GRID_PARAMETERS,        only: radial_distance_range, energy_range
      USE SET_VELOCITY_DIRECTION, only: gen_points
      USE MPI_MATE,               only: lon, lat, rad

      integer iE, rad2
      real*8, dimension(nvel,nEnergy,7) :: ptl
      real*8, dimension(nvel,3) :: vel_dir
      real*8 energy_to_speed, cos_lat,sin_lat, cos_lon,sin_lon

      ptl=0.d0

      call gen_points(vel_dir)
      sin_lat = sin(lat)  ;  cos_lat = cos(lat)
      sin_lon = sin(lon)  ;  cos_lon = cos(lon)

      do iE=1, nEnergy
         energy_to_speed = sqrt(energy_range(iE)*e*2.d0/mH)
!            ptl(:,iR,iE,2) = radial_distance_range(iR)*cos_lat*cos_lon        ! X
!            ptl(:,iR,iE,3) = radial_distance_range(iR)*cos_lat*sin_lon        ! Y
!            ptl(:,iR,iE,4) = radial_distance_range(iR)*sin_lat                ! Z
            ptl(:,iE,2) = rad*cos_lat*cos_lon        ! X
            ptl(:,iE,3) = rad*cos_lat*sin_lon        ! Y
            ptl(:,iE,4) = rad*sin_lat                ! Z
            ptl(:,iE,5) = vel_dir(:,1) * energy_to_speed       ! Vx
            ptl(:,iE,6) = vel_dir(:,2) * energy_to_speed       ! Vy
            ptl(:,iE,7) = vel_dir(:,3) * energy_to_speed       ! Vz
      enddo

      return
   End


