      Subroutine Init_Parameter(radial_distance_range, energy_range, &
                 longitude_range, latitude_range, latitudeNS_range,  radial_boundary, tmax)

            include "Setting.inc"
            integer iR, iE, ilon, ilat
            real*8 radial_distance_range(nRadial), energy_range(nEnergy)
            real*8 longitude_range(nLong), latitude_range(nLat), latitudeNS_range(nLat_NS)
            real*8 radial_boundary(2), tmax
            integer index_Emin

            do iR=1,nRadial
               radial_distance_range(iR) = RadialRange_min + (iR-1)*dR
!               radial_distance_range(iR) = 2.d0**(1.d0+(iR-1)/3.d0)   !  log-scale
            enddo
            radial_distance_range = radial_distance_range * Re

            if (nEnergy .eq. 121) then  ! for 0.0025 - 10 eV
!            index_Emin = -33              ! Emin at 0.001 eV
                  index_Emin = -20              ! Emin at 0.0025 eV
                  do iE=index_Emin,index_Emin+nEnergy-1
                        energy_range(iE+1-index_Emin) = 0.01 * 1000.d0 ** (iE/100.d0)
                  enddo
            else
                  if (nEnergy .eq. 61) then  ! for 0.001 - 1 eV
                        do iE=1,nEnergy 
                              energy_range(iE) = 10.d0 ** ((iE-1)/20.d0 - 3.d0)
                        enddo
                  else
                        print*, 'ERROR: Set a proper nEnergy'
                  endif
            endif
 
            do ilon=1,nLong
                  longitude_range(ilon) = (ilon-1)*pi/180.d0*geores
            enddo
            do ilat=1,nLat
                  latitude_range(ilat) = (ilat-1)*pi/180.d0*geores
            enddo
            do ilat=1,nLat_NS
                  latitudeNS_range(ilat) = ((ilat-1.d0)*geores-90.d0) *pi/180.d0
            enddo

            radial_boundary(1) = inner_boundary
            radial_boundary(2) = outer_boundary

            tmax = ntmax * 86400.d0          ! 60 days

            return
      End


      Subroutine Init_Particles(ptl, radial_distance_range, energy_range, lon,lat)

            use Module_for_NVelocityDirection
            include "Setting.inc"
            external gen_points

            integer iR, iE
            real*8, dimension(N_vel_directions,nRadial,nEnergy,7) :: ptl
            real*8, dimension(N_vel_directions,3) :: vel_dir
            real*8 radial_distance_range(nRadial), energy_range(nEnergy)
            real*8 energy_to_speed, lon,lat, cos_lat,sin_lat, cos_lon,sin_lon

            ptl=0.d0

            call gen_points(vel_dir)
            sin_lat = sin(lat) ;    cos_lat = cos(lat)
            sin_lon = sin(lon) ;    cos_lon = cos(lon)

            do iE=1, nEnergy
                  energy_to_speed = sqrt(energy_range(iE)*e*2.d0/mH)
                  do iR=1, nRadial
                        ptl(:,iR,iE,2) = radial_distance_range(iR)*cos_lat*cos_lon        ! X
                        ptl(:,iR,iE,3) = radial_distance_range(iR)*cos_lat*sin_lon        ! Y
                        ptl(:,iR,iE,4) = radial_distance_range(iR)*sin_lat                ! Z
                        ptl(:,iR,iE,5) = vel_dir(:,1) * energy_to_speed       ! Vx
                        ptl(:,iR,iE,6) = vel_dir(:,2) * energy_to_speed       ! Vy
                        ptl(:,iR,iE,7) = vel_dir(:,3) * energy_to_speed       ! Vz
                  enddo
            enddo

            return
      End


