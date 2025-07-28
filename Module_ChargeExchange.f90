Module ChargeExchange

   USE SETTING
   IMPLICIT NONE

   integer :: nx=201, ny=201, nz=201
   real*8, dimension(nx,ny,nz) :: nps, Tps
   real*8, dimension(nx) :: xps, yps, zps


   contains

   Subroutine Calculate_ChargeExchange(iE,iv,current_time, ICX)

      USE SETTING
      IMPLICIT NONE

      integer :: iE, iv
      real*8 :: current_time, ICX
      real*8, dimension(7) :: ptl0

      call Retrieve_initptl(iE,iv, ptl0)
      ptl0(1) = current_time

   End Subroutine Calculate_ChargeExchange


   Subroutine Retrieve_initptl(iE,iv, ptl0)

      real*8, dimension(nvel,3) :: vel_dir


      call gen_points(vel_dir)
      sin_lat = sin(lat)  ;  cos_lat = cos(lat)
      sin_lon = sin(lon)  ;  cos_lon = cos(lon)

      energy_to_speed = sqrt(energy_range(iE)*e*2.d0/mH)
      ptl0(2) = rad*cos_lat*cos_lon        ! X
      ptl0(3) = rad*cos_lat*sin_lon        ! Y
      ptl0(4) = rad*sin_lat                ! Z
      ptl0(5) = vel_dir(iv,1) * energy_to_speed       ! Vx
      ptl0(6) = vel_dir(iv,2) * energy_to_speed       ! Vy
      ptl0(7) = vel_dir(iv,3) * energy_to_speed       ! Vz


   End Subroutine Retrieve_initptl


   Subroutine Read_Plasmasphere(nps, Tps)

      integer :: i, j, k

      do i=1,nx
         xps(i) = -10.d0 + (i-1)*0.1
      enddo
      yps = xps ; zps = xps

      do k=1,nz
         do j=1,ny
            do i=1,nx
               nps(i,j,k) = 100.d0
            enddo
         enddo
      enddo


   End Subroutine Read_Plasmasphere


End Module ChargeExchange