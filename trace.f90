      Subroutine rk1(temp, dt, k_out, f0)

         include "Setting.inc"

         real*8 temp(6), dt
         real*8, dimension(3) :: pos, vel
         real*8, dimension(6) :: deriv_G, deriv_R, deriv_C, k_out
         real*8 rho, f0

         pos = temp(1:3)
         vel = temp(4:6)

         deriv = 0.d0
         deriv(1:3) = vel 

         deriv_G = 0.d0
         deriv_R = 0.d0
         deriv_C = 0.d0

         deriv_G(4:6) = -GM * pos / (pos(1)**2 + pos(2)**2 + pos(3)**2)**1.5   

         ! Add Radialtion Pressure
            ! No radiation by Earth's shadow
         rho = sqrt(pos(2)**2 + pos(3)**2)
         if (pos(1) .gt. 0 .or. rho .gt. Re) then
            deriv_R(4) = - arad*f0
         endif

         ! Add Coriolis Force
         deriv_C(4) = 2*Wrot*vel(2) 
         deriv_C(5) = -2*Wrot*vel(1)

         deriv = deriv_G*i_EarthGravity + deriv_R*i_SolarRadiationPressure + deriv_C*i_CoriolisForce_GSE
         k_out = dt * deriv

         return
      End


      Subroutine rk4(one,dt,f0)

         include "Setting.inc"
         external rk1

         real*8 one(7), temp(6), dt, f0
         real*8, dimension(6) :: k1,k2,k3,k4
         integer i

         temp = one(2:7)

         call rk1(temp, dt, k1, f0)

         temp = one(2:7) + 0.5d0*k1; 
         call rk1(temp, dt, k2, f0)

         temp = one(2:7) + 0.5d0*k2
         call rk1(temp, dt, k3, f0)

         temp = one(2:7) + k3
         call rk1(temp, dt, k4, f0)

         do i=1,6
            one(i+1) = one(i+1) + (k1(i)+ 2*k2(i) + 2*k3(i) + k4(i))/6.d0
         enddo
         one(1) = one(1) + dt
      
      End
        

      Subroutine calculate_final_timestep(old,new,radial_boundary,dt,f0)      ! Do interpolation for dt_final

         include "Setting.inc"
         external rk4

         real*8, dimension(7) :: old, new
         real*8 radial_boundary(2), dt, k1(6)
         real*8 r_old, r_new, fac, dt_final, f0

         r_old = sqrt(old(2)**2+old(3)**2+old(4)**2)
         r_new = sqrt(new(2)**2+new(3)**2+new(4)**2)
         if (r_old .lt. r_new) then
            print*, "R_old > R_new"
            stop
         endif

         fac = (r_old - radial_boundary(1)) / (r_old - r_new)
         dt_final = dt * fac

         call rk4(old,dt_final,f0)
         new = old

!         call rk1(old(2:7),dt_final,k1)     !   Just apply Euler method for the last step.
!         new(2:7) = old(2:7) + k1
!         new(1) = old(1) + dt_final

         return
      End


      Subroutine Trace_particle(ptl,flags, radial_boundary, tmax, Lya, current_time)

!         use omp_lib
         include "Setting.inc"
         external rk4, calculate_final_timestep

         real*8, dimension(N_vel_directions,nRadial,nEnergy,7) :: ptl
         integer :: flags(N_vel_directions,nRadial,nEnergy)
         real*8, dimension(7) :: one, old
         integer :: flag     ! 0: orbiting Earth t<tmax;   1: into exobase;  2: out of outer boundary;  3: orbiting but t>tmax
         real*8 radial_boundary(2), radial_distance, radial_distance_old
         integer :: iR, iE, iv, i
         real*8 tmax,dt, vt,vt_old,dv, ds
         real*8, parameter :: max_ds = 1.d6
         real*8 Emec0, Emec1, Eerr
         real*8 x0, f0, current_time, trace_time
         real*8, dimension(start_ydoy_index:end_ydoy_index) :: Lya
         integer ydoy, ii

         do iE=1,nEnergy
          do iR=1,nRadial
           do iv=1,N_vel_directions
!         do iE=41,41
!          do iR=17,17
!           do iv=3418,3418
!print*, iv
            flag=0; 
            do i=1,7;   one(i) = ptl(iv,iR,iE,i);   enddo
!open(file='traj_41_17_3418_GR.out',unit=11)
!write(11,*) real(one)

            radial_distance_old = sqrt(one(2)**2+one(3)**2+one(4)**2)
            vt = sqrt(one(5)**2 + one(6)**2 + one(7)**2)
            trace_time = current_time - 1e-5
!print*, 't0', trace_time, current_time

            do while (abs(one(1)) < tmax)

               radial_distance_old = sqrt(one(2)**2+one(3)**2+one(4)**2)
               vt_old = sqrt(one(5)**2 + one(6)**2 + one(7)**2)

               dt = -1.d0*max_ds / vt_old     ! -1e6 or 4e6 is a "factor" in python code. The maximum distance jump at single time step.
               if (mod(trace_time-1,1000.0) .gt. 500) then               ! eg. trace_time=2010000.98, then it should be 2009365.98, 
                  ii=1000-mod(int(trace_time-1),1000)                  ! eg. trace_time-1 = 2009999.98, ii=1000-999=1
                  if (mod(int((trace_time-1)/1000),4) .eq. 0) then
                     trace_time = int((trace_time-1)/1000)*1000 + (367-ii) + mod(trace_time,1.0)      ! For leap years (400-year period is not applied)
                  else
!                     print*, "trace_time1", trace_time
                     trace_time = int((trace_time-1)/1000)*1000 + (366-ii) + mod(trace_time,1.0)      ! eg. trace_time = 2009000+365+0.98 = 2009365.98
!                     print*, "trace_time2", trace_time
                  endif
               endif

               ydoy = int(trace_time)
!print*, 'dt = ', dt, ydoy, current_time, trace_time
!print*, trace_time + dt/86400, int(trace_time + dt/86400), int(trace_time)
               f0 = Lya(ydoy)
               if (int(trace_time + dt/86400) .ne. int(trace_time)) then
                  dt = (int(trace_time)-trace_time)*86400 - 1e-5    ! trace_time always hits the time (00:00:00) for daily-varying Lya.
!print*, "modified dt = ", dt
                  if (abs(dt) .lt. 1e-6) then
                     print*, "ERROR: dt is too small"
                     stop
                  endif
               endif

               old = one
   100 continue
               call rk4(one,dt,f0)
               trace_time = current_time + one(1)/86400  ! one(2) < 0
               ! FIX ME (if time cross year)
!print*, "after push", trace_time

               radial_distance = sqrt(one(2)**2+one(3)**2+one(4)**2)
               vt = sqrt(one(5)**2 + one(6)**2 + one(7)**2)
               ds = sqrt((old(2)-one(2))**2+(old(3)-one(3))**2+(old(4)-one(4))**2)
               dv = vt-vt_old

               ! If the solution is diverging ...
               if (abs(ds/radial_distance_old) .gt. 1e-1 .or. abs(ds)/max_ds .gt. 1.2 .or. dv/vt_old .gt. 10) then
                   dt=dt/2
                   one=old
                   goto 100
               endif

               radial_distance = sqrt(one(2)**2+one(3)**2+one(4)**2)
               if (radial_distance .lt. radial_boundary(1)) then
                  flag = 1
                  call calculate_final_timestep(old,one,radial_boundary,dt,f0)
!if (abs(one(1)) .gt. 15.d0*86400) then
!   print*, abs(one(1))/86400, iE,iR,iv
!   stop
!endif
               else if (radial_distance .gt. radial_boundary(2)) then
                  flag = 2
               endif

!write(11,*) real(one)

               if (flag > 0) then
                  exit
               endif

            enddo ! end while

            do i=1,7;   ptl(iv,iR,iE,i) = one(i);   enddo
            flags(iv,iR,iE) = flag

         enddo ! iv

!print*, int(iE,1), int(iR,1)

         enddo ! iR

!print*, int(iE,1)

         enddo ! iE

!   close(11)
!   stop
!   print*, flag, one/Re

         return
      End

