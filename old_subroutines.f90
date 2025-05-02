!      Subroutine Calculate_Density_Hodge(fin, flags, nH_BC,TH_BC, number_density_at_single_LON_LAT, rank)
      Subroutine Calculate_Escaping_Flux_constBC_Related_to_Forward_Tracing(b_ptl, f_ptl, f_flags, PSD2, energy_range, rank)

         USE SETTING
         use, intrinsic :: ieee_arithmetic

         external calculate_Velocity_Volume_Element
         external GSE2SPH

         real*8, dimension(nvel,nRadial,nEnergy,7) :: b_ptl, f_ptl
         integer, dimension(nvel,nRadial,nEnergy) :: f_flags
         real*8, dimension(nbx,nby) :: number_density_2D
         real*8, dimension(nEnergy,nvel) :: dV2
         real*8, dimension(:,:,:), allocatable :: each_n
         real*8 esc_flux, cexo2
         real*8 pos(3), vel(3), vel2
         real*8 temp_BC, n_BC, vel_BC(3), fac, PSD1
         integer iR,iE,iv, i
         real*8, dimension(nbx,nby,nbtperday,start_ydoy-nt_bwd_bc:end_ydoy) :: nH_BC, TH_BC
         real*8 finlon, finlat, vr
         real*8 current_time, t0, t1 
         integer iflon, iflat, it, rank, quotient
         character*30 fn2D, fn3D
         integer idoy, iday, iE1
         real*8, dimension(nEnergy) :: PSD2, energy_range
         real*8 Enr, w1, w2, deltaE

         vel_BC = 0.d0;
         call calculate_Velocity_Volume_Element(dV2)
!         call calculate_Velocity_Volume_Element_For_Flux(dV2)

         allocate(each_n(nvel,nRadial,nEnergy))
         each_n = 0.d0
         fac = 2.d0*kb/mH
         do iR=1,nRadial      ! Outermost iR-loop
            do iE=1,nEnergy
               do iv=1,nvel
!                  t0 = current_time + fin(iv,iR,iE,1)/86400.    ! unit day
!                  idoy = int(t0)                               ! yyyy+doy
!                  t1 = (t0 - idoy)*86400.                       ! hms in seconds
!                  it = floor(t1/tb_res)+1
!                  if (idoy .lt. start_ydoy-nt_bwd_bc) then ; idoy=start_ydoy-nt_bwd_bc ; it=1 ; endif
                  if (f_flags(iv,iR,iE) .eq. 12) then
                     do i=1,3
                        pos(i) = b_ptl(iv,iR,iE,i+1)
                        vel(i) = b_ptl(iv,iR,iE,i+4)
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

!                     n_BC    = nH_BC(iflon,iflat,it,idoy)
!                     temp_BC = TH_BC(iflon,iflat,it,idoy)
                     n_BC    = 1.2d5
                     temp_BC = 1.d3

!                     vel = (vel - vel_BC)
                     cexo2 = fac*temp_BC
                     vel2 = sum(vel*vel)
                     PSD1 = n_BC * exp(-vel2/cexo2) / (pi*cexo2)**1.5
!                     each_n = PSD1 * dV2(iE,iv) 

                     do i=1,3
                        pos(i) = f_ptl(iv,iR,iE,i+1)
                        vel(i) = f_ptl(iv,iR,iE,i+4)
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

                     ! vr = (v.r)/|r|
                     vr = (pos(1)*vel(1)+pos(2)*vel(2)+pos(3)*vel(3))/sqrt(pos(1)*pos(1)+pos(2)*pos(2)+pos(3)*pos(3))
!                  jE = integer value of enr
!                  PSD2(jE) = PSD2(jE) + PSD1



                  Enr = 0.5*mH*vr*vr/e ! in eV
                  do iE1 = 1, nEnergy - 1
                     if (Enr >= energy_range(iE1) .and. Enr <= energy_range(iE1 + 1)) then
                        deltaE = energy_range(iE1 + 1) - energy_range(iE1)
                        w1 = (energy_range(iE1 + 1) - Enr) / deltaE
                        w2 = (Enr - energy_range(iE1)) / deltaE
                        ! Assign the interpolated density to the grid points
                        PSD2(iE1) = PSD2(iE1) + PSD1 * w1
                        PSD2(iE1 + 1) = PSD2(iE1 + 1) + PSD1 * w2
                        exit
                     end if
                  enddo


!                  PSD2(iflon,iflat,jE) = PSD2(iflon,iflat,jE) + PSD1
!                     number_density_2D(iflon,iflat) = number_density_2D(iflon,iflat) + PSD1*vr

                  endif
               enddo
            enddo
            density_R = sum(each_n(:,iR,:))
         enddo

         return
      End



      Subroutine MSIS_averaged_over_exobase(MSIS_nH,MSIS_TH)

         USE SETTING
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



      Subroutine outptl(ptl,filename)

         use SET_VELOCITY_DIRECTION
         USE SETTING
         IMPLICIT NONE

         real*8, dimension(nvel,nRadial,nEnergy,7) :: ptl
         real one(7)
         integer iR,iE,iv
         character*12 filename

         open(file=filename,unit=21)
         do iE=1,nEnergy
            do iR=1,nRadial
               do iv=1,nvel
                  one = real(ptl(iv,iR,iE,:))
                  write(21,*) one
               enddo
            enddo
         enddo
         close(21)

         return
      End


      Subroutine outind(flags,filename)

         use SET_VELOCITY_DIRECTION
         USE SETTING
         IMPLICIT NONE

         integer, dimension(nvel,nRadial,nEnergy) :: flags
         integer iR,iE,iv
         character*12 filename

         open(file=filename,unit=22)
         do iE=1,nEnergy
            do iR=1,nRadial
               do iv=1,nvel
                  write(22,*) flags(iv,iR,iE)
               enddo
            enddo
         enddo
         close(22)

         return
      End


      Subroutine out_init_binary(init,input_dir,tag)

         use SET_VELOCITY_DIRECTION
         USE SETTING
         IMPLICIT NONE

         real*8, dimension(nvel,nRadial,nEnergy,7) :: init
         real, dimension(:,:,:,:), allocatable :: real_init
         integer nlen
         character*70 input_dir
         character*30 tag
         character*100 filename

         allocate(real_init(nvel,nRadial,nEnergy,7))
         real_init = real(init)

         inquire(iolength=nlen) real_init 
         filename = trim(input_dir) // "EXO_ini" // trim(tag) // ".data"
      print*, "Generate init file: ", filename
         open(file=filename,unit=21,form='unformatted',access='direct',recl=nlen,status='replace')
         write(21,rec=1) real_init
         close(21)

         init = real_init * 1.d0
         deallocate(real_init)

         return
      End


      Subroutine out_fin_binary(fin,input_dir,tag)

         use SET_VELOCITY_DIRECTION
         USE SETTING
         IMPLICIT NONE

         real*8, dimension(nvel,nRadial,nEnergy,7) :: fin
         real, dimension(:,:,:,:), allocatable :: real_fin
         integer nlen
         character*70 input_dir
         character*30 tag
         character*100 filename

         allocate(real_fin(nvel,nRadial,nEnergy,7))
         real_fin = real(fin)

         inquire(iolength=nlen) real_fin 
         filename = trim(input_dir) // "EXO_fin" // trim(tag) // ".data"
      print*, "Generate fin file: ", filename
         open(file=filename,unit=22,form='unformatted',access='direct',recl=nlen,status='replace')
         write(22,rec=1) real_fin
         close(22)

         return
      End


      Subroutine out_ind_binary(flags,input_dir,tag)

         use SET_VELOCITY_DIRECTION
         USE SETTING
         IMPLICIT NONE

         integer, dimension(nvel,nRadial,nEnergy) :: flags
         integer nlen
         character*70 input_dir
         character*30 tag
         character*100 filename

         inquire(iolength=nlen) flags
         filename = trim(input_dir) // "EXO_ind" // trim(tag) // ".data"
      print*, "Generate ind file: ", filename
         open(file=filename,unit=23,form='unformatted',access='direct',recl=nlen,status='replace')
         write(23,rec=1) flags
         close(23)

         return
      End


     Subroutine outRuntime(Runtime_dist,filename)

         USE SETTING
         IMPLICIT NONE
         real, dimension(nRadial,nEnergy) :: Runtime_dist
         integer iR,iE
         character*16 filename

         open(file=filename,unit=39)
         do iE=1,nEnergy
            do iR=1,nRadial
               write(39,*) Runtime_dist(iR,iE)
            enddo
         enddo
         close(39)

         return
      End
      

      
      Subroutine read_ind_binary(flags,input_dir,tag,thread_num)

         use SET_VELOCITY_DIRECTION
         USE SETTING
         IMPLICIT NONE

         integer, dimension(nvel,nRadial,nEnergy) :: flags
         integer nlen, iexist, thread_num, IO_unit
         character*70 input_dir
         character*30 tag
         character*100 filename

!         thread_num = omp_get_thread_num()
         IO_unit=thread_num+1
         print*, 'IO_unit at ind = ', IO_unit

         filename = trim(input_dir) // 'EXO_ind' // trim(tag) // '.data'
         print*, "Read ind file: ", filename
         inquire(file=filename, exist=iexist)
         if (iexist .eq. 0) then
            print*, "*** ERROR!! FILE IS NOT EXIST!! ***"
            flags=-1
            stop
         else
            inquire(iolength=nlen) flags
            open(file=filename,unit=IO_unit,form='unformatted', &
               access='direct',action='read',recl=nlen,status='old')
            read(IO_unit,rec=1) flags
            close(IO_unit)
         endif

         return
      End


      Subroutine read_fin_binary(fin,input_dir,tag,thread_num)

         use SET_VELOCITY_DIRECTION
         USE SETTING
         IMPLICIT NONE
         real*8, dimension(nvel,nRadial,nEnergy,7) :: fin
         real, dimension(:,:,:,:), allocatable :: real_fin
         integer nlen, thread_num, IO_unit, iexist
         character*70 input_dir
         character*30 tag
         character*100 filename

         allocate(real_fin(nvel,nRadial,nEnergy,7))

!         thread_num = omp_get_thread_num()
         IO_unit=thread_num+512
         print*, 'IO_unit at fin = ', IO_unit

         filename = trim(input_dir) // 'EXO_fin' // trim(tag) // '.data'
         print*, "Read fin file: ", filename
         inquire(file=filename, exist=iexist)
         if (iexist .eq. 0) then
            print*, "*** ERROR!! FILE IS NOT EXIST!! ***"
            real_fin = 0.d0
            stop
         else
            inquire(iolength=nlen) real_fin
            open(file=filename,unit=IO_unit,form='unformatted', &
               access='direct',action='read',recl=nlen,status='old')
            read(IO_unit,rec=1) real_fin
            close(IO_unit)
         endif

         fin = real_fin * 1.d0
         deallocate(real_fin)

         return
      End


      Subroutine write_density_1D(density_1D,tag)

         USE SETTING
         IMPLICIT NONE

         real*8 density_1D(nRadial)
         real, dimension(:), allocatable :: real_density_1D
         integer nlen
         character*30 tag
         character*40 filename

         allocate(real_density_1D(nRadial))
         real_density_1D = real(density_1D)

         filename = 'EXO_Density_1D' // trim(tag) // '.data' 
         inquire(iolength=nlen) real_density_1D
         open(file=filename,unit=41,form='unformatted',access='direct',recl=nlen,status='replace')
         write(41,rec=1) real_density_1D
         close(41)

         deallocate(real_density_1D)

         return
      End


      Subroutine write_density_3D(density_3D,tag)

         USE SETTING
         IMPLICIT NONE
         
         real*8 density_3D(nRadial,nLong,nLat_NS)
         real, dimension(:,:,:), allocatable :: real_density_3D
         integer nlen
         character*30 tag
         character*100 filename

         allocate(real_density_3D(nRadial,nLong,nLat_NS))
         real_density_3D = real(density_3D)

         filename = trim(outdir) // 'EXO_Density_3D' // trim(tag) //    '.data' 
         inquire(iolength=nlen) real_density_3D
         open(file=filename,unit=43,form='unformatted',access='direct',recl=nlen,status='replace')
         write(43,rec=1) real_density_3D
         close(43)

         deallocate(real_density_3D)

         return
      End


      Subroutine Write_2D_Real(fn2D, arr2D, nx,ny)

         USE SETTING
         IMPLICIT NONE

         integer nx, ny, nlen
         real*8, dimension(nx,ny) :: arr2D
         real, dimension(:,:), allocatable :: real_arr2D
         character*30 fn2D
         character*100 filename

         allocate(real_arr2D(nx,ny))
         real_arr2D = real(arr2D)

         filename = trim(outdir) // trim(fn2D) // '.data'
         inquire(iolength=nlen) real_arr2D
         open(file=filename, unit=45, form='unformatted',access='direct',recl=nlen,status='replace')
         write(45,rec=1) real_arr2D
         close(45)

         deallocate(real_arr2D)

         return
      End


      Subroutine Write_3D_Real(fn3D, arr3D, nx,ny,nz)

         USE SETTING
         IMPLICIT NONE

         integer nx, ny, nz, nlen
         real*8, dimension(nx,ny,nz) :: arr3D
         real, dimension(:,:,:), allocatable :: real_arr3D
         character*30 fn3D
         character*100 filename

         allocate(real_arr3D(nx,ny,nz))
         real_arr3D = real(arr3D)

         filename = trim(outdir) // trim(fn3D) // '.data'
            print*, filename
         inquire(iolength=nlen) real_arr3D
         open(file=filename, unit=46, form='unformatted',access='direct',recl=nlen,status='replace')
         write(46,rec=1) real_arr3D
         close(46)

         deallocate(real_arr3D)

         return
      End




      Subroutine Write_ESC_FLUX_2D(density_2D)

         USE PHYSICS_TAG
         USE SETTING
         IMPLICIT NONE
         real*8 density_2D(nbx,nby)
         real, dimension(:,:), allocatable :: real_density_2D
         integer nlen
         character*100 filename

         call Physics_tag()
         allocate(real_density_2D(nbx,nby))
         real_density_2D = real(density_2D)

         filename = trim(outdir) // 'ESC_FLUX_2D_' // trim(tag_phys) // '_' // trim(tag0) // '.data' 
         inquire(iolength=nlen) real_density_2D
         open(file=filename,unit=42,form='unformatted',access='direct',recl=nlen,status='replace')
         write(42,rec=1) real_density_2D
         close(42)

         deallocate(real_density_2D)

         return
      End


      Subroutine Write_ESC_FLUX_1D(density_1D)

         USE PHYSICS_TAG
         USE SETTING
         IMPLICIT NONE

         real*8 density_1D(nEnergy)
         real, dimension(:), allocatable :: real_density_1D
         integer nlen
         character*100 filename

         call Physics_tag()
         allocate(real_density_1D(nEnergy))
         real_density_1D = real(density_1D)

         filename = trim(outdir) // 'ESC_FLUX_1D_' // trim(tag_phys) // '_' // trim(tag0) // '.data' 
         inquire(iolength=nlen) real_density_1D
         open(file=filename,unit=42,form='unformatted',access='direct',recl=nlen,status='replace')
         write(42,rec=1) real_density_1D
         close(42)

         deallocate(real_density_1D)

         return
      End



      Subroutine Generate_tag(lon,lat, tag)
         ! Generate "tag" in format "i3.3" considering negative latitudes.
         ! Example: tag = "_lon270_lat000"
         USE SETTING
         real*8 lon, lat
         character*30 tag

         if (lat .ge. 0) then
            write(tag,'(A, i3.3, A, i3.3, A)') "_lon", nint(lon*180/pi), "_lat", nint(lat*180/pi)
         else
            write(tag,'(A, i3.3, A, i3.2, A)') "_lon", nint(lon*180/pi), "_lat", nint(lat*180/pi)
         endif

         return
      End



   Subroutine Initialize(input_dir, init, fin, flags, tag, thread_num)

      USE SETTING
      external Init_Parameter, read_ind_binary, read_fin_binary

      real*8, dimension(nvel,nRadial,nEnergy,7) :: init, fin
      integer, dimension(nvel,nRadial,nEnergy) :: flags

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






   Subroutine Forward_Tracing_particle(ptl,flags, init_flags, tmax, Lya, current_time)

!         use omp_lib
      USE SETTING
      USE GRID_PARAMETERS, only: radial_boundary
      external rk4, calculate_final_timestep

      real*8, dimension(nvel,nRadial,nEnergy,7) :: ptl
      integer, dimension(nvel,nRadial,nEnergy) :: flags, init_flags
      real*8, dimension(7) :: one, old
      integer :: flag     ! 0: orbiting Earth t<tmax;   1: into exobase;  2: out of outer boundary;  3: orbiting but t>tmax
      real*8 radial_distance, radial_distance_old
      integer :: iR, iE, iv, i
      real*8 tmax,dt, vt,vt_old,dv, ds
      real*8, parameter :: max_ds = 1.d6
      real*8 x0, f0, current_time, trace_time
      real*8, dimension(start_ydoy_index:end_ydoy_index) :: Lya
      integer ydoy, ii

      do iE=1,nEnergy
         do iR=1,nRadial
         do iv=1,nvel
            if(init_flags(iv,iR,iE) .eq. 1) then

            flag=10; 
            do i=1,7;   one(i) = ptl(iv,iR,iE,i);   enddo

            radial_distance_old = sqrt(one(2)**2+one(3)**2+one(4)**2)
            vt = sqrt(one(5)**2 + one(6)**2 + one(7)**2)
            trace_time = current_time - 1e-5

            do while (abs(one(1)) < tmax)

               radial_distance_old = sqrt(one(2)**2+one(3)**2+one(4)**2)
               vt_old = sqrt(one(5)**2 + one(6)**2 + one(7)**2)

               dt = +1.d0*max_ds / vt_old     ! -1e6 or 4e6 is a "factor" in python code. The maximum distance jump at single time step.
               if (mod(trace_time-1,1000.0) .gt. 500) then               ! eg. trace_time=2010000.98, then it should be 2009365.98, 
                  ii=1000-mod(int(trace_time-1),1000)                  ! eg. trace_time-1 = 2009999.98, ii=1000-999=1
                  if (mod(int((trace_time-1)/1000),4) .eq. 0) then
                     trace_time = int((trace_time-1)/1000)*1000 + (367-ii) + mod(trace_time,1.0)      ! For leap years (400-year period is not applied)
                  else
                     trace_time = int((trace_time-1)/1000)*1000 + (366-ii) + mod(trace_time,1.0)      ! eg. trace_time = 2009000+365+0.98 = 2009365.98
                  endif
               endif

               ydoy = int(trace_time)
               f0 = Lya(ydoy)
               if (int(trace_time + dt/86400) .ne. int(trace_time)) then
!                     dt = (int(trace_time)-trace_time)*86400 - 1e-5    ! trace_time always hits the time (00:00:00) for daily-varying Lya.
                  dt = (int(trace_time)+1-trace_time)*86400 + 1e-5  ! Forward Tracing
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
                  flag = 11
                  call calculate_final_timestep(old,one,dt,f0)
               else if (radial_distance .gt. radial_boundary(2)) then
                  flag = 12
               endif

               if (flag > 10) then
                  exit
               endif

            enddo ! end while

            do i=1,7;   ptl(iv,iR,iE,i) = one(i);   enddo
            flags(iv,iR,iE) = flag

         endif ! init_flag=1

      enddo ! iv

!print*, int(iE,1), int(iR,1)

      enddo ! iR

!print*, int(iE,1)

      enddo ! iE

      return
   End


