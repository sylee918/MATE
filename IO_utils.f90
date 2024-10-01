      Subroutine outptl(ptl,filename)

         use Module_for_NVelocityDirection
         include "Setting.inc"

         real*8, dimension(N_vel_directions,nRadial,nEnergy,7) :: ptl
         real one(7)
         integer iR,iE,iv
         character*12 filename

         open(file=filename,unit=21)
         do iE=1,nEnergy
            do iR=1,nRadial
               do iv=1,N_vel_directions
                  one = real(ptl(iv,iR,iE,:))
                  write(21,*) one
               enddo
            enddo
         enddo
         close(21)

         return
      End


      Subroutine outind(flags,filename)

         use Module_for_NVelocityDirection
         include "Setting.inc"

         integer, dimension(N_vel_directions,nRadial,nEnergy) :: flags
         integer iR,iE,iv
         character*12 filename

         open(file=filename,unit=22)
         do iE=1,nEnergy
            do iR=1,nRadial
               do iv=1,N_vel_directions
                  write(22,*) flags(iv,iR,iE)
               enddo
            enddo
         enddo
         close(22)

         return
      End


      Subroutine out_init_binary(init,input_dir,tag)

         use Module_for_NVelocityDirection
         include "Setting.inc"
         real*8, dimension(N_vel_directions,nRadial,nEnergy,7) :: init
         real, dimension(:,:,:,:), allocatable :: real_init
         integer nlen
         character*70 input_dir
         character*30 tag
         character*100 filename

         allocate(real_init(N_vel_directions,nRadial,nEnergy,7))
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

         use Module_for_NVelocityDirection
         include "Setting.inc"
         real*8, dimension(N_vel_directions,nRadial,nEnergy,7) :: fin
         real, dimension(:,:,:,:), allocatable :: real_fin
         integer nlen
         character*70 input_dir
         character*30 tag
         character*100 filename

         allocate(real_fin(N_vel_directions,nRadial,nEnergy,7))
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

         use Module_for_NVelocityDirection
         include "Setting.inc"
         integer, dimension(N_vel_directions,nRadial,nEnergy) :: flags
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

         include "Setting.inc"

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

         use Module_for_NVelocityDirection
         include "Setting.inc"
         integer, dimension(N_vel_directions,nRadial,nEnergy) :: flags
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

         use Module_for_NVelocityDirection
         include "Setting.inc"
         real*8, dimension(N_vel_directions,nRadial,nEnergy,7) :: fin
         real, dimension(:,:,:,:), allocatable :: real_fin
         integer nlen, thread_num, IO_unit, iexist
         character*70 input_dir
         character*30 tag
         character*100 filename

         allocate(real_fin(N_vel_directions,nRadial,nEnergy,7))

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

         include "Setting.inc"
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

         include "Setting.inc"
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

         include "Setting.inc"
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

         include "Setting.inc"
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



      Subroutine write_density_4D(density_4D,iday)

         use Module_Physics_tag
         include "Setting.inc"
         
         real*8 density_4D(nRadial,nLong,nLat_NS,ntperday)
         real, dimension(:,:,:,:), allocatable :: real_density_4D
         integer iday, nlen
         character*10 dayst
         character*100 filename

         allocate(real_density_4D(nRadial,nLong,nLat_NS,ntperday))
         real_density_4D = real(density_4D)

         write(dayst, '(I7.7)') iday
         filename = trim(outdir) // 'MATE_nH_' // trim(tag_phys) // '_' // trim(tag0) // '_' // trim(dayst) // '.data'
         inquire(iolength=nlen) real_density_4D
         open(file=filename,unit=45,form='unformatted',access='direct',recl=nlen,status='replace')
         write(45,rec=1) real_density_4D
         close(45)

         deallocate(real_density_4D)

         return
      End


      Subroutine Write_ESC_FLUX_2D(density_2D)

         use Module_Physics_tag
         include "Setting.inc"
         real*8 density_2D(nbx,nby)
         real, dimension(:,:), allocatable :: real_density_2D
         integer nlen
         character*100 filename

         allocate(real_density_2D(nbx,nby))
         real_density_2D = real(density_2D)

         filename = trim(outdir) // 'ESC_FLUX_2D' // '_' // trim(tag_phys) // '_' // trim(tag0) // '.data' 
         inquire(iolength=nlen) real_density_2D
         open(file=filename,unit=42,form='unformatted',access='direct',recl=nlen,status='replace')
         write(42,rec=1) real_density_2D
         close(42)

         deallocate(real_density_2D)

         return
      End


      Subroutine read_exobaseBC(filename, nH_temp,TH_temp, thread_num)

         include "Setting.inc"
         real*8, dimension(nbx,nby,nbtperday) :: nH_temp, TH_temp
         real, dimension(:,:,:), allocatable :: nH_real, TH_real
         integer nlen, thread_num, IO_unit, iexist
         character*100 filename

         allocate(nH_real(nbx,nby,nbtperday),TH_real(nbx,nby,nbtperday))
         IO_unit=thread_num+600

         print*, "Read exobase BC file: ", filename
         inquire(file=filename, exist=iexist)
         if (iexist .eq. 0) then
            print*, "File is not exist: ", filename
         else
            inquire(iolength=nlen) nH_real
            nlen=nlen*2

            open(file=filename,unit=IO_unit,form='unformatted', &
               access='direct',action='read',recl=nlen,status='old')
            read(IO_unit,rec=1) nH_real, TH_real
            close(IO_unit)
         endif

         nH_temp = nH_real*1.d0
         TH_temp = TH_real*1.d0
         deallocate(nH_real,TH_real)

         return
      End


      Subroutine read_Lya_Bph(Lya, bph)

         include "Setting.inc"
         real*8, dimension(start_ydoy_index:end_ydoy_index) :: Lya, bph
         character(len=80) :: line
         integer :: year, doy, i, yyyydoy
         real :: f10_7, f107a, ap, lyman_alpha, beta_ph, factor

         open(unit=101, file=trim(Lya_dir), status='old', action='read')
         read(101, '(A)', iostat=i) line   ! Skip the header line
         do while (.true.)
            read(101, '(A)', iostat=i) line
            if (i /= 0) exit
            read(line, *, iostat=i) year, doy, f10_7, f107a, ap, lyman_alpha, beta_ph
            yyyydoy = year * 1000 + doy
            if (yyyydoy >= start_ydoy_index .and. yyyydoy <= end_ydoy_index) then
!               print *, "Year:", year, "DOY:", doy, "Lyman-alpha:", lyman_alpha
               Lya(yyyydoy) = lyman_alpha
               bph(yyyydoy) = beta_ph
            end if
         end do
         close(101)

         ! Convert line-integrated Lya to line-centered Lya [Emerich et al., 2005]
         factor = (h*c/121.6d-9)*1e11*1e4       !  121.6e-9 m for the wavelength of Lyman-alpha, 1e4 for m2->cm2, and 1e12 from Emmerich et al. (2005)
         Lya = 0.64*(Lya/factor)**1.21        

         return
      End


      Subroutine Make_Parameters_OutFile()
         use Module_for_NVelocityDirection
         use Module_Physics_tag
         include "Setting.inc"
         character*100 filename

         filename = 'MATE_Parameters_' // trim(Runname_in_10char) // '.in'
         open(file=filename,unit=123,status='replace')
         write(123,*) N_vel_directions, nRadial, nEnergy
         write(123,*) nRadial, nLon, nLat_NS, ntperday

         write(123,*) "Above paramters are ..."
         write(123,*) "    [N_vel_directions, nRadial, nEnergy]"
         write(123,*) "    [nRadial, nLon, nLat_NS, ntperday]"
         write(123,*) "Start_Time_in_YYYYDOY = ", Start_Time_in_YYYYDOY
         write(123,*) "End_Time_in_YYYYDOY   = ", End_Time_in_YYYYDOY
         write(123,*) "Output interval       = ", Output_Time_Interval_in_Minute, " [minutes]"
         write(123,*) "Exobase BC:             ", ExobaseBC_Model_Name
         write(123,*) "Physics:                ", tag_phys
         write(123,*) ""

         close(123)


      End