   Module EXOBASE_BC

      USE SETTING
      USE MPI_MATE, only: rank
      IMPLICIT NONE

      real*8, dimension(nbx,nby,nbtperday) :: nH_temp, TH_temp
      real*8, dimension(nbx,nby,nbtperday,start_ydoy-nt_bwd_bc:end_ydoy) :: nH_BC, TH_BC
      character*100 filename_BC

      contains

      Subroutine Get_exobaseBC

         integer iday, maxdoy
         character*7 ydoy_str, yearst

         if (ExobaseBC_Model_Name .eq. "CONST") then
            nH_BC = 1.2e5
            TH_BC = 1e3
            return
         endif


         if (start_ydoy/1000 .eq. end_ydoy/1000) then
            write(yearst, '(I4.4)') start_ydoy/1000

            do iday=start_ydoy-nt_bwd_bc,end_ydoy
!               nH_temp = 0.d0 ; TH_temp=0.d0
               write(ydoy_str,'(I7.7)') iday
               write(yearst, '(I4.4)') start_ydoy/1000
               filename_BC = trim(BC_dir) // trim(yearst) // "/" // trim(ExobaseBC_Model_Name) // "_" // trim(ydoy_str) //  ".bc"
               call read_exobaseBC
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
!            nH_temp = 0.d0 ; TH_temp=0.d0
            write(ydoy_str,'(I7.7)') iday
            write(yearst, '(I4.4)') start_ydoy/1000
            filename_BC = trim(BC_dir) // trim(yearst) // "/" // trim(ExobaseBC_Model_Name) // "_" // trim(ydoy_str) //  ".bc"
            call read_exobaseBC()
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
!            nH_temp = 0.d0 ; TH_temp=0.d0
            write(ydoy_str,'(I7.7)') iday
            write(yearst, '(I4.4)') end_ydoy/1000
            filename_BC = trim(BC_dir) // trim(yearst) // "/" // trim(ExobaseBC_Model_Name) // "_" // trim(ydoy_str) //  ".bc"
            call read_exobaseBC()
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


      Subroutine read_exobaseBC

         IMPLICIT NONE

         real, dimension(:,:,:), allocatable :: nH_real, TH_real
         integer nlen, thread_num, IO_unit
         logical iexist

         allocate(nH_real(nbx,nby,nbtperday),TH_real(nbx,nby,nbtperday))

         if (rank .eq. 0) print*, "Read exobase BC file: ", filename_BC

         inquire(file=filename_BC, exist=iexist)
         if (.not. iexist) then
            print*, "File is not exist: ", filename_BC
         else
            inquire(iolength=nlen) nH_real
            nlen=nlen*2

            open(file=filename_BC,newunit=IO_unit,form='unformatted', &
               access='direct',action='read',recl=nlen,status='old')
            read(IO_unit,rec=1) nH_real, TH_real
            close(IO_unit)
         endif

         nH_temp = nH_real*1.d0
         TH_temp = TH_real*1.d0
         deallocate(nH_real,TH_real)

         return
      End

   END MODULE EXOBASE_BC
