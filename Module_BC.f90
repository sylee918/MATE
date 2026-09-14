   Module EXOBASE_BC

      USE SETTING
      USE TIME_UTILS, only: ydoy_diff_days, ydoy_add_days_int
      USE MPI_MATE, only: rank
      IMPLICIT NONE

      real*8, dimension(nbx,nby,nbtperday) :: nH_temp, TH_temp
      real*8, dimension(:,:,:,:), allocatable :: nH_BC, TH_BC
      integer :: total_bc_days
      character*100 filename_BC

      contains

      Subroutine Get_exobaseBC

         integer :: iday_idx, cur_day
         character*7 :: ydoy_str
         character*4 :: yearst

         total_bc_days = ydoy_diff_days(end_ydoy, BC_Start_Time_in_YYYYDOY) + 1

         if (allocated(nH_BC)) deallocate(nH_BC)
         if (allocated(TH_BC)) deallocate(TH_BC)
         allocate(nH_BC(nbx, nby, nbtperday, total_bc_days))
         allocate(TH_BC(nbx, nby, nbtperday, total_bc_days))

         if (ExobaseBC_Model_Name .eq. "CONST") then
            nH_BC = 1.2e5
            TH_BC = 1e3
            return
         endif

         cur_day = BC_Start_Time_in_YYYYDOY
         do iday_idx = 1, total_bc_days
            write(ydoy_str, '(I7.7)') cur_day
            write(yearst, '(I4.4)') cur_day / 1000
            filename_BC = trim(BC_dir) // trim(yearst) // "/" // trim(ExobaseBC_Model_Name) // "_" // trim(ydoy_str) // ".bc"
            call read_exobaseBC()
            nH_BC(:,:,:,iday_idx) = nH_temp
            TH_BC(:,:,:,iday_idx) = TH_temp
            if (minval(nH_temp) .lt. 1e-15 .or. minval(TH_temp) .lt. 1e-15) then
               print*, 'extra_tools', minval(nH_temp), minval(TH_temp)
               print*, "ERROR: BC has zero or invalid values in: ", trim(filename_BC)
               stop
            endif
            cur_day = ydoy_add_days_int(cur_day, 1)
         enddo

         return
      End Subroutine Get_exobaseBC


      Subroutine read_exobaseBC()

         IMPLICIT NONE

         real, dimension(:,:,:), allocatable :: nH_real, TH_real
         integer :: nlen, IO_unit
         logical :: iexist

         allocate(nH_real(nbx,nby,nbtperday), TH_real(nbx,nby,nbtperday))

         if (rank .eq. 0) print*, "Read exobase BC file: ", trim(filename_BC)

         inquire(file=filename_BC, exist=iexist)
         if (.not. iexist) then
            print*, "ERROR: File does not exist: ", trim(filename_BC)
            stop
         endif

         inquire(iolength=nlen) nH_real
         nlen = nlen * 2

         open(file=filename_BC, newunit=IO_unit, form='unformatted', &
            access='direct', action='read', recl=nlen, status='old')
         read(IO_unit, rec=1) nH_real, TH_real
         close(IO_unit)

         nH_temp = nH_real * 1.d0
         TH_temp = TH_real * 1.d0
         deallocate(nH_real, TH_real)

         return
      End Subroutine read_exobaseBC

   END MODULE EXOBASE_BC
