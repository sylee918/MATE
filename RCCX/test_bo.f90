Program test_bo
   use constants
   use cread1
   use cread2
   use cgrid
   use cfield
   implicit none
   character(len=80) :: CIMI_flux_file
   integer :: is, it2
   CIMI_flux_file = '/home/sylee/CIMI/output/June2008/2008164_MATE_h.fls'
   is = 1
   do it2 = 1, 70
      call read_CIMI_flux(CIMI_flux_file, it2, is)
   enddo
   print*, 'bo min/max (it=70):', minval(bo), maxval(bo)
   print*, 'fl min/max (it=70):', minval(fl), maxval(fl)
   print*, 'ro min/max (it=70):', minval(ro), maxval(ro)
   print*, 'gridy (sinPA):', gridy
   print*, 'iba (max):', maxval(iba)
End Program test_bo
