Program test_xmlto
   use constants
   use cread1
   use cread2
   use cgrid
   use cfield
   character(len=80) :: CIMI_flux_file
   integer :: is, it2, i, j
   is = 1
   CIMI_flux_file = '/home/sylee/CIMI/output/June2008/2008164_MATE_h.fls'
   call readInputData()
   do it2 = 1, 70
      t = real(it2 - 1) * 3600.0
      call calculate_Parmod
      call read_CIMI_flux(CIMI_flux_file, it2, is)
   enddo
   print*, 'xmlto(1:5, 1):', xmlto(1:5, 1)
   print*, 'xmlto(1:5, 2):', xmlto(1:5, 2)
   print*, 'ro(1:5, 1):', ro(1:5, 1)
   print*, 'ro(1:5, 2):', ro(1:5, 2)
End Program test_xmlto
