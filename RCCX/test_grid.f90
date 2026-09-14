Program test_grid
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
   print*, 'ip (MLT sectors):', ip
   print*, 'ir (radial steps):', ir
   print*, 'iba min/max:', minval(iba), maxval(iba)
   print*, 'xmlto(1, 1:8):', xmlto(1, 1:8)
   print*, 'xmlto(1, 41:48):', xmlto(1, 41:48)
   print*, 'ro(1, 1), ro(iba(1), 1):', ro(1, 1), ro(iba(1), 1)
End Program test_grid
