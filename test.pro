Pro test

  res=5 

  nRadial=17
  nLat=90/res+1
  nLat2=180/res+1
  nLon=360/res

  dir = '/nobackup/slee122/MATE/0728/'
  f = dir + 'MATE_nH_GRC_RCCX1_2008164.data'

  openu, 1, f   &   readu, 1, nH   &   close, 1

  print, total(nH)
  
  stop


End

