
Subroutine read_CIMI_flux(CIMI_flux_file, it, is)
   use constants
   use cread1
   use cread2
   use cgrid
   use cfield
   use useless
   character(len=80) :: CIMI_flux_file
   integer :: rc_Re1, ir1, ip1, je1, ig1, i, j, k, m, n
   integer :: it, is, IO_unit
   real :: hour, ro1, xmlto1

   if (it == 1) then
      open(unit=IO_unit,file=trim(CIMI_flux_file),status='old')
      read(IO_unit,*) rc_Re1, ir1, ip1, je1, ig1

      if (ir1.ne.ir.or.ip1.ne.ip.or.je1.ne.je.or.ig1.ne.ig) then
         write(*,*) 'Error: ir, ip, je, ig do not match'
         stop
      endif

      read(IO_unit,*) (varL(i),i=1,ir)
      read(IO_unit,*) (mphi(j),j=1,ip)
      read(IO_unit,*) (gride(is,k),k=1,je)
      read(IO_unit,*) (gridy(m),m=1,ig)
   endif

   ! Read hour, parmod, and Lstar_max
   read(IO_unit,*) hour, parmod, Lstar_max(0)
   ihour = int(hour) + 1

   print*, hour, parmod, Lstar_max(0)
   
   ! Read fluxes @ fixed E & y grids
   do i=1,ir
      do j=1,ip
!         read(11,*) xlati(i,j), xmlt(j), xlatiS(i,j), xmltS(i,j), ro1, xmlto1, &
         read(IO_unit,*) xlati(i,j), xmlt(j), xlatiS(i,j), xmltS(i,j), ro(i,j), xmlto(i,j), &
                    BriN(i,j), BriS(i,j), bo(i,j), iba(j)
         read(IO_unit,*) density(i,j), ompe(i,j), CHpower(i,j), HIpower(i,j), &
                    denWP(n,i,j), TparaWP(n,i,j), TperpWP(n,i,j), &
                    HRPee(n,i,j), HRPii(n,i,j), rppa(j), Lstar(i,j,0), volume(i,j)
         do k=1,je
            read(IO_unit,*) (fl(i,j,k,m), m=1,ig)
         enddo
      enddo
   enddo

   if (it == nt) close(IO_unit)


End Subroutine read_CIMI_flux


Subroutine ChargeExchangeCrossSection(CXsigma)
   use constants
   use cgrid
   real, dimension(ns,je) :: CXsigma
   real, parameter :: p1=4.15, p2=0.531, p3=67.3
   real, parameter :: o1=3.13, o2=0.170, o3=87.5

   ! H-H+ charge exchange cross section from Linsay and Stebing, 2005.
   CXsigma(1,:)= (p1 - p2*log(gride(1,:)))**2 * (1.-exp(-p3/gride(1,:)))**4.5
   CXsigma(2,:)= (o1 - o2*log(gride(2,:)))**2 * (1.-exp(-o3/gride(2,:)))**0.8

   CXsigma = CXsigma*1e-16

End Subroutine ChargeExchangeCrossSection



Subroutine calculate_Parmod
  use constants
  use cread1
  use cread2
  use cfield
  integer IO_unit
  COMMON /GEOPACK1/ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI, &
          SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33, &
          E11,E21,E31,E12,E22,E32,E13,E23,E33

!  Setup parmod0 at tstart 
   thalf=t+0.5*tstep
   if (t.eq.tstart) then
      if (ires.eq.0.and.itype.eq.2) then
         open(newunit=IO_unit,file=trim(outname)//'.le',status='old')
         read(IO_unit,'(a80)') header
         read(IO_unit,*) parmod0
         close(IO_unit)
      else
         call TsyParmod(thalf,tsw,xnswa,vswa,nsw,tdst,Dsta,ndst,timf,byw,bzw,&
                        nimf,imod,parmod0,Dst,DstRC)
      endif
   endif
      
!  Determine parmod 
   parmod(1:10)=parmod0(1:10)
   if (t.gt.tstart.and.ires.gt.0) call TsyParmod(thalf,tsw,xnswa,vswa,nsw, &
                      tdst,Dsta,ndst,timf,byw,bzw,nimf,imod,parmod,Dst,DstRC)
   write(*,'(11f9.3)') t/3600.,parmod

!  Call recalc_08 to calculate the dipole tilt
      vgseX=-400.
      vgseY=0.
      vgseZ=0.
      isec=mod(ifix(thalf),60)
      min1=ifix(thalf)/60
      minu=mod(min1,60)
      ihour1=ifix(thalf)/3600
      ihour=mod(ihour1,24)
      iday2=iday+ifix(thalf)/86400
      call recalc_08(iyear,iday2,ihour,minu,isec,vgseX,vgseY,vgseZ)
      psi1=psi

! Calculate Kp, F107 and Ap
  call lintp(tKp,zkpa,nKp,thalf,zkp)
  call lintp(tF10,F107a,nf10,thalf,F107)
  call lintp(tAp,Api,nAp,thalf,Apt)  



End Subroutine calculate_Parmod



! ************************************************************************
!                             readInputData
! Read parameters: dt, tmax, species, storm, ..., from cimi.dat
! Solar wind, geomagnetic data from *.level file
! Ap, AE data can be found at http://omniweb.gsfc.nasa.gov/form/dx1.html
! ************************************************************************
subroutine readInputData

   use constants
   use cread1
   use cread2
   use tsy_plasma
   use dub_plasma
!   use ModCurvScatt, only: iflc,tflc
!   use ModMate, only: MateFileName,iDoyStartMate,read_mate_data
   integer Kp8(8),iAp(8),iDst24(24),indx(24)
   real bxw1(ndmax),byw1(ndmax),bzw1(ndmax),xnsw1(ndmax),vsw1(ndmax)
   real AE1(ndmax),AL1(ndmax),coef1(6),coef2(2),coef4(4)
   real,allocatable,dimension(:) :: geMLT
   real,allocatable,dimension(:,:) :: geflux,giflux
   character pmz(8)*1,header*80,tmp*3,tmp2*2,tmp4*4
   COMMON /GEOPACK2/G(105),H(105),REC(105)

! open a file to read parameters which control the speices, dt and so on
  open(unit=4,file='cimi.dat',status='old')
  read(4,*) itype              ! itype:1=initialRun,2=continuousRun
  read(4,*) tstart         
  read(4,*) dt
  read(4,*) tmax
  read(4,*) tint               ! time resolution in printing result
  read(4,*) imod,intB          ! imod:0=no Bext,1=t96,2=t04; intB:0=dip,1=IGRF
  read(4,*) ires               ! 0=fixed B, 1=changing B
  read(4,*) tstep,tf1          ! tstep: time updating B; tf1: time smoooth data
  read(4,*) ijs                ! no. of species
  read(4,*) (js(i),i=1,ijs)    ! species: 1=H+, 2=O+, 3=He+, 4==HiE e-
  read(4,*) (init(i),i=1,ijs)  ! ion:0=0,1=data; e-:0=AE8MIN,1=AE8MAX,-ve=+local
  read(4,*) (ibset(i),i=1,ijs) ! at boundary:(99) Maxwellian,(<99) Kappa
  read(4,*) iplsh,icom         ! iPS:1=E-E,2=T-M,4=geo; icom:1=Young,2=1+Pandya
  read(4,*) eplsh              ! ePS: (1)ne=ni,TiTe, (2)RBE, (3)Dubyagin's model
  read(4,*) TiTe               ! Ti/Te @ booundary. TiTe=-1. if RBE bc for e-
  read(4,*) rb                 ! model boundary in RE at equator
  read(4,*) hlosscone          ! loss cone altitude in km
  read(4,*) ipot      ! convection: 1=Weimer, 2=SCE+Hardy, 3=SCE+j||, 4=SCE+Prec
  read(4,*) ihigh,icon         ! hiLatPot:0=Weimer,1=Boyle,2=Hill,3=hybrid; icon
  read(4,*) iain,ExAC          ! iain: potent polar bound, ExAC: Sigma expansion
  read(4,'(1x,a13)') storm
  read(4,'(1x,a13)') outname
  read(4,*) iplsp,tpls         ! 0=no PP, 1=PP by n, 2=PP by gradient; tpls 
  read(4,*) icoul              ! 0=no Coulomb collision, 1=w/ Coulomb collision
  read(4,*) ichor              ! 0=no Chorus, 1=LB_BU low-latitude chorus
  read(4,*) icP      ! ChorusPower: 1=GaussianFit(GF),2=BU,3=Aryan+GF,4=Aryan+BU
  read(4,*) ihiss              ! 0=no Hiss, 1=use D_hiss_BU
  read(4,*) ihP                ! hiss power: 1=Gaussian fit, 3=Aryan
  read(4,*) iEMIC              ! 0=no EMIC, 1=with EMIC
  read(4,*) iEMICdiff          ! 0=no EMIC, 1=with EMIC diffusion
  read(4,*) iflc,tflc          ! 0=no FLC, 1=with FC; tflc=time to apply FLC
  read(4,*) igeo               ! geocorona model: 1=Rairden, 2=Hodges, 3=MATE
!  if (igeo==3) read(4,*) MateFileName    ! MATE file name
!  if (igeo.eq.1) close(4)

  rc=(re_m+hlosscone*1000.)/re_m       ! losscone in Re 

  if (ijs.gt.ns) then
     write(*,*) 'Error: ijs.gt.ns'
     stop
  endif
  iplsp1=0
  if (ichor.ne.0.or.ihiss.ne.0.or.icoul.eq.1) iplsp1=1
  if (iplsp1.gt.0.and.iplsp.eq.0) then
     write(*,*) 'Error: forgot to turn on the plasmasphere'
     stop
  endif
       
  if (iplsp.eq.0) tpls=tmax+1.           ! if iplsp=0, never turnon plasmasphere
  if (iflc.eq.0) tflc=tmax+1.            ! if iflc=0, never turnon FLC 
  if (ipot.eq.1) ihigh=0                 ! All Weimer

! Setup time steps
  nstep=ifix((tmax-tstart)/dt/2.)
  nprint=ifix((tmax-tstart)/tint)+1
  tf2=tf1/2.                ! half time range in sec of smoothing input data

  do n=1,ijs
     if (js(n).eq.1) st2(n)='_h'
     if (js(n).eq.2) st2(n)='_o'
     if (js(n).eq.3) st2(n)='he'
     if (js(n).eq.4) st2(n)='_n'
     if (js(n).eq.5) st2(n)='_e'
  enddo

! Open storm.level to read iyear, iday, f107, Ap, Kp, Dst and SW data
      open(unit=14,file=trim(storm)//'.level',status='old')
      read(14,*) iyear,iday

      ! Read F107 data
      read(14,*) nf10 
      do i=1,nf10 
         read(14,*) iyr,iday1,ihr,F107a(i)
         if (iyr.ne.iyear) call NewYear(iyear,iyr,iday1)
         tF10(i)=(iday1-iday)*86400.+43200.
      enddo

      ! Read Kp, Ap data 
      read(14,*) nline
      nKp=8*nline
      nAp=nKp
      if (nKp.gt.nKpMax) then
         write(*,*) 'Error: nKp.gt.nKpMax'
         stop
      endif
      j=1
      do i=1,nline
         read(14,'(i4,2i2,1x,8(i1,a1),3x,8i3)') iyr,month,idy,Kp8(1), &
             pmz(1),Kp8(2),pmz(2),Kp8(3),pmz(3),Kp8(4),pmz(4),Kp8(5),&
             pmz(5),Kp8(6),pmz(6),Kp8(7),pmz(7),Kp8(8),pmz(8),iAp
         call modd_dayno(iyr,month,idy,iday1,j)
         if (iyr.ne.iyear) call NewYear(iyear,iyr,iday1)
         tp0=(iday1-iday)*86400.
         do k=1,8
            m=(i-1)*8+k
            tKp(m)=tp0+k*10800.-5400.
            dKp=0.
            if (pmz(k).eq.'-') dKp=-0.33
            if (pmz(k).eq.'+') dKp=0.33
            zkpa(m)=float(Kp8(k))+dKp
            tAp(m)=tKp(m)     
            Api(m)=float(iAp(k))
         enddo
      enddo

      ! Read Dst data
      read(14,*) nline
      if (nline.gt.0) ndst=24*nline
      if (nline.le.0) ndst=-1*nline
      if (nline.lt.0) read(14,'(a80)') header
      if (ndst.gt.nDstmax) then
         write(*,*) 'Error: ndst.gt.nDstmax'
         stop
      endif
      j=1
      do i=1,abs(nline)
         if (nline.gt.0) then
            read(14,*) iyr,month,idy,iDst24(1:24)
            call modd_dayno(iyr,month,idy,iday1,j)
            if (iyr.ne.iyear) call NewYear(iyear,iyr,iday1)
            tDst0=(iday1-iday)*86400.
            do k=1,24
               m=(i-1)*24+k
               tDst(m)=tDst0+k*3600.-1800.
               Dsta(m)=float(iDst24(k))
            enddo
         else
            read(14,*) month,idy,ihr,iDst
            call modd_dayno(iyear,month,idy,iday1,j)
            tDst(i)=(iday1-iday)*86400.+ihr*3600.-1800.
            Dsta(i)=float(iDst)
         endif
      enddo

      ! Read and smooth Nsw and Vsw    
      read(14,* ) swlag   ! time in sec for sw travel from s/c to subsolar pt
      read(14,*) nsw
      read(14,'(a80)') header
      j=1
      do i=1,nsw
         read(14,*) idy,month,iyr,ihr,minute,sec,xnsw1(i),vsw1(i)         !ACE
         call modd_dayno(iyr,month,idy,iday1,j)
         if (iyr.ne.iyear) call NewYear(iyear,iyr,iday1)
         tsw(i)=swlag+(iday1-iday)*86400.+ihr*3600.+minute*60.+sec
      enddo
      do i=1,nsw                   ! smooth solar wind data
         tti=tsw(i)-tf2
         ttf=tsw(i)+tf2
         call locate1(tsw,nsw,tti,j1)
         call locate1(tsw,nsw,ttf,j_2)
         j2=j_2+1
         if (j1.eq.0) j1=1
         if (j2.gt.nsw) j2=nsw
         xnswa(i)=0.
         vswa(i)=0.
         do j=j1,j2
            xnswa(i)=xnswa(i)+xnsw1(j)/(j2-j1+1)
            vswa(i)=vswa(i)+vsw1(j)/(j2-j1+1)
         enddo
      enddo
      ! read and smooth IMF data
      read(14,*) nimf
      if (nsw.gt.ndmax.or.nimf.gt.ndmax) then
         print *,'Error: nsw.gt.ndmax.or.nimf.gt.ndmax'
         stop
      endif
      read(14,'(a80)') header
      do i=1,nimf
         read(14,*) idy,month,iyr,ihr,minute,sec,bxw1(i),byw1(i),bzw1(i)
         call modd_dayno(iyr,month,idy,iday1,j)
         if (iyr.ne.iyear) call NewYear(iyear,iyr,iday1)
         timf(i)=swlag+(iday1-iday)*86400.+ihr*3600.+minute*60.+sec
      enddo
      do i=1,nimf                  ! smooth IMF data
         tti=timf(i)-tf2
         ttf=timf(i)+tf2
         call locate1(timf,nimf,tti,j1)
         call locate1(timf,nimf,ttf,j_2)
         j2=j_2+1
         if (j1.eq.0) j1=1
         if (j2.gt.nimf) j2=nimf
         bxw(i)=0.
         byw(i)=0.
         bzw(i)=0.
         do j=j1,j2
            bxw(i)=bxw(i)+bxw1(j)/(j2-j1+1)
            byw(i)=byw(i)+byw1(j)/(j2-j1+1)
            bzw(i)=bzw(i)+bzw1(j)/(j2-j1+1)
         enddo
      enddo

      ! Read and smooth AE, AL data
      read(14,*) nAE
      if (nAE.gt.0) read(14,'(a80)') header
      do i=1,nAE
         read(14,*) iyr,month,idy,ihr,minute,sec,iday1,AE1(i),AU,AL1(i)
         if (iyr.ne.iyear) call NewYear(iyear,iyr,iday1)
         tAE(i)=(iday1-iday)*86400.+ihr*3600.+minute*60.
      enddo
      do i=1,nAE
         tti=tAE(i)-tf2
         ttf=tAE(i)+tf2
         call locate1(tAE,nAE,tti,j1)
         call locate1(tAE,nAE,ttf,j_2)
         j2=j_2+1
         if (j1.eq.0) j1=1
         if (j2.gt.nAE) j2=nAE
         AEa(i)=0.
         ALa(i)=0.
         do j=j1,j2
            AEa(i)=AEa(i)+AE1(j)/(j2-j1+1)
            ALa(i)=ALa(i)+AL1(j)/(j2-j1+1)
         enddo
      enddo
      close(14)


! Find dipole moment, xme
  vgseX=-400.
  vgseY=0.
  vgseZ=0.
  call recalc_08(iyear,iday,i_zero,i_zero,i_zero,vgseX,vgseY,vgseZ)
  DIPMOM=SQRT(G(2)**2+G(3)**2+H(3)**2)   ! DIPMOM in (nT RE^3)
  xme=DIPMOM*re_m**3*1.e-9               ! dipole moment in (T m^3)
!  write(*,*) 'xme ',xme

! MATE exosphere model
!  if (igeo==3) call read_mate_data
!  TimeMate=TimeMate+float(iday-iDoyStartMate)*86400.

! Find location of dipole north in geographic coordinates: elon,ctp,stp
  xmag=0.
  ymag=0.
  zmag=rc
  call geomag_08(xgeo,ygeo,zgeo,xmag,ymag,zmag,m_one)
  ctp=zgeo/rc                  ! cosine of angle between geographic N & dipole N
  stp=sqrt(1.-ctp*ctp)         ! sine of angle between geographic N and dipole N
  elon=atan2(ygeo,xgeo)*180./pi  ! geographic longitude of diple N in degree 
!  write(*,*) 'elon ',elon
!  write(*,*) 'ctp,stp ',ctp,stp 

! Determine whether output from Yu Lin's hybrid code is needed
  ihy=0
  if (imod.eq.3.or.iplsh.eq.3.or.ihigh.eq.3) ihy=1

! FLC or no FLC
!  if (iflc.eq.1) write(*,*) ' FLC scattering is used'
!  if (iflc.eq.0) write(*,*) ' no FLC scattering'

end subroutine readInputData




!-----------------------------------------------------------------------------
  subroutine traceF(imod,intB,xi,yi,zi,dir,dsmax,err,rlim,rmn,parmod,psi, &
                    np,xf,yf,zf,xa,ya,za,ra,ba,npf,iout)
!-----------------------------------------------------------------------------
! Routine does field line tracing from (xi,yi,zi) to rmn
!
! xi,yi,zi,xf,yf,zf,xa,ya,za are in RE and in sm coordinates
! ba is in nT

  implicit none
  external dip_08,IGRF_GSW_08,t04_s,zeroB
  integer,parameter :: i_one=1,m_one=-1
  integer np,np1,npf,iout,m,imod,intB
  real xi,yi,zi,dir,rlim,rmn,parmod(10),psi,xf,yf,zf,rf,xa(np),ya(np),za(np), &
       xg,yg,zg,ra(np),ba(np),xa1(np),ya1(np),za1(np),bxint,byint,bzint, &
       bxext,byext,bzext,bx,by,bz,dsmax,err,rmn1
   
! Initial setup 
  np1=np-1   
  iout=0 
  call smgsw_08(xi,yi,zi,xg,yg,zg,i_one)      ! sm to gsm

   xa1=0. ; ya1=0. ; za1=0.
  
! Start fieldline tracing
     if (imod.eq.0.and.intB.eq.0) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
                rmn,imod,parmod,zeroB,dip_08,xf,yf,zf,xa1,ya1,za1,npf,np)
!     if (imod.eq.1.and.intB.eq.0) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
!                rmn,imod,parmod,t96_01,dip_08,xf,yf,zf,xa1,ya1,za1,npf,np)
     if (imod.eq.2.and.intB.eq.0) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
                rmn,imod,parmod,t04_s,dip_08,xf,yf,zf,xa1,ya1,za1,npf,np)
     if (imod.eq.0.and.intB.eq.1) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
                rmn,imod,parmod,zeroB,IGRF_GSW_08,xf,yf,zf,xa1,ya1,za1,npf,np)
!     if (imod.eq.1.and.intB.eq.1) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
!                rmn,imod,parmod,t96_01,IGRF_GSW_08,xf,yf,zf,xa1,ya1,za1,npf,np)
     if (imod.eq.2.and.intB.eq.1) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
                rmn,imod,parmod,t04_s,IGRF_GSW_08,xf,yf,zf,xa1,ya1,za1,npf,np)
!   print*, 'traceF', xf,yf,zf

! Check iout
  rmn1=rmn+err
  rf=sqrt(xf*xf+yf*yf+zf*zf)
  if (rf.gt.rmn1.or.npf.ge.np1) iout=1

! Calculate ra, ba and convert points to SM if iout=0
  if (iout.eq.0) then
     do m=1,npf
        bxext=0.     
        byext=0.     
        bzext=0.     
        ra(m)=sqrt(xa1(m)*xa1(m)+ya1(m)*ya1(m)+za1(m)*za1(m))
        if (intB.eq.0) call dip_08(xa1(m),ya1(m),za1(m),bxint,byint,bzint)
        if (intB.eq.1) call IGRF_GSW_08(xa1(m),ya1(m),za1(m),bxint,byint,bzint)
!        if (imod.eq.1) call t96_01(imod,parmod,psi,xa1(m),ya1(m),za1(m), &
!                                   bxext,byext,bzext)
        if (imod.eq.2) call t04_s(imod,parmod,psi,xa1(m),ya1(m),za1(m), &
                                   bxext,byext,bzext)
        bx=bxint+bxext
        by=byint+byext
        bz=bzint+bzext
        ba(m)=sqrt(bx*bx+by*by+bz*bz)
        call smgsw_08(xa(m),ya(m),za(m),xa1(m),ya1(m),za1(m),m_one)  ! gsm to sm
     enddo
     ! reset xf,yf,zf in SM coordinates
     xf=xa(npf)
     yf=ya(npf)
     zf=za(npf)
  endif

  end subroutine traceF



!-----------------------------------------------------------------------------
  subroutine zeroB(iopt,parmod,psi,x,y,z,bx,by,bz)
!-----------------------------------------------------------------------------
! A subroutine with the same parameters as t96 and t04 but giving zero B.

  implicit none
  integer iopt
  real parmod(10),psi,x,y,z,bx,by,bz,dummy

! dummy statements to avoid warning at compilation
  dummy=float(iopt)
  dummy=parmod(1)
  dummy=psi
  dummy=x
  dummy=y
  dummy=z

! Zero B
  bx=0.
  by=0.
  bz=0.

  end subroutine zeroB


  !*******************************************************************************
!                             TsyParmod
!  Rountine calculates the parmod in Tsyganenko model.
!*******************************************************************************
      subroutine TsyParmod(thalf,tsw,xnswa,vswa,nsw,tdst,Dsta,ndst,&
                           timf,byw,bzw,nimf,imod,parmod,Dst,DstRC)
  use constants
  real tsw(nsw),xnswa(nsw),vswa(nsw),tdst(ndst),Dsta(ndst),timf(nimf),&
       byw(nimf),bzw(nimf),parmod(10),w04(6),rr(6),xlamb(6),beta1(6),gamm(6)
      
! Parameters for T04_S model
      data rr/0.39,0.7,0.031,0.58,1.15,0.88/     ! relaxation rate in hour^-1
      data xlamb/0.39,0.46,0.39,0.42,0.41,1.29/
      data beta1/0.8,0.18,2.32,1.25,1.6,2.4/
      data gamm/0.87,0.67,1.32,1.29,0.69,0.53/

      parmod(1:10)=0.             ! initial values
      
!  parmod(1): solar wind pressure in nPa     
      call lintp(tsw,xnswa,nsw,thalf,xnsw)
      call lintp(tsw,vswa,nsw,thalf,vsw)
      v2n=xnsw*vsw*vsw
      parmod(1)=xmp*v2n*1.e12/1.e-9    ! Pdyn in nPa
      if (parmod(1).lt.2.0) parmod(1)=2.0  ! set min parmod(1) to 2.0

!  parmod(2): Dst      
      call lintp(tdst,Dsta,ndst,thalf,dst)
      parmod(2)=dst

!  Calculate Dst* (DstRC) from Burton et al.
   c1=0.2        ! constant in nT/(eV cm-3)^0.5
   c2=20.        ! constant in nT
   Psw=0.01*v2n
   DstRC=Dst-c1*sqrt(Psw)+c2

!  parmod(3:4): IMF By, Bz in nT      
      call lintp(timf,byw,nimf,thalf,byimf)
      call lintp(timf,bzw,nimf,thalf,bzimf)
      parmod(3)=byimf
      parmod(4)=bzimf

!  parmod(5:10) for t04_s: w04(1:6) defined in Tsyganenko and Sitnov, 2005
      tti=thalf-100.*3600.              ! 100 hours before thalf 
      call locate1(tsw,nsw,tti,j1)
      if (j1.eq.0) j1=1
      call locate1(tsw,nsw,thalf,j2)
      if (j2.eq.0) j2=1
      w04(1:6)=0.
      do j=j1,j2      ! average over preceding hours
        tk=tsw(j)
        tdiff=(tk-thalf)/3600.   ! time difference in hour
        call lintp(timf,bzw,nimf,tk,bz1)
        if (bz1.lt.0.) Bs1=-bz1
        if (bz1.ge.0.) goto 1       ! +ve Bz, no contribution to w04
        xnsw_n=xnswa(j)/5.          ! normalized sw density
        vsw_n=vswa(j)/400.          ! normalized sw velocity
        Bs_n=Bs1/5.                 ! normalized Bs
        do m=1,6
            ert=exp(rr(m)*tdiff)
            Sk=xnsw_n**xlamb(m)*vsw_n**beta1(m)*Bs_n**gamm(m)
            w04(m)=w04(m)+Sk*ert
        enddo
1           continue
      enddo
      del_t=(tsw(j2)-tsw(j1))/(j2-j1+1)/3600.      ! delta t in hour
      if (del_t.le.0.) del_t=1./12.
      do m=1,6
        w04(m)=w04(m)*rr(m)*del_t
        parmod(m+4)=w04(m)
      enddo


! Set limit to parmod for t04_s
      if (parmod(1).gt.18.) parmod(1)=18.         ! limit solar wind pressure
      if (parmod(2).lt.-300.) parmod(2)=-300.     ! limit Dst
      if (parmod(4).lt.0.) then                   ! limit By, Bz when Bz<0
        Bmax=16.0
        Bmag=sqrt(parmod(3)*parmod(3)+parmod(4)*parmod(4))
        if (Bmag.gt.Bmax) then
            parmod(3)=parmod(3)*Bmax/Bmag
            parmod(4)=parmod(4)*Bmax/Bmag
        endif
      endif
      if (parmod(3).lt.-12.0) parmod(3)=-12.0
      if (parmod(3).gt.12.0) parmod(3)=12.0
      if (parmod(4).lt.-12.) parmod(4)=-12.
      if (parmod(4).gt.12.) parmod(4)=12.
      if (parmod(8).gt.8.0) parmod(8)=8.0         ! partial ring current
      if (parmod(10).gt.100.) parmod(10)=100.     ! region 2 current

      end subroutine TsyParmod


!-----------------------------------------------------------------------
      subroutine lintp(xx,yy,n,x,y)
!-----------------------------------------------------------------------
!  Routine does 1-D interpolation.  xx must be increasing or decreasing
!  monotonically. If x is beyound xx, it will be forced inside xx range.

      implicit none
      integer n,i,j,jl,ju,jm
      real xx(n),yy(n),x,x1,minxx,maxxx,y,d

!  Make sure xx is increasing or decreasing monotonically
      do i=2,n
         if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
            write(*,*) ' lintp: xx is not increasing monotonically '
            write(*,*) n,(xx(j),j=1,n)
            stop
          endif
         if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
            write(*,*) ' lintp: xx is not decreasing monotonically '
            write(*,*) n,(xx(j),j=1,n)
            stop
          endif
      enddo

!  Make sure x is inside xx range
      minxx=minval(xx)
      maxxx=maxval(xx)
      x1=x
      if (x1.lt.minxx) x1=minxx 
      if (x1.gt.maxxx) x1=maxxx 

!    initialize lower and upper values
!
      jl=1
      ju=n
!
!    if not dne compute a midpoint
!
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
!
!    now replace lower or upper limit
!
        if((xx(n).gt.xx(1)).eqv.(x1.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
!
!    try again
!
      go to 10
      endif
!
!    this is j
!
      j=jl      ! if x.le.xx(1) then j=1
!                 if x.gt.xx(j).and.x.le.xx(j+1) then j=j
!                 if x.gt.xx(n) then j=n-1
      d=xx(j+1)-xx(j)
      y=(yy(j)*(xx(j+1)-x1)+yy(j+1)*(x1-xx(j)))/d

      end subroutine lintp


!-------------------------------------------------------------------------------
        subroutine lintp2(x,y,v,nx,ny,x1,y1,v1)
!-------------------------------------------------------------------------------
!  Routine does 2-D interpolation.  x and y must be increasing or decreasing
!  monotonically
!
        implicit none
        integer nx,ny,i,i1,j,j1
        real x(nx),y(ny),v(nx,ny),x1,x2,y1,y2,v1,a,b,q00,q01,q10,q11
        real minx,maxx,miny,maxy

        minx=minval(x)
        maxx=maxval(x)
        x2=x1
        if (x2.lt.minx) x2=minx      ! force x2 inside the x range
        if (x2.gt.maxx) x2=maxx      !
        call locate1(x,nx,x2,i)
        if (i.gt.(nx-1)) i=nx-1      
        if (i.lt.1) i=1               
        i1=i+1
        a=(x2-x(i))/(x(i1)-x(i))

        miny=minval(y)
        maxy=maxval(y)
        y2=y1
        if (y2.lt.miny) y2=miny      ! force y2 inside the y range
        if (y2.gt.maxy) y2=maxy      !
        call locate1(y,ny,y2,j)
        if (j.gt.(ny-1)) j=ny-1      
        if (j.lt.1) j=1               
        j1=j+1
        b=(y2-y(j))/(y(j1)-y(j))

        q00=(1.-a)*(1.-b)
        q01=(1.-a)*b
        q10=a*(1.-b)
        q11=a*b
        v1=q00*v(i,j)+q01*v(i,j1)+q10*v(i1,j)+q11*v(i1,j1)

        end subroutine lintp2


!--------------------------------------------------------------------------
      subroutine locate1(xx,n,x,j)
!--------------------------------------------------------------------------
!  Routine return a value of j such that x is between xx(j) and xx(j+1).
!  xx must be increasing or decreasing monotonically. If not, the locate will
!  stop at the turning point.
!  If xx is increasing:
!     If x=xx(m), j=m-1 so if x=xx(1), j=0  and if x=xx(n), j=n-1
!     If x < xx(1), j=0  and if x > xx(n), j=n
!  If xx is decreasing:
!     If x=xx(m), j=m so if x=xx(1), j=1  and if x=xx(n), j=n
!     If x > xx(1), j=0  and if x < xx(n), j=n

      real xx(n)

!  Make sure xx is increasing or decreasing monotonically
      nn=n
      monoCheck: do i=2,n
         if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
            nn=i-1
            exit monoCheck
         endif
         if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
            nn=i-1
            exit monoCheck
         endif
      enddo monoCheck
      if (nn.ne.n) then
         write(*,*)'locate1: xx is not increasing or decreasing monotonically'
         write(*,*)'n,x ',n,x
         write(*,*)'xx ',xx
         stop
      endif

      jl=0
      ju=nn+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(nn).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 10
      endif
      j=jl

      end subroutine locate1
        

!--------------------------------------------------------------------------
  subroutine NewYear(iyear,iyr,iday1)
!--------------------------------------------------------------------------
! Routine calculates new iday1 at new year transition from iyear to iyr 
! or from iyr to iyear.
! input: iyear,iyr
! input/output: iday1

  implicit none
  integer iyear,iyr,iyear1,iday1,nday0

! determine iyear1
  iyear1=min(iyear,iyr)

! find number of day (nday0) in iyear1
  nday0=365
  if (mod(iyear1,4).eq.0) nday0=366     ! leap year

! Calculate new iday1
  iday1=iday1+(iyr-iyear)*nday0

  end subroutine NewYear



!-----------------------------------------------------------------------------
      subroutine modd_dayno(iyy,imo,idy,iday,j)
!-----------------------------------------------------------------------------
!  Routine finds day number in a year for given month and day of the month,
!  and vice versa
!
!   imo: month
!   idy: day number in the month
!  iday: day number in the year
!
!  When j>0, find day number in a year for given month and day of the month.
!  When j<0, find month and day of the month for given day number in a year.
 
      parameter (nm=12)
      integer imv(nm),imv_r(nm),imv_l(nm)
      data imv_r/0,31,59,90,120,151,181,212,243,273,304,334/ !days in each month
      data imv_l/0,31,60,91,121,152,182,213,244,274,305,335/  ! leap year

!  Determine regular or leap year
      leap=0
      if (mod(iyy,4).eq.0) leap=1
      do i=1,nm
         if (leap.eq.0) then
            imv(i)=imv_r(i)
            iday_max=365
         else
            imv(i)=imv_l(i)
            iday_max=366
         endif
      enddo

!  Find iday when j>0 and imo,idy when j<0
      if (j.ge.0) then
         iday=idy+imv(imo)
      else
         if (iday.gt.iday_max) then       ! year boundary
            iday=iday-iday_max
            iyy=iyy+1
         endif
         call ilocate(imv,nm,iday,imo)
         idy=iday-imv(imo)
      endif

      end subroutine modd_dayno


!--------------------------------------------------------------------------
      subroutine ilocate(ixx,n,ix,j)
!--------------------------------------------------------------------------
!  Routine modified from locate1.  ilocate find the location of an integer
!  in an integer array.

      integer ixx(n)

      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((ixx(n).gt.ixx(1)).eqv.(ix.gt.ixx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 10
      endif
      j=jl

      end subroutine ilocate



   Subroutine write_4D(array4D,iday,is,tag)
      use MATEgrid

      IMPLICIT NONE
      
      real, dimension(nRadial,nLon,nLat,nt) :: array4D
      real, dimension(:,:,:,:), allocatable :: real_array4D
      integer iday, nlen, is, IO_unit
      character(len=7) dayst
      character(len=*), intent(in) :: tag
      character(len=100) filename

      allocate(real_array4D(nRadial,nLon,nLat,nt))
      real_array4D = real(array4D)
      print*, 'maxval(real_array4D)', maxval(real_array4D)

      write(dayst, '(I7.7)') iday
!      filename = trim(outdir) // 'MATE_beta_' // trim(dayst) // '.data'
!      filename = 'RCCX_' // trim(dayst) // '.data'
      if (is==1) filename = trim(tag) // '_p_' // trim(dayst) // '.data'
      if (is==2) filename = trim(tag) // '_o_' // trim(dayst) // '.data'
      inquire(iolength=nlen) real_array4D
      open(file=filename,newunit=IO_unit,form='unformatted',access='direct',recl=nlen,status='replace')
      write(IO_unit,rec=1) real_array4D
      close(IO_unit)

      deallocate(real_array4D)

      return
   End