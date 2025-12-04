!*******************************************************************************
!
!                               cimi.f90  
!
! CIMI - Comprehensive Inner Magnetosphere-Ionosphere model
!
! cimi.f90 is a merge of crcm_06.f90 and rbe_v02.f90. The code calculates 
! ion (0.1-400keV) and electron (1keV-4MeV) fluxes, plasmasphere density,
! region 2 field aligned current and subauroral ionospheric potential.
!
! Magetic field model is t96_01 or t04_s. 
!
! Hardy's model is used for the auroral conductance. Calculated e- precipitation
! is not used to modify auroral conductance.
!  
! Input files : cimi.dat
!               quiet_*.fin            ! * represent species
!               storm.level
!               w2k_ascii.dat          ! if Weimer model is used
!               TSYGANENKO_MUKAI.dat   ! if Tsy-Mukai plasmasheet model is used
!               Hodges_data            ! if Hodges [H] model is used
!               storm.rcmcon1a         ! RCM file, created by setup_cond_rev.f90
!               elecoef.dat            ! RCM file, coeff. for the Hardy model
!          D_LBchorus_BUMLAT010.dat    ! e- LB chorus diff coeffs 
!          D_hissD_hiss_BU.dat         ! e- hiss diffusion coefficients 
!               hybrid.dat             ! data from Yu Lin's hybrid code
! Files needed for continuous run: outname_*_c.f2, outname_*.le, outname_c.Nion
!
! Main output files: outname_*.fls, outname_*.preci, outname.pot, outname.nps 
!                    outname_*.ece
!            
! To compile with double precision: 
!    -o cimi.out ModCurvScatt.f90 cimi.f90 plasmasphere.f90 w2k.f T96_01.f 
!                TS04c.f geopack_2008.f Lstar2.f90 rcm6.f90 gmresm.f90
!
! Created on 23 September 2011 by Mei-Ching Fok, Code 673, NASA GSFC
!
! Modification History
! October 4, 2011
!   - Add chorus and hiss diffusion for electrons
! November 21, 2011
!   - When using T04, limit parmod(3)(IMFBy) and parmod(4)(IMFBz) in case
!     parmod(4) is negative
! December 19, 2011
!   - Add the calculation of precipitation.
! December 21, 2011
!   - When using T04, limit parmod(1) to 20.
!   - Add O+/H+ calculation in subroutine boundary
! December 27, 2011
!   - Set ichdim=4, same as ns
! December 20, 2012
!   - Use AE to parameterize wave powers.
! December 21, 2012
!   - Use finer grid for plasmasphere model. Drive the plasmasphere and wave
!     diffusion few hours later from an empty ion ring current.
! January 3, 2013
!   - Include upper-band chorus diffusion.
! January 31, 2013
!   - Add a subroutine NewYear to handle new year transition.
! February 5, 2013
!   - Use new fitting parameters for LB and UB chorus wave powers.
! February 13, 2013
!   - In subroutine boundary, electron flux at rb is 0 if E < 0.5*gride(1)
! March 27,2013
!   - Add a localized source of particle.
! May 6, 2013
!   - Extend the energy in lower-band chorus diffusion from 10 to 0.1 keV.
! May 14, 2013
!   - Bin energy and particle gain/loss from each process in energy and output
!     results in .ece file
! December 23, 2013
!   - In .ece file, also bin energy and particle gain/loss in 3 bands:
!     inner belt+slot, core of outer belt and outer part of outer belt.
! March 24, 2014
!   - Add the option (iplsh) to use Tsyganenko-Mukai plasmasheet model
! April 25, 2014
!   - Add the calculation of EMIC convective growth rate with bi-Maxwellian
!     and real distribution
! June 19, 2014 
!   - Add the capability to use output of the hybrid code of Yu Lin.
! July 10, 2014
!   - Add the capability of calculating He+ flux
!   - Add option to use alternated Lower-band chorus diffusion coeffs.
! August 18, 2014
!   - Add an option to use D_hiss_UCLA.dat
! September 23, 2014
!   - Add the capability that xlati can be changing with time. 
! January 7, 2015
!   - Use 1/Bo to scale wave diffusion coefficients instead of L^3
! January 20, 2015
!   - Add one more species 'e1', js=4: low-energy (10eV - 40keV) electrons 
! January 21, 2015
!   - Add calculation of EMIC convective growth rate with ring distribution
! September 12, 2016
!   - Add Coulomb collisions 
! December 6, 2016
!   - Add ipot=3 to specific conductances with field-aligned current
!     (Robinson's model) 
! December 20, 2016
!   - Add ipot=4 to specific conductances with electron precipitation
!     (Robinson's model) 
! February 17, 2017
!   - Add George's formula to enhance Pedersen and Hall conductance.
! March 6, 2017
!   - Add George's formula to enhance Pedersen and Hall conductance when E field
!     is large
! May 23, 2017
!   - Identify particle with VarV. VarV=xmm*(K+a)^b. 
! July 18, 2017
!   - Diffusion coefficients take the closest ompe bin instead of interpolation
!   - magnitude of DaE is limit to sqrt(Daa*DEE)
! August 1, 2017
!   - Solve the cimi equation in varL (1/cos(xlati)), instead of xlati
! August 17, 2017
!   - B field below L=3 is constant over time to fix bugs in Tsy model during
!     high geomagnetic activities.
! September 21, 2017
!   - Kappa value is input from cimi.dat - ibset: 99=Maxwellian, <99=Kappa value
! October 24, 2017
!   - Add alternate way to define plasmapause, rppa
! June 14, 2018
!   - Change coordinates from (V,K) to (lnp,lnK)
! November 23, 2018
!   - Calculate birk in cimi instead of in RCM.
! February 1, 2019
!   - Add option for internal field: dipole or IGRF
! September 27, 2019
!   - Add high order scheme to solve advection
! June 11, 2020
!   - Set iplsh=4 for using geosyn data as boundary condition at rb
!   - Generalize varL to 1/(cos(xlati))^n
! December 23, 2020
!   - add an option to use chorus diffusion coefficients and wave amplitude 
!     from Qianli Ma at BU.
! January 27, 2021
!   - The azimuth grid is changed from MLT to MLONG (geomagnetic longitude)
! February 19, 2021
!   - xlati in 2D (xlati(i,j)) and guided by magnetic flux
! April 9, 2021
!   - xlati in 2D (xlati(i,j)) and guided by flux tube volume per unit B flux
! April 21, 2021
!   - add an option (icom) to estimate the O+/H+ and He+/H+ at the outer
!     boundary: 1=Young et al, 2,3=Pandya Bz, Psw
! July 26, 2021
!   - xlati in 2D (xlati(i,j)) and guided by Euler potential alpha 
! February 9, 2022
!   - replace geopack_2005.f by geopack_2008.f
!   - use igrf routines in geopack_2008.f instead of those in IGRF12syn.f90
! February 20, 2022
!   - use D_hiss_BU.dat instead of D_hiss_UCLA.dat. D_hiss_BU.dat covers
!     Fpe/Fce range from 1.2 to 30, while D_hiss_UCLA.dat covers only up to 20.
! September 20, 2022
!   - Add an option (icP ge 3) to use Homayon Aryan's chorus wave model
!     https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JA028403
! October 18, 2022
!   - Add an option (ihP=3) to use Homayon Aryan's hiss wave model
! November 15, 2022
!   - Add an option (eplsh=3) to use Dubyagin's electron plasma sheet model
! February 16, 2023
!   ~ Add fieldline-curvature scattering 
! April 24, 2023
!   - Add the calculation of corotation potential (potentc) and output in *.pot
! April 25, 2023
!   - Include ion precipitation in conductance calculation
!*******************************************************************************

  module constants
      integer,parameter :: i_zero=0,i_one=1,m_one=-1
      real,parameter :: pi=3.14159265358979
      real,parameter :: re_m=6.3712e6	            ! earth's average radius (m)
      real,parameter :: xmp=1.67262192369e-27       ! mass of H+ in kg
      real,parameter :: e_mass=9.1093837015e-31     ! electron mass in kg
      real,parameter :: echarge=1.6e-19	            ! electron charge
      real,parameter :: EM_speed=2.998e8            ! speed of light (m/s)
      real,parameter :: epsilon0=8.8542e-12         ! permittivity of free space
  end module

  module cimigrid_dim
  integer,parameter :: ir=70,ip=48,ik=40,je=24,ig=18,ns=5,irh=20
  integer,parameter :: iw=2*je  
             ! ir=no. of grids in latitude in the ionosphere
             ! ip=no. of magnetic local grids ionosphere
             ! iw = no. of grids in magnetic moment
             ! ik = no. of grids in invariant K
             ! je = no. of fixed energy grids
             ! ig = no. of fixed y grids
             ! ns = no. of species
    parameter (nKpMax=150,nDstmax=400,ndmax=100000) ! max Kp,F107,Dst,SW data pt
    integer ipc,iwc,iph,iwh            ! dimension of LB chorus & hiss diff coef
  end module

  module cread1
        use cimigrid_dim
        character outname*13,storm*13,st2(ns)*2
        real rc,xme
  end module

  module cread2
        use cimigrid_dim
        integer iyear,iday,ibset(ns),ipot,ihigh,icon,iain,ijs,js(ns),itype, &
                nstep,nprint,imod,intB,ires,init(ns),ichor,ihiss,iEMIC,iEMICdiff,&
                icP,ihP,iplsp,igeo,icoul,iplsh,eplsh,icom,ndst,nf10,nAp,nKp, &
                nAE,nsw,nimf,ihy
        integer ntg,nMLTg                        ! dimension of geosyn flux
        real,allocatable,dimension(:,:,:) :: eflux66,iflux66   ! arrays for
        real,allocatable,dimension(:,:) :: MLTgeo              ! geosyn flux
        real,allocatable,dimension(:) :: tgeo                  !
        real tstart,tmax,dt,tint,tstep,TiTe,rb,hlosscone,F107a(nKpMax), &
             tF10(nKpMax),Api(nKpMax),zkpa(nKpMax),tKp(nKpMax),Dsta(nDstmax), &
             tAp(nKpMax),tdst(nDstmax),AEa(ndmax),ALa(ndmax),tAE(ndmax), &
             xnswa(ndmax),vswa(ndmax),tsw(ndmax),bxw(ndmax),byw(ndmax), &
             bzw(ndmax),timf(ndmax),tpls,ExAC,elon,ctp,stp
  end module

  module cgrid
      use cimigrid_dim
      real xlati(ir,ip),dlati(ir+1,ip),xlath(irh),xlatd(0:ir+1,ip),varL(ir+1), &
           ksai(ir+1,ip),xkb(ik),mphi(ip),xmass(ns),xk(0:ik+1),&
           gride(ns,je),ebound(ns,0:je),gridp(ns,0:iw+1),lnp(ns,0:iw+1),dlnk, &
           dlnp(ns),gridy(ig),cosSq(ig),sinSq(ig),xmass1(5),gride_e(je), &
           d4(ns),dvarL,ekev(ns,0:iw+1),dkeV(ns,iw),vel(ns,0:iw+1), &
           pcEo(ns,iw),potentc(ir,ip),dphi,mlon(ip),gride_i(je)
  end module

  module cfield
    use cimigrid_dim
    integer iba(ip),iday2,ihour,jnoon,mcN(ir,ip),mcS(ir,ip)
    real Dst,DstRC,bo(ir,ip),ro(ir,ip),xmlto(ir,ip),xo(ir,ip),yo(ir,ip), &
         volume(ir,ip),y(ir,ip,0:ik+1),Hdens(ir,ip,ik),xL1,xL2,psi1, &
         xkcN(ir,ip),fcone(ns,ir,ip,iw,ik),bm(ir,ip,0:ik+1),dmu(ir,ip,ik), &
         xkcS(ir,ip),xlatiS(ir,ip),phiS(ir,ip),Tbounce(ns,ir,ip,iw,ik), &
         xjac(ns,ir+1,ip,0:iw+1,0:ik+1),dlnBmdL(ir,ip,ik),dlnBmdp(ir,ip,ik), &
         dlnBmdL1(ir,ip,ik),dlnBmdp1(ir,ip,ik),lnbm(ir,ip,ik),parmod0(10), &
         phi(ip),xmlt(ip),rsb(ip),xmltb(ip),bob(ip),tya(ir,ip,0:ik+1), &
         BiN(ir,ip),BiS(ir,ip),parmod(10),mlonS(ir,ip),sini_f(ir,ip), &
         BriN(ir,ip),BriS(ir,ip),zkp,F107,Apt,dsmax,err
  end module

  module ccepara
        use cimigrid_dim
        real ceSigma(ns,iw,2),achar(ns,ir,ip,iw,ik)
  end module

  module cPlasmasphere
     use cimigrid_dim
     real,parameter :: Tcold=1.       ! plasmasphere temperature in eV
     real,parameter :: densityP=1.e8  ! density to define plasmapause
     real,parameter :: Rcold(4)=[0.77,0.03,0.2,1.] !plasmasphere H,O,He,e
     real rppa(ip)                    ! plasmapause radial distance
     real density(ir,ip)              ! cold plasma densities in m^-3
  end module

  module cPlasmasphere_new
     integer,parameter :: nlp=209,npp=192    ! grid size for plasmasphere
     integer ibp(npp)
     real Nion(nlp,npp) ! Nion=ion per magnetic flux
     real xlatp(nlp,npp),phip(npp),rp(nlp,npp),Brip(nlp,npp)
     real xlatpS(nlp,npp),phipS(npp),volp(nlp,npp)
     real varLp(nlp),mphip(npp),vlEp(nlp,npp),vpEp(nlp,npp),ksaip(nlp,npp)
  end module

  module cCoulpara
     use cimigrid_dim     ! below: velocity of 1eV H+,O+,He+,e-
     real,parameter :: coullog=21.5      ! Coulomb logarithm       
     real coulii(ns,0:iw),coulee(ns,0:iw)
  end module

  module convect
        use cimigrid_dim
        real birk_f(ir,ip),potent(ir,ip),potenth(0:irh,ip),AE,AL,SWBz,SWVel
  end module

  module cVdrift
       use cimigrid_dim
       integer ihol,ihop      ! recommend 7 for diople & 5 for IGRF
       real vlB(ns,ir,ip,iw,ik),vpB(ns,ir,ip,iw,ik)
       real vlE(ir,ip),vpE(ir,ip),vlnp(ir,ip,ik)
  end module

  module cInterFlux
        use cimigrid_dim
        real fbL0(ip),fbL1(ip),cl(0:ir,ip),cp(ir,ip), &
             faL(0:ir,ip),fap(ir,ip),fupL(0:ir,ip),fupp(ir,ip) 
  end module

  module cinitial
        use cimigrid_dim
        integer ib0(ip)
        real f2(ns,ir,ip,iw,ik),esum(ns,ir,ip,je+2),psum(ns,ir,ip,je+2), &
             rbsum(ns),rcsum(ns),eout(ns,2,3,je+2), &
             esum3(ns,3,ip,je+2),psum3(ns,3,ip,je+2),E2D(ns,ir,ip), &
             Lstar_max(0:ik),Lstar(ir,ip,0:ik)
  end module

  module closs
        use cimigrid_dim
        real xled(ns,3,ip,je+2),xlee(ns,3,ip,je+2),xlel(ns,3,ip,je+2), &
             xlec(ns,3,ip,je+2),xleb(ns,3,ip,je+2), &
             pled(ns,3,ip,je+2),plee(ns,3,ip,je+2),plel(ns,3,ip,je+2), &
             plec(ns,3,ip,je+2),pleb(ns,3,ip,je+2), &
             xPreN(ns,ir,ip,je+2),pPreN(ns,ir,ip,je+2),xPre0N(ns,ir,ip,je+2), &
             pPre0N(ns,ir,ip,je+2),xPreS(ns,ir,ip,je+2),pPreS(ns,ir,ip,je+2), &
             xPre0S(ns,ir,ip,je+2),pPre0S(ns,ir,ip,je+2), &
             preFe(ir,ip,je),prePe(ir,ip,je),preFi(ir,ip,je),prePi(ir,ip,je)
        real HRPe(ns,ir,ip),HRPi(ns,ir,ip),HRPe0(ns,ir,ip),HRPi0(ns,ir,ip)
  end module

  module cbound
        use cimigrid_dim
        integer ntg0
        real parE(ns,ip),perE(ns,ip),dena(ns,ip),fb(ns,ip,iw,ik)  ! fb is psd
  end module

  module waveDiffCoef
     use cimigrid_dim
     integer,parameter :: iDaa=1,iDpp=1,iDap=0
     integer,parameter :: kLat=3 
     integer ipa,jpa,iLc,iLh
     real,allocatable,dimension(:) :: cOmpe,hOmpe,cPA,hPA,xLc,xLh,BLc,BLh
     real,allocatable,dimension(:,:,:,:,:) :: cDaa,cDpp,cDap
     real,allocatable,dimension(:,:,:,:) :: hDaa,hDpp,hDap
     ! variables for EMIC waves
     integer :: ipEmic,&   ! size of fpe/fce array for EMIC wave diffusion   
                iePA       ! size of energy arrary for "
     real,allocatable,dimension(:) :: &
                EmicOmpe,& ! fpe/fce array
                BEmic,&    ! EMIC wave band array (H+,He+,O+)
                ePA        ! pitch-angle array
     real,allocatable,dimension(:,:,:,:) :: &
          EmicHDaa,EmicHDpp,EmicHDap,&    ! H-band EMIC waves' Daa, Dap, Dpp
          EmicHeDaa,EmicHeDpp,EmicHeDap,& ! He-band
          EmicODaa,EmicODpp,EmicODap      ! O-band
                    !dimension: ns, fpe/fce, iw, ik  
  end module

  module cWpower
     use cimigrid_dim
     use waveDiffCoef, only: kLat
     integer kAE,kBs,kVs,kLA,kLB,kltA,kltB
     real CHpower(ir,ip),HIpower(ir,ip),ompe(ir,ip),CHpowerL(ir,ip,kLat)
     real Cpower0,Hpower0
     real,allocatable,dimension(:,:,:,:,:,:) :: CpowerAr      ! from Homayan
     real,allocatable,dimension(:,:,:,:,:) :: HpowerAr     ! Homayan's hiss
     real,allocatable,dimension(:,:,:,:) :: CpowerBU
     real,allocatable,dimension(:) :: LcAr,MLTAr,LcBU,MLTBU,Lhh,MLTh
     ! EMIC waves amplitude (nT) 
     !  inside plasmasphere
     real EmicPower0
     real EmicHPowerInPs(ir,ip),&  ! H band
          EmicHePowerInPs(ir,ip),& ! He band
          EmicOPowerInPs(ir,ip),&  ! O band
     !  outside plasmasphere
          EmicHPowerOutPs(ir,ip),&  ! H band
          EmicHePowerOutPs(ir,ip),& ! He band
          EmicOPowerOutPs(ir,ip)    ! O band
  end module

  module WaveGrowth       ! EMIC growth
         use cimigrid_dim
         integer,parameter :: nw=100 ! number of normalized frequency, Xpar
         real Xpar(nw)               ! normalized freq by H+ gyro-freq
         real denWP(ns,ir,ip),TparaWP(ns,ir,ip),TperpWP(ns,ir,ip) ! m^-3,keV
         real fWP(ns,ir,ip,iw,ik)    ! psd in m^-6 s^3
         real Vpar(ns,ir,ip,iw,ik),Vper(ns,ir,ip,iw)
  end module

  module crtp
     use cimigrid_dim, only: ir,ip
     integer,parameter :: nth=59
     real Rrtp(ir),Trtp(nth),Prtp(ip+1)
  end module

  module cStDiTime
         use cimigrid_dim
         real SDtime(ir,ip,iw,ik)
  end module

  module tsy_plasma
         integer, parameter :: n_points_tsy=16
         real a_tsy(6,n_points_tsy)
  end module

  module dub_plasma
         use cimigrid_dim
         integer, parameter :: n_points_dub=9
         integer, parameter :: n_const_dub=8
         real a_dub(4,n_points_dub),t_dub(2,n_const_dub)
         real BsAveNa(ndmax),BsAveTa(ndmax),BnAvea(ndmax)
         real NswAvea(ndmax),VswAvea(ndmax)
  end module

  module hodges_data
       integer, parameter :: nr=40 ! Number of radial distances 
       ! Alm and Blm are coeffs from the table; Ylm's are harmonic expansion
       real,dimension(:,:,:) :: Alm(0:3,0:3,1:nr),Blm(0:3,0:3,1:nr),Ylm(0:3,0:3)
       character(40) :: hodges_file
       character(6) :: season
       real, dimension(:) :: tabular_rad(1:nr),tabular_n(1:nr),rad_deriv(1:nr-1)
  end module

  module hybrid
         integer :: nxp,nyp,nxb,nyb,nzb,mih,nih
         real,pointer,dimension(:) :: xpp,ypp,xbb,ybb,zbb,thetah,phih
         real,pointer,dimension(:,:) :: potenthy
         real,pointer,dimension(:,:,:) :: bxx,byy,bzz,eqdata
  end module

! Modules for rcm6
      module rcmgrid_dim
             use cimigrid_dim
             parameter(jdimf=ip ,idim=ir,jdim=ip+3,ncoeff=5,signbe=1.)
             parameter(iwdim=iw,ikdim=ik,ichdim=ns)
             parameter (n_kp=6)       ! Hardy model: kp from 0 to n_kp
      end module
      module rcmgrid
             use rcmgrid_dim,only: idim,jdim,jdimf,ncoeff,signbe,iwdim,ikdim, &
                                   ichdim,n_kp
             real dlam,dpsi,ri,re,offset,alpha(idim,jdim),beta(idim,jdim),&
                  colat(idim,jdim),aloct(idim,jdim),sini(idim,jdim),&
                  vcorot(idim,jdim),bir(idim,jdim),ain(jdim),&
                  xmin(idim,jdim),ymin(idim,jdim)
             integer ibndy(jdim),imin,imax,jmax,ii1,ii2,jj1,jj2,iint,jint, &
                     jwrap,jdif
      end module
      module plasma
             use rcmgrid_dim,only: idim,jdim
             real, allocatable, dimension(:,:,:,:,:) :: wk,eta
             real vm(idim,jdim)
      end module
      module potential
             use rcmgrid_dim,only: idim,jdim
             real vrcm(idim,jdim),vpob(jdim)
      end module
      module coefficients
             use rcmgrid_dim,only: idim,jdim,ncoeff
             real c(ncoeff,idim,jdim)
      end module
      module conductance
             use rcmgrid_dim,only: idim,jdim
             real pedpsi(idim,jdim),pedlam(idim,jdim),hall(idim,jdim),&
                  qtped(idim,jdim),qtplam(idim,jdim),qthall(idim,jdim),ss(jdim)
      end module
      module qtconductance
             integer idim1,jdim1
             real,allocatable,dimension(:,:) :: qtped0,qtplam0,qthall0
             real,allocatable,dimension(:) :: ss0,colat0,aloct0
      end module
      module current
             use rcmgrid_dim,only: idim,jdim
             real birk(idim,jdim)
      end module

! Main program
      program cimi   
      use constants
      use cread1
      use cread2
      use cgrid
      use cfield
      use ccepara
      use cPlasmasphere
      use cPlasmasphere_new
      use cCoulpara
      use convect
      use cVdrift
      use cinitial
      use closs   
      use cbound  
      use waveDiffCoef
      use cWpower     
      use WaveGrowth
      use crtp
      use cStDiTime
      use hodges_data
      use potential
      use hybrid
      use ModCurvScatt, only: iflc,init_flc

! Initial setup 
      vrcm=0.            ! initial values for RCM potenitals (Volt)
      call readInputData
      t=tstart
      ihol=7                 ! recommend 7 for diople
      if (intB.eq.1) ihol=5  ! recommend 5 for IGRF
      ihop=ihol
      call grids(ijs,js)
      write(*,'(5x,a13)') outname
      if (ipot.ge.2) call qt_conductance
      if (ihy.eq.1) call hybridsetup
      call init_flc(ns,ir,ip,iw,ik,echarge,re_m)
      call fieldpara(t)
      if (ipot.ge.2) call RCMsetup(t,xlati,xmlt)
      call find_ceSigma
      call cepara
      if (icoul.eq.1) call coulpara
      call StDiTime(ijs,js,dt,volume,xlati,rc,re_m,xme,iba)
      call initial(itype,t,outname,st2,rb,gride,xmass,js,ijs,init)
      write(*,*)'rbsum ',sum(rbsum(1:ijs))
      call convection(t)
      if (iplsh.le.2) call boundary(t,f2)
      if (iplsh.eq.4) call geoboundary(t,f2)
      call VdriftE(potent,ijs,js)
      call VdriftB(ijs,js)

! Setup for the plasmasphere calculation 
  rppa(:)=0.             ! iniriL plasmapause location
  density(:,:)=0.        ! initial plasmasphere density
  Nion(:,:)=0.           ! start w/ an empty plasmasphere
  initPL=0               ! before plasmasphere initialization
  if (iplsp.ne.0) then
     call setgrid(nlp,npp,varLp,mphip)
     if (t.ge.tpls) then
        initPL=1
        call mapgrid
        if (itype.eq.1) then    ! start with a saturated plasmasphere
           delt=10.*86400.
           write(*,*) 'delt(day) ',delt/86400.
           call mapPotent
           call plasmasphere_new(iyear,iday2,ihour,rc,delt)
        else
           open(unit=31,file=trim(outname)//'_c.Nion',status='old', &
                form='unformatted')
           read(31) Nion
           close(31)
        endif
        call get_density(nlp,npp,Nion,volp,varLp,mphip,ir,ip,varL,mphi)
     endif
  endif

! Read diffusion coeffs, BU wave pwoer, and calculate wave power 
      iwave=0
      if (ichor.ne.0.or.ihiss.ne.0.or.iEMICdiff) iwave=1
      if (iwave.eq.1) write(*,*) 'iDaa,iDpp,iDap ',iDaa,iDpp,iDap
      if (ichor.ne.0) write(*,*) 'klat ',klat
      call readDiffCoef(ichor,ihiss)
      CHpower(:,:)=0.
      HIpower(:,:)=0.
      ompe(:,:)=0.           
      if (iwave.eq.1) then
         call ReadWavePower
         call WavePower(density,Bo,ro,xmlto,iba,icP,ihP,ichor,ihiss)
      endif

! Print initial results if itype=1.
      if (itype.eq.1) then
         call p_result(t,tstart,parE,perE,dena,ompe,CHpower, &
                       HIpower,nprint,js,ijs,ibset,tpls,ipot)
         call p_rtp(t)
      endif

! Start the calculation, the time loop
  timeloop: do i3=1,nstep

     call Edrift(t,dt,f2,eout,ib0,ijs)
     call Bdrift(t,dt,f2,eout,ib0,ijs)
     call sum3band(xL1,xL2,xled,pled)

     if (initPL.eq.1.and.icoul.eq.1) then
        call energyDep(HRPe,m_one)  ! Calculate RC energy before Coulomb Col
        call Coulomb(f2,coulee,density) ! RC Coulomb interact w plasmasphere e-
        call energyDep(HRPe,i_one)  ! find RC heat rate to plasmasphere e-
        call Coulomb(f2,coulii,density) ! RC Coulomb interact w plasmasphere ion
        call energyDep(HRPi,i_one)  ! find RC heat rate to plasmasphere ions
        call sum3band(xL1,xL2,xlec,plec)
     endif

     call charexchange(f2,achar,iba,js,ijs)
     call StrongDiff(f2,iba,js,ijs)
     if (initPL.eq.1.and.iwave.eq.1) call diffuse_VLF(f2) 
     if (initPL.eq.1.and.iEMICdiff.eq.1) then
        call EmicPower(t)
        call diffuse_UV_ions 
     endif
     if (iflc.gt.1) call diffuse_flc(f2) 
     call sum3band(xL1,xL2,xlee,plee)
  
     call sume(xPreN,pPreN,ijs,m_one)   
     call lossconeN(f2,ijs)
     call sume(xPreN,pPreN,ijs,i_one)   
     call sume(xPreS,pPreS,ijs,m_one)   
     call lossconeS(f2,ijs)
     call sume(xPreS,pPreS,ijs,i_one)   
     call sum3band(xL1,xL2,xlel,plel)

     t=t+dt

     ! Update E field, drift velocity and boundary distributions
     call convection(t)
     call VdriftE(potent,ijs,js)
     if (iplsh.le.2) call boundary(t,f2)
     if (iplsh.eq.4) call geoboundary(t,f2)

     if (iflc.gt.1) call diffuse_flc(f2) 
     if (initPL.eq.1.and.iEMICdiff.eq.1) then
        call EmicPower(t)
        call diffuse_UV_ions 
     endif
     if (initPL.eq.1.and.iwave.eq.1) call diffuse_VLF(f2)
     call StrongDiff(f2,iba,js,ijs)
     call charexchange(f2,achar,iba,js,ijs)
     call sum3band(xL1,xL2,xlee,plee)

     call sume(xPreS,pPreS,ijs,m_one)   
     call lossconeS(f2,ijs)
     call sume(xPreS,pPreS,ijs,i_one)   
     call sume(xPreN,pPreN,ijs,m_one)   
     call lossconeN(f2,ijs)
     call sume(xPreN,pPreN,ijs,i_one)   
     call sum3band(xL1,xL2,xlel,plel)

     if (initPL.eq.1.and.icoul.eq.1) then
        call energyDep(HRPi,m_one)  ! Calculate RC energy before Coulomb Col
        call Coulomb(f2,coulii,density) ! RC Coulomb interact w plasmasphere ion
        call energyDep(HRPi,i_one)  ! find RC heat rate to plasmasphere ions
        call Coulomb(f2,coulee,density) ! RC Coulomb interact w plasmasphere e-
        call energyDep(HRPe,i_one)  ! find RC heat rate to plasmasphere e-
        call sum3band(xL1,xL2,xlec,plec)
     endif

     call Bdrift(t,dt,f2,eout,ib0,ijs)
     call Edrift(t,dt,f2,eout,ib0,ijs)
     call sum3band(xL1,xL2,xled,pled)

     t=t+dt

     ! Update B field every tstep seconds
     if (mod(t,tstep).eq.0.) then
        call fieldpara(t)
        call sum3band(xL1,xL2,xleb,pleb)
        write(*,*)'rbsum ',sum(rbsum(1:ijs))
        if (ipot.ge.2) call RCMsetup(t,xlati,xmlt)
        call cepara
        call VdriftB(ijs,js)
        call StDiTime(ijs,js,dt,volume,xlati,rc,re_m,xme,iba)
        if (t.ge.tpls.and.iplsp.ne.0) call mapgrid
     endif

     ! Update E field, drift velocity and boundary distributions
     call convection(t)
     call VdriftE(potent,ijs,js)
     if (iplsh.le.2) call boundary(t,f2)
     if (iplsh.eq.4) call geoboundary(t,f2)

     ! update the plasmasphere density
     if (t.ge.tpls.and.iplsp.ne.0) then
        delt=2.*dt
        if (initPL.eq.0) then
           delt=10.*86400.
           initPL=1
           write(*,*) 'delt(day) ',delt/86400.
        endif
        call mapPotent
        call plasmasphere_new(iyear,iday2,ihour,rc,delt)
        call get_density(nlp,npp,Nion,volp,varLp,mphip,ir,ip,varL,mphi)
     endif

     ! Update wave power at t.ge.tpls
     if (t.ge.tpls.and.iwave.eq.1) &
        call WavePower(density,Bo,ro,xmlto,iba,icP,ihP,ichor,ihiss)

     ! print result at every tint sec
     if (mod(t,tint).eq.0.) then
         call p_result(t,tstart,parE,perE,dena,ompe,CHpower, &
                       HIpower,nprint,js,ijs,ibset,tpls,ipot)
         call p_rtp(t)
     endif

     if (t.ge.tmax) exit timeloop

  enddo timeloop     ! end of the time loop

  end


! **************************  end of main  *****************************



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
      use ModCurvScatt, only: iflc,tflc
      use ModMate, only: MateFileName,iDoyStartMate,read_mate_data
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
  if (igeo==3) read(4,*) MateFileName    ! MATE file name
  if (igeo.eq.1) close(4)

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

! Read data for Tsyganenko-Mukai plasmasheet model
   a_tsy(:,:)=0.
   if (iplsh.eq.2) then
      open(14,file='TSYGANENKO_MUKAI.dat',status='old')
      do i=1,n_points_tsy
         read(14,*) tmp,coef1
         a_tsy(1:6,i)=coef1(1:6)
      enddo
      close(14)
   endif

! Read data for Dubyagin et al. 2016 electron plasmasheet model and calculate
! averaged IMF, Nsw and Vsw.
   a_dub(:,:)=0.
   t_dub(:,:)=0.
   if (eplsh.eq.3) then
      open(15,file='dubyagin_2016.dat',status='old')
      do i=1,n_points_dub
         read(15,*) tmp2,coef4
         a_dub(1:4,i)=coef4(1:4)
      enddo
      do i=1,n_const_dub
         read(15,*) tmp4,coef2
         t_dub(1:2,i)=coef2(1:2)
      enddo
      t_dub(:,:)=t_dub(:,:)*3600.     ! convert hour to second
      close(15)
      ! Calculate average Nsw
      do i=1,nsw
         ttf=tsw(i)-t_dub(1,1)
         tti=ttf-t_dub(1,2)
         call locate1(tsw,nsw,tti,j1)
         call locate1(tsw,nsw,ttf,j_2)
         j2=j_2+1
         if (j1.eq.0) j1=1
         if (j2.gt.nsw) j2=nsw
         NswAvea(i)=0.
         do j=j1,j2
            NswAvea(i)=NswAvea(i)+xnsw1(j)
         enddo
         NswAvea(i)=NswAvea(i)/(j2-j1+1)
      enddo
      ! Calculate average BsN
      do i=1,nimf
         ttf=timf(i)-t_dub(1,3)
         tti=ttf-t_dub(1,4)
         call locate1(timf,nimf,tti,j1)
         call locate1(timf,nimf,ttf,j_2)
         j2=j_2+1
         if (j1.eq.0) j1=1
         if (j2.gt.nimf) j2=nimf
         BsAveNa(i)=0.
         do j=j1,j2
            Bs1=-bzw1(j)
            if (Bs1.lt.0.) Bs1=0.
            BsAveNa(i)=BsAveNa(i)+Bs1     
         enddo
         BsAveNa(i)=BsAveNa(i)/(j2-j1+1)
      enddo
      ! Calculate average BsT
      do i=1,nimf
         ttf=timf(i)-t_dub(2,3)
         tti=ttf-t_dub(2,4)
         call locate1(timf,nimf,tti,j1)
         call locate1(timf,nimf,ttf,j_2)
         j2=j_2+1
         if (j1.eq.0) j1=1
         if (j2.gt.nimf) j2=nimf
         BsAveTa(i)=0.
         do j=j1,j2
            Bs1=-bzw1(j)
            if (Bs1.lt.0.) Bs1=0.
            BsAveTa(i)=BsAveTa(i)+Bs1     
         enddo
         BsAveTa(i)=BsAveTa(i)/(j2-j1+1)
      enddo
      ! Calculate average Vsw
      do i=1,nsw
         ttf=tsw(i)-t_dub(2,5)
         tti=ttf-t_dub(2,6)
         call locate1(tsw,nsw,tti,j1)
         call locate1(tsw,nsw,ttf,j_2)
         j2=j_2+1
         if (j1.eq.0) j1=1
         if (j2.gt.nsw) j2=nsw
         VswAvea(i)=0.
         do j=j1,j2
            VswAvea(i)=VswAvea(i)+vsw1(j)
         enddo
         VswAvea(i)=VswAvea(i)/(j2-j1+1)
      enddo
      ! Calculate average Bn
      do i=1,nimf
         ttf=timf(i)-t_dub(2,7)
         tti=ttf-t_dub(2,8)
         call locate1(timf,nimf,tti,j1)
         call locate1(timf,nimf,ttf,j_2)
         j2=j_2+1
         if (j1.eq.0) j1=1
         if (j2.gt.nimf) j2=nimf
         BnAvea(i)=0.
         do j=j1,j2
            Bn1=bzw1(j)
            if (Bn1.lt.0.) Bn1=0.
            BnAvea(i)=BnAvea(i)+Bn1
         enddo
         BnAvea(i)=BnAvea(i)/(j2-j1+1)
      enddo
   endif

! Read geosyn data if iplsh eq 4. Make sure electron end ion energy bins are 
! consistent with gride for electrons and ions.
  if (iplsh.eq.4) then
     open(unit=18,file=trim(storm)//'.geosyn',status='old')
     do i=1,14
        read(18,'(a80)') header
     enddo
     read(18,*) ntg           ! number of geosyn data block
     read(18,*) nMLTg         ! number of MLT/spacecraft
     allocate (tgeo(ntg),MLTgeo(ntg,nMLTg))
     allocate (eflux66(ntg,nMLTg,iw),iflux66(ntg,nMLTg,iw))
     allocate (geMLT(nMLTg),geflux(nMLTg,je),giflux(nMLTg,je))
     read(18,'(a80)') header
     j=1
     do i=1,ntg
        read(18,*) iyr,month,idy,ihr,minute,isec
        call modd_dayno(iyr,month,idy,iday1,j)
        if (iyr.ne.iyear) call NewYear(iyear,iyr,iday1)
        tgeo(i)=(iday1-iday)*86400.+ihr*3600.+minute*60.+isec*1.
        do m=1,nMLTg
           read(18,*) geMLT(m)
           read(18,*) geflux(m,1:je)    ! flux in 1/s-cm2-sr-keV
           read(18,*) giflux(m,1:je)    ! flux in 1/s-cm2-sr-keV
           ! Remove -ve flux
           do k=1,je
              if (geflux(m,k).lt.1.e-30) geflux(m,k)=1.e-30
              if (giflux(m,k).lt.1.e-30) giflux(m,k)=1.e-30
           enddo
        enddo
        ! arrange flux in accending local time and expand at ekev
        call indexx(nMLTg,geMLT,indx)
        do m=1,nMLTg
           m1=indx(m)
           MLTgeo(i,m)=geMLT(m1)
           do k=1,je-1
              k2=k*2
              eflux66(i,m,k2)=(geflux(m1,k)**0.75)*(geflux(m1,k+1)**0.25)
              eflux66(i,m,k2+1)=(geflux(m1,k)**0.25)*(geflux(m1,k+1)**0.75)
              iflux66(i,m,k2)=(giflux(m1,k)**0.75)*(giflux(m1,k+1)**0.25)
              iflux66(i,m,k2+1)=(giflux(m1,k)**0.25)*(giflux(m1,k+1)**0.75)
           enddo
           eflux66(i,m,1)=geflux(m1,1)
           eflux66(i,m,iw)=geflux(m1,je)
           iflux66(i,m,1)=giflux(m1,1)
           iflux66(i,m,iw)=giflux(m1,je)
        enddo  
     enddo
     close(18)
  endif

! Read data files for Hodges geocorona model
  if (igeo.eq.2) then
     call hodges_init
  endif

! Find dipole moment, xme
  vgseX=-400.
  vgseY=0.
  vgseZ=0.
  call recalc_08(iyear,iday,i_zero,i_zero,i_zero,vgseX,vgseY,vgseZ)
  DIPMOM=SQRT(G(2)**2+G(3)**2+H(3)**2)   ! DIPMOM in (nT RE^3)
  xme=DIPMOM*re_m**3*1.e-9               ! dipole moment in (T m^3)
  write(*,*) 'xme ',xme

! MATE exosphere model
  if (igeo==3) call read_mate_data
  TimeMate=TimeMate+float(iday-iDoyStartMate)*86400.

! Find location of dipole north in geographic coordinates: elon,ctp,stp
  xmag=0.
  ymag=0.
  zmag=rc
  call geomag_08(xgeo,ygeo,zgeo,xmag,ymag,zmag,m_one)
  ctp=zgeo/rc                  ! cosine of angle between geographic N & dipole N
  stp=sqrt(1.-ctp*ctp)         ! sine of angle between geographic N and dipole N
  elon=atan2(ygeo,xgeo)*180./pi  ! geographic longitude of diple N in degree 
  write(*,*) 'elon ',elon
  write(*,*) 'ctp,stp ',ctp,stp 

! Determine whether output from Yu Lin's hybrid code is needed
  ihy=0
  if (imod.eq.3.or.iplsh.eq.3.or.ihigh.eq.3) ihy=1

! FLC or no FLC
  if (iflc.eq.1) write(*,*) ' FLC scattering is used'
  if (iflc.eq.0) write(*,*) ' no FLC scattering'

  end subroutine readInputData


!***********************************************************************
!                        grids
!                  set up all the grids
!***********************************************************************
  subroutine grids(ijs,js)

  use constants
  use cread1
  use cread2, only: iyear,intB,rb
  use cgrid
  use crtp, only: nth,Rrtp,Trtp,Prtp
  implicit none
  integer js(ns),k2,i,k,m,j,n,ijs
  integer,parameter :: nlg=300    ! no. of pt along long in vol & corPot calc.
  integer,parameter :: npg=97     ! no. of azimuthal pt in corPot calc.
  real volg(nlg),xlatg(nlg,ip)
  real thetag(nlg),xlongg(npg),potentg(nlg,npg),cor,sint2,potentc1
  real y_data(ig),lngridp1,lnk,lnk1,varL1,varL2,varLi,rc_m,delR,th1,th2,dth, &
       dlath,gride1,gride2,Eo,gridp1,gridp2,c2m2,c1m1,r_mass,d3, &
       xlatmn,xlatmx,xlatdd,pcsq,coslat,sinlat,vfactor,volgi, &
       xgeo,ygeo,zgeo,xi,yi,zi,theta,xlong,Bri1,Btheta,Bphi,xlatu(0:ir+1,ip)
  character header*80

  data y_data/0.009417,0.019070,0.037105,0.069562,0.122536,0.199229, &
              0.296114,0.405087,0.521204,0.638785,0.750495,0.843570, &
              0.910858,0.952661,0.975754,0.988485,0.995792,0.998703/

  cor=2.*pi/86400.   ! corotation speed in rad/s
  rc_m=rc*re_m       ! ionosphere distance in meter
  xlatmn=18.         ! for dipole internal field (intB=0)
  xlatmx=70.5        !

! Set up geomagnetic longitude grid at the ionosphere.
  dphi=2.*pi/ip    
  do j=1,ip
     mphi(j)=(j-1)*dphi      ! Euler potential beta (mag longitude in radian)
     mlon(j)=mphi(j)*180./pi    ! mlong in degree 
     if (mlon(j).gt.180.) mlon(j)=mlon(j)-360.  ! mlon ranging -180 to 180 deg
  enddo

! Find flux tube volume per unit B flux and corresponding xlatg if intB>0
  if (intB.gt.0) call make_igrfvol(nlg,mlon,volg,xlatg)

! Find corotation potential at a fixed geographic grid if intB>0
  if (intB.gt.0) call corPotentgeo(cor,nlg,npg,thetag,xlongg,potentg)

! Setup xlati, dlati, varL, ksai, potentc(corotation potential)
  if (intB.eq.0) then
     varL1=rc/cosd(xlatmn)/cosd(xlatmn)
     varL2=rc/cosd(xlatmx)/cosd(xlatmx)
  else
     vfactor=1.3e-3
     varL1=vfactor*volg(1)**0.25
     varL2=vfactor*volg(nlg)**0.25
  endif
  dvarL=(varL2-varL1)/(ir*1.+1.5)
  do i=0,ir+1
     ! setup varL and xlatd
     varLi=varL1+i*dvarL
     if (i.ge.1) varL(i)=varLi
     if (i.ge.1.and.i.le.ir) write(*,*) 'i,varL ',i,varL(i)
     if (intB.eq.0) xlatd(i,:)=acos(sqrt(rc/varLi))
     if (intB.gt.0) then
        volgi=(varLi/vfactor)**4
        do j=1,ip
           call lintp(volg,xlatg(:,j),nlg,volgi,xlatdd)
           xlatd(i,j)=xlatdd*pi/180.
        enddo
     endif
     ! setup xlatu, the upper grid boundary of xlatd
     varLi=varLi+0.5*dvarL
     if (intB.eq.0) xlatu(i,:)=acos(sqrt(rc/varLi))
     if (intB.gt.0) then
        volgi=(varLi/vfactor)**4
        do j=1,ip
           call lintp(volg,xlatg(:,j),nlg,volgi,xlatdd)
           xlatu(i,j)=xlatdd*pi/180.
        enddo
     endif
  enddo
  do j=1,ip
     do i=1,ir+1
        if (i.le.ir) xlati(i,j)=xlatd(i,j)
        dlati(i,j)=xlatu(i,j)-xlatu(i-1,j)
        if (j.eq.1.and.i.le.ir) write(*,*) 'i,xlati,dlati ', &
                                i,xlati(i,j)*180./pi,dlati(i,j)*180./pi
     enddo
     write(*,'(a,i6,2f14.5)') 'j,xlati(1),xlati(ir) ', &
                               j,xlati(1,j)*180./pi,xlati(ir,j)*180./pi
     do i=1,ir+1  
        coslat=cos(xlatd(i,j))
        sinlat=sin(xlatd(i,j))
        xi=rc*coslat*cos(mphi(j))
        yi=rc*coslat*sin(mphi(j))
        zi=rc*sinlat
        call GEOMAG_08(xgeo,ygeo,zgeo,xi,yi,zi,m_one)    ! mag to geo
        theta=acos(zgeo/rc)
        sint2=(sin(theta))**2
        if (intB.eq.0) then
           Bri1=2.*xme*sinlat/rc_m**3
           if (i.le.ir) potentc(i,j)=-xme*cor*sint2/rc_m
        else
           xlong=atan2(ygeo,xgeo)
           call IGRF_GEO_08(rc,theta,xlong,Bri1,Btheta,Bphi)
           Bri1=abs(Bri1)*1.e-9           ! Br in Tesla
           call lintp2(thetag,xlongg,potentg,nlg,npg,theta,xlong,potentc1)
           if (i.le.ir) potentc(i,j)=potentc1
        endif
        ksai(i,j)=Bri1*rc_m*rc_m*coslat*dlati(i,j)/dvarL
     enddo
  enddo

! Setup xlath, latitude grid beyond xlati 
  xlatdd=maxval(xlatd(ir+1,:))*180./pi
  dlath=(90.-xlatdd)/(irh-1)
  do i=1,irh
     xlath(i)=xlatdd+(i-1)*dlath
  enddo

! mass no. of H+, O+, He, N+, and e-
  xmass1(1:5)=[1.,16.,4.,14.,5.4462e-4]

! Set up gride, gridp, dlnp, ekev, vel, gride
  do n=1,ijs
     gride1=0.1 
     gride2=500.
     if (js(n).ge.5) gride1=gride1*10.
     if (js(n).ge.5) gride2=gride2*10.
     ! find dlnp
     xmass(n)=xmp*xmass1(js(n))    ! rest mass of each species (kg)
     Eo=511.*xmass(n)/e_mass       ! rest E in keV
     gridp1=sqrt(gride1*(gride1+2.*Eo))*1.6e-16/EM_speed
     gridp2=sqrt(gride2*(gride2+2.*Eo))*1.6e-16/EM_speed
     lngridp1=log(gridp1)
     dlnp(n)=(log(gridp2)-lngridp1)/iw
     ! find lnp, gridp, ekev, vel, pcEo
     do k=0,iw+1
        lnp(n,k)=lngridp1+(k*1.-0.5)*dlnp(n)
        gridp(n,k)=exp(lnp(n,k))
        pcsq=(gridp(n,k)*EM_speed/1.6e-16)**2
        c2m2=pcsq+Eo*Eo
        c1m1=sqrt(c2m2)
        ekev(n,k)=c1m1-Eo
        r_mass=Eo/sqrt(c2m2)                              ! mo/m
        vel(n,k)=gridp(n,k)*r_mass/xmass(n)               ! v in m/s
        if (k.ge.1.and.k.le.iw) pcEo(n,k)=pcsq/2./c1m1    ! in keV
     enddo
     ! find gride and ebound
     ebound(n,0)=sqrt(ekev(n,0)*ekev(n,1))
     do k=1,je   
        k2=k*iw/je
        gride(n,k)=sqrt(ekev(n,k2)*ekev(n,k2-1))
        ebound(n,k)=sqrt(ekev(n,k2)*ekev(n,k2+1))
     enddo
     if (js(n).eq.1) gride_i(:)=gride(n,:)
     if (js(n).eq.5) gride_e(:)=gride(n,:)
     ! find dkeV
     do k=1,iw
        k2=k*je/iw
        if (mod(k,2).eq.0) dkeV(n,k)=ebound(n,k2)-gride(n,k2)
        if (mod(k,2).ne.0) dkeV(n,k)=gride(n,k2+1)-ebound(n,k2)
     enddo
  enddo

! Grids in K (T^0.5 m)
      lnk1=log(58.)
      dlnk=0.267
      do m=0,ik+1
         lnk=lnk1+(m-1)*dlnk  
         xk(m)=exp(lnk) 
      enddo

! find K at grid boundary
  do m=1,ik
     xkb(m)=sqrt(xk(m)*xk(m+1))
  enddo

! Set up y grid of output f and flux
  gridy(1:ig)=y_data(1:ig)
  do m=1,ig
     sinSq(m)=gridy(m)*gridy(m)
     cosSq(m)=1.-sinSq(m)
  enddo

! Setup d4 
  d3=dvarL*dphi*dlnk
  do n=1,ijs
     d4(n)=d3*dlnp(n)
  enddo 

! Construct the spherical grid, rtp
  delR=(rb-rc)/(ir-1)
  do i=1,ir
     Rrtp(i)=rc+(i-1)*delR
  enddo
  th1=pi/2.-maxval(xlati(ir,:))
  th2=pi-th1
  dth=(th2-th1)/(nth-1)
  do m=1,nth
     Trtp(m)=th1+(m-1)*dth
  enddo
  Prtp(1:ip)=mphi(1:ip)
  Prtp(ip+1)=Prtp(1)+2.*pi

  end subroutine grids


!*****************************************************************************
!                           readDiffCoef
!  Routine reads chorus and hiss diffusion coeff.
!*****************************************************************************
  subroutine readDiffCoef(ichor,ihiss)

  use constants
  use cread1
  use cread2, only: ijs,js,iEMICdiff
  use cgrid, only: ekeV
  use waveDiffCoef
  use cWpower
  implicit none

 integer i,j,k,m,n,ichor,ihiss,iband
 real pa,daa,dap,dpp,EEo,daadpp,DD0,ekeVlog(0:iw+1),pkeVlog(0:iw+1)
 real,allocatable,dimension(:,:) :: wDpp,wDaa,wDap
 real,allocatable,dimension(:) :: Daa1,Dap1,Dpp1,ckeV,ckeVlog,hkeV,hkeVlog
 character header*80,fhead*33
 ! variables for EMIC
 integer :: iSpec,ipEmic1,iwe,ipa1,iHband=0,iHeband=0,iReadPA=0
 real  :: LEmic1
 real,allocatable,dimension(:) :: EmicKeV,EmicKeVlog,EmicOmpe1,ePA1
 logical :: DoReadEmic=.false.
 logical :: DoAllocateEmic=.false.

 fhead='D_LBchorus_BUMLAT'
! Setup ekeVlog
  do n=1,ijs
     if (js(n).eq.5) ekeVlog(:)=log10(ekev(n,:))     ! electron energy grid
  enddo

! Read LB chorus ckeV, Daa, Dap, Dpp at L = 6.5 and map to cimi grid
  if (ichor.gt.0) then
     do n=1,kLat
        if (n.eq.1) open(unit=48,file=trim(fhead)//'010.dat',status='old')
        if (n.eq.2) open(unit=48,file=trim(fhead)//'1020.dat',status='old')
        if (n.eq.3) open(unit=48,file=trim(fhead)//'2030.dat',status='old')
        read(48,*) Cpower0
        read(48,*) iLc,ipc,iwc,ipa
        if (n.eq.1) then     ! allocate arrays
           allocate (cPA(ipa),cOmpe(ipc),ckeV(iwc),ckeVlog(iwc),xLc(iLc), &
                     BLc(iLc))
           allocate (wDaa(iwc,ipa),wDap(iwc,ipa),wDpp(iwc,ipa))
           allocate (Daa1(iwc),Dap1(iwc),Dpp1(iwc))
           allocate (cDaa(ipc,0:iw+1,ipa,iLc,kLat), &
                    cDap(ipc,0:iw+1,ipa,iLc,kLat),cDpp(ipc,0:iw+1,ipa,iLc,kLat))
        endif
        read(48,'(a80)') header
        read(48,*) xLc  
        read(48,'(a80)') header
        read(48,*) cOmpe
        read(48,'(a80)') header
        read(48,*) ckeV
        read(48,'(a80)') header
        read(48,*) cPA 
        ckeVlog(:)=log10(ckeV(:))
        do i=1,iLc
        BLc(i)=xme/(xLc(i)*re_m)**3      ! dipole B at xLc(i)
        do j=1,ipc
           ! Read coefficients
           do k=1,iwc
              read(48,'(a80)') header
              read(48,'(a80)') header
              do m=1,ipa
                 read(48,*) pa,daa,dap,dpp
                 wDaa(k,m)=daa    ! coeff in 1/sec
                 wDap(k,m)=dap    !
                 wDpp(k,m)=dpp    !
              enddo
           enddo
           ! map coefficients to ekev grid
           do m=1,ipa
              Daa1(:)=wDaa(:,m)
              Dap1(:)=wDap(:,m)
              Dpp1(:)=wDpp(:,m)
              do k=0,iw+1
                 cDaa(j,k,m,i,n)=0.
                 cDap(j,k,m,i,n)=0.
                 cDpp(j,k,m,i,n)=0.
                 if(ekevlog(k).ge.ckeVlog(1).and.ekevlog(k).le.ckeVlog(iwc))then
                    call lintp(ckevlog,Daa1,iwc,ekevlog(k),DD0)
                    cDaa(j,k,m,i,n)=DD0
                    call lintp(ckevlog,Dpp1,iwc,ekevlog(k),DD0)
                    cDpp(j,k,m,i,n)=DD0
                    DaaDpp=sqrt(cDaa(j,k,m,i,n)*cDpp(j,k,m,i,n))
                    call lintp(ckevlog,Dap1,iwc,ekevlog(k),DD0)
                    if (DaaDpp.lt.abs(DD0)) DD0=DD0*DaaDpp/abs(DD0)
                    cDap(j,k,m,i,n)=DD0
                 endif
              enddo
           enddo
        enddo          ! end of do j=1,ipc
        enddo          ! end of do i=1,iLc
        close(48)
     enddo             ! end of do n=1,kLat
     deallocate (wDaa,wDap,wDpp,Daa1,Dap1,Dpp1)
  endif     ! if (ichor.gt.0)

! Read hiss Daa, Dap, Dpp at L = 5.5
  if (ihiss.gt.0) then
     open(unit=48,file='D_hiss_BU.dat',status='old')
     read(48,*) Hpower0
     read(48,*) iph,iwh,jpa,iLh
     allocate (hPA(jpa),hOmpe(iph),hkeV(iwh),hkeVlog(iwh),xLh(iLh),BLh(iLh))
     allocate (wDaa(iwh,jpa),wDap(iwh,jpa),wDpp(iwh,jpa))
     allocate (Daa1(iwh),Dap1(iwh),Dpp1(iwh))
     allocate (hDaa(iph,0:iw+1,jpa,iLh),hDap(iph,0:iw+1,jpa,iLh), &
               hDpp(iph,0:iw+1,jpa,iLh))
     hDaa(:,:,:,:)=0.
     hDap(:,:,:,:)=0.
     hDpp(:,:,:,:)=0.
     read(48,'(a80)') header
     read(48,*) hOmpe
     read(48,'(a80)') header
     read(48,*) hkeV
     read(48,'(a80)') header
     read(48,*) hPA 
     read(48,'(a80)') header
     read(48,*) xLh 
     hkeVlog(:)=log10(hkeV(:))
     do i=1,iLh
        BLh(i)=xme/(xLh(i)*re_m)**3      ! dipole B at xLh(i)
        do j=1,iph
           read(48,'(a80)') header
           read(48,'(a80)') header
           ! read coefficients
           do m=1,jpa
              do k=1,iwh
                 read(48,*) pa,EEo,Daa,Dap,Dpp
                 wDaa(k,m)=Daa
                 wDap(k,m)=Dap
                 wDpp(k,m)=Dpp
              enddo
           enddo
           ! map coefficients to ekev grid
           do m=1,jpa
              Daa1(:)=wDaa(:,m)
              Dap1(:)=wDap(:,m)
              Dpp1(:)=wDpp(:,m)
              do k=0,iw+1
                 hDaa(j,k,m,i)=0.
                 hDap(j,k,m,i)=0.
                 hDpp(j,k,m,i)=0.
                 if(ekevlog(k).ge.hkeVlog(1).and.ekevlog(k).le.hkeVlog(iwh))then
                    call lintp(hkevlog,Daa1,iwh,ekevlog(k),DD0)
                    hDaa(j,k,m,i)=DD0
                    call lintp(hkevlog,Dpp1,iwh,ekevlog(k),DD0)
                    hDpp(j,k,m,i)=DD0
                    DaaDpp=sqrt(hDaa(j,k,m,i)*hDpp(j,k,m,i))
                    call lintp(hkevlog,Dap1,iwh,ekevlog(k),DD0)
                    if (DaaDpp.lt.abs(DD0)) DD0=DD0*DaaDpp/abs(DD0)
                    hDap(j,k,m,i)=DD0
                 endif
              enddo
           enddo      ! end of mapping coefficients to ekev grid
        enddo         ! end of do j=1,iph
     enddo            ! end of do i=1,iLh
     deallocate (wDaa,wDap,wDpp,Daa1,Dap1,Dpp1)
     close(48)
  endif        ! if (ihiss.gt.0)

! Read H and He band EMIC for proton Daa, Dap, Dpp at L = 4.
  write(*,*) ' ... reading He band EMIC waves Daa, Dap, Dpp'
  ipEmic=0
  EmicPower0=1.
  allocate(BEmic(3)) ! H,He,and O band 
  BEmic(1)=xme/(re_m*LEmic1)**3
  BEmic(2)=xme/(re_m*LEmic1)**3
  if (iHeband==0) then
     iHeband=1
  endif
  if (iHband==0) then
     iHband=1
  endif
  if (iEMICdiff>0) then
     do iband=1,2   ! 1=H-band, 2=He-band
        do n=1,ijs 
           if (iband==1.and.js(n)==1) then
              open(unit=48,file='D_EMIC_H_p.dat',status='old')
              iSpec=1 ! H+ ring current
              DoReadEmic=.true.
           else if (iband==1.and.js(n)==2) then
              open(unit=48,file='D_EMIC_H_o.dat',status='old')
              iSpec=1 ! H+ ring current
              DoReadEmic=.true.
           else if (iband==2.and.js(n)==1) then
              open(unit=48,file='D_EMIC_He_p.dat',status='old')
              iSpec=2 ! H+ ring current
              DoReadEmic=.true.
           else if (iband==2.and.js(n)==2) then
              open(unit=48,file='D_EMIC_He_o.dat',status='old')
              iSpec=2 ! H+ ring current
              DoReadEmic=.true.
           endif 
     !if (iEMICdiff>0.and.any(js(1:ijs)==1)) then
     !   open(unit=48,file='D_EMIC_He_p.dat',status='old')
     !   iSpec=1 ! H+ ring current
           if (DoReadEmic) then
              read(48,*) header
              read(48,*) LEmic1
              read(48,*) ipEmic1,iwe,ipa1
              ! allocate Emic Daa, Dap, Dpp
              if (.not.allocated(EmicHeDaa)) &
                 allocate (EmicHeDaa(ijs,ipEmic1,0:iw+1,ipa1),&
                           EmicHeDap(ijs,ipEmic1,0:iw+1,ipa1), &
                           EmicHeDpp(ijs,ipEmic1,0:iw+1,ipa1))
              if (.not.allocated(EmicHDaa)) &
                 allocate (EmicHDaa(ijs,ipEmic1,0:iw+1,ipa1),&
                           EmicHDap(ijs,ipEmic1,0:iw+1,ipa1), &
                           EmicHDpp(ijs,ipEmic1,0:iw+1,ipa1))
              allocate (ePA1(ipa1),EmicOmpe1(ipEmic1),EmicKeV(iwe),EmicKeVlog(iwe))
              allocate (wDaa(iwe,ipa1),wDap(iwe,ipa1),wDpp(iwe,ipa1))
              allocate (Daa1(iwe),Dap1(iwe),Dpp1(iwe))
              read(48,*) EmicOmpe1
              read(48,*) EmicKeV
              EmickeVlog(:)=log10(EmicKeV(:))
              do j=1,ipEmic1
                 do k=1,iwe
                    read(48,*) header
                    read(48,*) header
                    do m=1,ipa1
                       read(48,*) pa,Daa,Dap,Dpp
                       wDaa(k,m)=Daa
                       wDap(k,m)=Dap
                       wDpp(k,m)=Dpp
                       if (iReadPA==0) ePA1(m)=pa
                    enddo ! end of m
                    iReadPA=1
                 enddo    ! end of k, iwe
                 ! map coefficients to ekev grid
                 do m=1,ipa1
                    Daa1(:)=wDaa(:,m)
                    Dap1(:)=wDap(:,m)
                    Dpp1(:)=wDpp(:,m)
                    do k=0,iw+1
                       if      (iband==1) then
                          EmicHDaa(iSpec,j,k,m)=0.
                          EmicHDap(iSpec,j,k,m)=0.
                          EmicHDpp(iSpec,j,k,m)=0.
                       else if (iband==2) then
                          EmicHeDaa(iSpec,j,k,m)=0.
                          EmicHeDap(iSpec,j,k,m)=0.
                          EmicHeDpp(iSpec,j,k,m)=0.
                       endif
                       if(pkevlog(k).ge.EmicKeVlog(1).and.&
                          pkevlog(k).le.EmickeVlog(iwe))then
                          call lintp(EmicKevlog,Daa1,iwe,pkevlog(k),DD0)
                          !EmicHeDaa(iSpec,j,k,m)=DD0
                          if (iband==1) EmicHDaa(n,j,k,m)=DD0
                          if (iband==2) EmicHeDaa(n,j,k,m)=DD0
                          call lintp(EmicKevlog,Dpp1,iwe,pkevlog(k),DD0)
                          !EmicHeDpp(iSpec,j,k,m)=DD0
                          if (iband==1) EmicHDpp(n,j,k,m)=DD0
                          if (iband==2) EmicHeDpp(n,j,k,m)=DD0
                          if (iband==1) DaaDpp=sqrt(EmicHDaa(iSpec,j,k,m)*EmicHDpp(iSpec,j,k,m))
                          if (iband==2) DaaDpp=sqrt(EmicHeDaa(iSpec,j,k,m)*EmicHeDpp(iSpec,j,k,m))
                          call lintp(EmicKevlog,Dap1,iwe,pkevlog(k),DD0)
                          if (DaaDpp.lt.abs(DD0)) DD0=DD0*DaaDpp/abs(DD0)
                          !EmicHeDap(iSpec,j,k,m)=DD0
                          if (iband==1) EmicHDap(n,j,k,m)=DD0
                          if (iband==2) EmicHeDap(n,j,k,m)=DD0
                       endif
                    enddo ! end if k,iw
                 enddo    ! end of mapping coefficients to ekev grid
              enddo       ! end of j
              if (ipEmic==0.and.iband==1) then
                 ipEmic=ipEmic1
                 iePA=ipa1
                 allocate(EmicOmpe(ipEmic))
                 allocate(ePA(iepa))
                 EmicOmpe(:)=EmicOmpe1(:)
                 ePA(:)=ePA1(:)
              else
                 if (ipEmic/=ipEmic1) write(*,*) ' **WARNING: ipEmic /= ipEmic1'      
                 if (iepa/=ipa1) write(*,*) ' **WARNING: iepa /= ipa1'      
              endif       
              close(48)
              deallocate (EmicOmpe1,ePA1,EmicKeV,EmicKeVlog)
              deallocate (wDaa,wDap,wDpp,Daa1,Dap1,Dpp1)
           endif
        enddo
     enddo
  endif

  end subroutine readDiffCoef


!**************************************************************************
  subroutine qt_conductance
!**************************************************************************
! Routine reads quiet time or background conductance from storm.rcmcon1a

  use constants
  use cread1
  use cread2,only: icon
  use qtconductance
  implicit none
  integer iday,j
  real f107r,Ap,colat1,colatidim

  if (icon.eq.1) open(unit=15,file=trim(storm)//'.rcmcon1a',status='old')
  if (icon.eq.2) open(unit=15,file='CIMI.rcmcon1a',status='old') ! default
  read(15,*) iday,f107r,Ap,idim1,jdim1,colat1,colatidim
  allocate (qtplam0(idim1,jdim1),qthall0(idim1,jdim1),qtped0(idim1,jdim1))
  allocate (ss0(jdim1),colat0(idim1),aloct0(jdim1))
  read(15,*) qtplam0,qthall0,qtped0,ss0
  read(15,*) colat0,aloct0  ! in radian, zero aloct0 --> noon
  close(15)
  aloct0(1)=aloct0(1)-2.*pi
  aloct0(2)=aloct0(2)-2.*pi
  aloct0(jdim1)=aloct0(jdim1)+2.*pi

  write(*,*) ' '
  write(*,*) 'in *.rcmcon1a, iday = ',iday
  write(*,*) '             F10.7 = ',f107r,', Ap = ',Ap
  write(*,*) '             idim1 = ',idim1,', jdim1 = ',jdim1
  write(*,*) '             colat(1) = ',colat1
  write(*,*) '             colat(idim) = ',colatidim
  write(*,*) '  '

  end subroutine qt_conductance


!**************************************************************************
!                            RCMsetup
!  Routine setup the rcm potential calculation mechinery and other things
!**************************************************************************
      subroutine RCMsetup(t,xlati,xmlt)
   
! Input: xlati,xmlt,qtplam0,qthall0,qtped0,ss0,colat0,aloct0
! Output: 'rcmcrd11','rcmcrd21','rcmcon1a'

      use constants
      use rcmgrid
      use cread2
      use cfield,only: Apt,f107,sini_f
      use coefficients
      use qtconductance
      use conductance

      real xlati(ir,ip),xmlt(ip)

!  Create files rcmcrd11 and rcmcrd21
      call fok2rcm_grd(xlati,xmlt,hlosscone,sini_f)

!  Make rcm grid.
      open(unit=10,file='rcmcrd11',status='old')
         read(10,*) id,a2,b2,aval,offset,dlam,dpsi,ri,re
      close(10)
      imin=1
      iint=1
      jint=1
      open(unit=19,file='rcmcrd21',status='old')
         read(19,*) idummy,jdummy       ! MCF added on 19 June 2003 
         read(19,*) colat
         read(19,*) aloct
         read(19,*) sini
         read(19,*) alpha
         read(19,*) beta
         read(19,*) vcorot
         read(19,*) bir
      close(19)
      colat1=0.5*pi-maxval(xlati)
      colatidim=0.5*pi-minval(xlati)
      if (t.eq.tstart) then
         print *,' ' 
         print *,'in cimi, iday = ',iday
         print *,'         F10.7 = ',f107,', Ap = ',Apt   
         print *,'         idim = ',idim,', jdim = ',jdim
         print '(a,f9.3)','         colat(1) = ',colat1
         print '(a,f9.3)','         colat(idim) = ',colatidim
      endif

! Calculate background conductances at (colat,aloct) and write to rcmcon1a
  do j=1,jdim
     aloct1=aloct(idim,j)       ! at equatorial boundary
     call lintp(aloct0,ss0,jdim1,aloct1,ss(j))
     do i=1,idim
        colat1=colat(i,j)
        aloct1=aloct(i,j)
        call lintp2(colat0,aloct0,qtplam0,idim1,jdim1,colat1,aloct1,qtplam(i,j))
        call lintp2(colat0,aloct0,qthall0,idim1,jdim1,colat1,aloct1,qthall(i,j))
        call lintp2(colat0,aloct0,qtped0,idim1,jdim1,colat1,aloct1,qtped(i,j))
     enddo
  enddo

! Add 2mho to qtplam, qthall, qtped
  qtplam(:,:)=qtplam(:,:)+8.
  qthall(:,:)=qthall(:,:)+8.
  qtped(:,:)=qtped(:,:)+8.

! Write qt conductances
  open(unit=16,file='rcmcon1a')
  write(16,*) qtplam,qthall,qtped,ss
  close(16)

  end subroutine RCMsetup


!***********************************************************************
!                             fieldpara
! Routine calculates kinetic energy, velocity, y, latitude and altitude
! at mirror point, etc, for given magnetic moment, K and position for a
! given magnetic field configuration.             
!**************************************************************************
  subroutine fieldpara(t)

  use constants
  use cgrid
  use cread1
  use cread2
  use cfield
  use cinitial,only : f2
  use rcmgrid,only : jdim,jdif
  use ModCurvScatt, only: iflc,tflc,calc_flc_para,calc_Daa_flc,Daa_flc
  use ModMate, only: mate_geocorona 
  parameter (np=1000)
  integer indx(np),iba0(ip),mstop(4)
  real xk3(np),bm1(np),rm(np),xa(np),ya(np),za(np),dss(np),dssi(np),yint(np),&
       yint1(np),yinth(np),h3(np),bba(np),xa2(np),ya2(np),za2(np),zabs(np), &
       rs(np),xs(np),ys(np),zs(np),bs(np),tya3(np),xa1(np),ya1(np),za1(np), &
       bm0(ir,ip,0:ik+1),cosPit(0:ik+1),MLT1D(nMLTg+1),dmltn(ip), &
       cww(0:iw),fww(0:iw+1),fw0(iw),faw(0:iw),ro1(ir)
  character header*80
  COMMON /GEOPACK1/ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI, &
          SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33, &
          E11,E21,E31,E12,E22,E32,E13,E23,E33

!  Parameters for fieldline tracing
   dir=1.
   dir1=-1.
   rlim=1.2*rb
   rmin=0.07
   err=0.0001
   dsmax=0.2          

!  Save previous bm and iba
      if (t.gt.tstart) then
         bm0(:,:,:)=bm(:,:,:)
         iba0(:)=iba(:)
      endif

!  Setup parmod0 at tstart 
   thalf=t+0.5*tstep
   if (t.eq.tstart) then
      if (ires.eq.0.and.itype.eq.2) then
         open(unit=11,file=trim(outname)//'.le',status='old')
         read(11,'(a80)') header
         read(11,*) parmod0
         close(11)
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

! Divide the inner magnetosphere into 3 bands
  xL1=3.      ! ro = 1 - xL1, inner belt and slot region
  xL2=6.0     ! ro = xL1 -xL2, outer belt; L > xL2, outer belt to NEPB

! Setup 1D arrays of MLTgeo if iplsh=4
  if (iplsh.eq.4) then
     call locate1(tgeo,ntg,thalf,ntg1)
     if (ntg1.eq.0) ntg1=1
     if (ntg1.lt.ntg) then
        if ((tgeo(ntg1+1)-t).lt.(t-tgeo(ntg1))) ntg1=ntg1+1
     endif
     MLT1D(1:nMLTg)=MLTgeo(ntg1,1:nMLTg)
     MLT1D(nMLTg+1)=MLT1D(1)+24.
  endif
 
! Find phi, xmlt, and rsb and xmltb, model boundary on equatorial plane
  dr0=0.05
  dmlt0=2.
  xmltd1=7.0          ! dayside is between xmltd1 and xmltd2
  xmltd2=17.0         ! 
  xlatmp=60.*pi/180.  ! min lat for dayside magnetopause
  dlatr=0.01*pi/180.
  do j=1,ip
     ! find phi from mphi(MLONG in radian)
     sinlat=sin(xlati(ir,j))
     coslat=cos(xlati(ir,j))
     xi=rc*coslat*cos(mphi(j))
     yi=rc*coslat*sin(mphi(j))
     zi=rc*sinlat
     call MAGSM_08(xi,yi,zi,xsm,ysm,zsm,i_one)     ! mag to sm
     phi(j)=atan2(ysm,xsm)+pi       ! MLT in radian, zero is midnight
     ! find rsb and xmltb
     xmlt(j)=phi(j)*12./pi
     phi1=phi(j)-pi              ! zero is noon
     xlat0=50.*pi/180.           ! min xlati for rsb 
     xlat2=xlati(ir,j)
     xlat1=xlat2
     find_rsb: do i=1,300
        xi=rc*cos(xlat1)*cos(phi1)
        yi=rc*cos(xlat1)*sin(phi1)
        zi=rc*sin(xlat1)
        ! trace field line from the northern ionosphere to southern ionosphere
        call traceF(imod,intB,xi,yi,zi,dir,dsmax,err,rlim,rc,parmod,psi, &
                    np,xf,yf,zf,xa1,ya1,za1,rs,bs,npf,iout)
        if (iout.eq.0) then
           zabs(1:npf)=abs(za1(1:npf))
           imn=minloc(zabs(1:npf),dim=1)    ! the point near Zsm=0
           rsb(j)=rs(imn)
           xmltb(j)=atan2(ya1(imn),xa1(imn))*12./pi+12.   ! mlt in hr
           bob(j)=bs(imn)*1.e-9   ! B in Tesla
           if (iplsh.eq.4) then
              xmltb1=xmltb(j)
              if (xmltb1.lt.MLT1D(1)) xmltb1=xmltb1+24.
           endif
           if ((rsb(j)-rb).gt.dr0) iout=2
           if (xi.gt.0..and.abs(xmltb(j)-xmlt(j)).ge.dmlt0) iout=3  ! mlt wrap
           if (iout.eq.0.and.xlat1.gt.xlatmp.and.xmlt(j).gt.xmltd1.and. &
               xmlt(j).lt.xmltd2) then
              ic=0
              find_min: do m=2,npf-1
                 if (bs(m).le.bs(m-1).and.bs(m).le.bs(m+1)) ic=ic+1
                 if (ic.gt.1) then
                    iout=4               ! dayside multiple local B min
                    exit find_min
                 endif
              enddo find_min
           endif
           if (iout.eq.0) then
              call indexx(npf,za1,indx)
              if (imn.le.indx(npf).or.imn.ge.indx(1)) iout=5   ! close to iono
           endif
        endif
        if (iout.eq.0) then
           xlat0=xlat1
           if ((xlat2-xlat0).lt.dlatr) exit find_rsb
           if (abs(rsb(j)-rb).lt.dr0) exit find_rsb
        endif
        if (iout.ge.1) xlat2=xlat1
        xlat1=0.5*(xlat0+xlat2)
     enddo find_rsb
     ! find iba
     call locate1(xlati(:,j),ir,xlat1,iba(j))
     if (xlat1.eq.xlati(ir,j)) iba(j)=ir
  enddo

! Locate noon (jnoon) and jdif (difference in j index with RCM)
  do j=1,ip
     dmltn(j)=abs(12.-xmlt(j))
  enddo
  jnoon=minloc(dmltn(1:ip),dim=1)
  jdif=jdim-jnoon
  if (jdif.gt.ip) jdif=jdif-ip

! Calculate values at north and south footpoints with Bext=0 when t=tstart
  rc3=(rc*re_m)**3
  if (t.eq.tstart) then
     do i=1,ir
        do j=1,ip
           xlat1=xlati(i,j)
           sinlat=sin(xlat1)
           coslat=cos(xlat1)
           BriN(i,j)=2.*xme*sinlat/rc3            ! dipole Bri
           BiN(i,j)=xme*sqrt(1.+3.*sinlat)/rc3    ! dipole Bi
           if (intB.eq.1) then                ! IGRF internal field
              xi=rc*coslat*cos(mphi(j))
              yi=rc*coslat*sin(mphi(j))
              zi=rc*sinlat
              call GEOMAG_08(xgeo,ygeo,zgeo,xi,yi,zi,m_one)    ! mag to geo
              theta=acos(zgeo/rc)
              xlong=atan2(ygeo,xgeo)
              call IGRF_GEO_08(rc,theta,xlong,Bri1,Btheta,Bphi)
              BriN(i,j)=abs(Bri1)*1.e-9         ! Br in Tesla
              BiN(i,j)=sqrt(Bri1*Bri1+Btheta*Btheta+Bphi*Bphi)*1.e-9  ! in Tesla
           endif
           sini_f(i,j)=BriN(i,j)/BiN(i,j)    ! sine of the dip angle
        enddo
     enddo
  endif

! Start field line tracing.  Field line tracing from north to south.
  Bfactor=1.
  jloop: do j=1,ip
     iloop: do i=1,iba(j)
        xlat1=xlati(i,j)
 555    phi1=phi(j)+pi                       ! +x corresponing to noon 
        xi=rc*cos(xlat1)*cos(phi1)
        yi=rc*cos(xlat1)*sin(phi1)
        zi=rc*sin(xlat1)

        ! Add n5 points below rc at north assuming dipole 
        call traceF(i_zero,i_zero,xi,yi,zi,dir1,dsmax,err,rb,rmin,parmod,psi, &
                    np,xf,yf,zf,xa1,ya1,za1,rs,bs,npf,iout)
        if (intB.gt.0) then   ! non-dipole internal field
           call SMGSW_08(xi,yi,zi,xgsm,ygsm,zgsm,i_one)           ! sm to gsm
           call IGRF_GSW_08(xgsm,ygsm,zgsm,bxint,byint,bzint)
           Bfactor=sqrt(bxint*bxint+byint*byint+bzint*bzint)/bs(1)
        endif
        n5=npf-1
        do m=1,n5
           m1=npf+1-m 
           bm1(m)=Bfactor*bs(m1)*1.e-9   ! B in Tesla, scaled with Bi
           xa(m)=xa1(m1)
           ya(m)=ya1(m1)
           za(m)=za1(m1)
           rm(m)=rs(m1)
        enddo 

        ! trace field line from the northern ionosphere to southern ionosphere
        call traceF(imod,intB,xi,yi,zi,dir,dsmax,err,rlim,rc,parmod,psi, &
                    np,xf,yf,zf,xa1,ya1,za1,rs,bs,npf,iout)
        if (iout.eq.0) &
           call EquatorialXing(i,j,xa1,ya1,za1,rs,bs,npf,imn,im0,dr0,iout)
        if (iout.eq.1) then
           if ((rsb(j)-ro(i-1,j)).gt.3.) then
              xlat1=0.9*xlat1+0.1*xlati(i-1,j)
              goto 555        ! redo tracing with MLAT closer to xlati(i-1,j)
           endif
           write(*,*) 'iout=1: j,i,iba ',j,i,iba(j)        
           iba(j)=i-1   
           goto 888      ! do next j
        endif

        ! Check field line equatorial crossing
        if (i.gt.1) then
           if (rs(im0).lt.ro(i-1,j)) then
              write(*,*) 'i,j,ro(i-1,j)>ro(i,j),phi,xlati ', &
                          i,j,ro(i-1,j),rs(im0),phi(j)*180./pi,xlat1*180./pi
              iba(j)=i-1   
              goto 888    ! do next j
           endif 
           if (rs(im0).gt.rsb(j)) then
              rsb(j)=rs(im0)
              xmltb(j)=atan2(ya1(im0),xa1(im0))*12./pi+12.   ! mlt in hr
              bob(j)=bs(im0)*1.e-9   ! B in Tesla
           endif
        endif
        do m=1,npf
           bba(m)=bs(m)*1.e-9   ! B in Tesla
        enddo
        brcN=bba(1)             ! B at rc at North
        ro(i,j)=rs(im0)     
        xmlto(i,j)=atan2(ya1(im0),xa1(im0))*12./pi+12.   ! mlt in hr
        xo(i,j)=xa1(im0)
        yo(i,j)=ya1(im0)

        ! Find southern footpoint 
        brcS=bba(npf)           ! B at rc at South
        BiS(i,j)=brcS
        rf=rs(npf)
        sinlat=zf/rf
        xlatiS(i,j)=asin(sinlat)      ! foot point in sm at southern ionosphere
        phiS(i,j)=atan2(yf,xf)+pi          ! phiS=0 --> midnight
        call MAGSM_08(xmag,ymag,zmag,xf,yf,zf,m_one)     ! sm to mag 
        mlonS(i,j)=atan2(ymag,xmag)*180./pi          
        BriS(i,j)=2.*xme*abs(sinlat)/rc3   ! dipole Bri in Tesla
        if (intB.eq.1) then                ! IGRF internal field
           call SMGSW_08(xf,yf,zf,xgsm,ygsm,zgsm,i_one)           ! sm to gsm
           call GEOGSW_08(xgeo,ygeo,zgeo,xgsm,ygsm,zgsm,m_one)    ! gsm to geo
           theta=acos(zgeo/rc)
           xlong=atan2(ygeo,xgeo)
           call IGRF_GEO_08(rc,theta,xlong,Bri1,Btheta,Bphi)
           BriS(i,j)=abs(Bri1)*1.e-9         ! Br in Tesla
        endif

        ! find the lengths between trace points and the accumulative length
        dssi(1)=0.
        do m=2,npf
           d2=(xa1(m)-xa1(m-1))**2+(ya1(m)-ya1(m-1))**2+(za1(m)-za1(m-1))**2
           dss(m-1)=sqrt(d2)
           dssi(m)=dssi(m-1)+dss(m-1)
        enddo
           
        ! make sure B at rc decreases to bba(imn) and rises
        mstop(1:4)=0
        ims=0
        do mns=-1,1,2
           bgood=bba(imn)
           igood0=1
           if (mns.eq.-1) nsf=2    ! check the north part
           if (mns.eq.1) nsf=npf-1 ! check the south part
           do m=imn+mns,nsf,mns    ! start checking near equator
              igood=-1
              if (bba(m).gt.bgood) igood=1       ! B should be increasing
              if (igood.eq.1) bgood=bba(m)
              if ((igood0*igood).lt.0) then
                ims=ims+1
                mstop(ims)=m
             endif
             igood0=igood
           enddo
        enddo
        ! Do interpolation in B if ims gt 0
        if (ims.gt.0) then
           write(*,'(a,7i8)') 'local Bmin, i,j,imn,mstop: ',i,j,imn,mstop
           do mseg=1,2
              mseg1=mseg*2-1
              if (mstop(mseg1).gt.mstop(mseg1+1)) mns=-1   ! local min at north
              if (mstop(mseg1).lt.mstop(mseg1+1)) mns=1    ! local min at south
              m00=mstop(mseg1)-mns
              m11=mstop(mseg1+1)
              do m=m00+mns,m11-mns,mns
                 dsss=(dssi(m)-dssi(m00))/(dssi(m11)-dssi(m00))
                 bba(m)=bba(m00)+(bba(m11)-bba(m00))*dsss
              enddo
           enddo
        endif
        bo(i,j)=bba(im0)

        ! find the flux tube volume per magnetic flux
        npf1=npf-1   ! number of grid intervals from north to south end
        do m=1,npf1
           b_mid=0.5*(bba(m)+bba(m+1))
           yint(m)=1./b_mid
        enddo
        call closed(i_one,npf1,yint,dss,sso)      ! use closed form
        volume(i,j)=sso*re_m  ! volume per unit flux

        ! Populate arrays from below N ionosphere to S ionosphere
        im2=imn+n5     ! new imn
        do m=1,npf
           ic=m+n5
           bm1(ic)=bba(m)
           xa(ic)=xa1(m)
           ya(ic)=ya1(m)
           za(ic)=za1(m)
           rm(ic)=rs(m)
        enddo

        ! Add n5 points below rc at south 
        call traceF(i_zero,i_zero,xa(ic),ya(ic),za(ic),dir,dsmax,err,rb,rmin, &
                    parmod,psi,np,xf,yf,zf,xa1,ya1,za1,rs,bs,npf,iout)
        if (intB.gt.0) then   ! non-dipole internal field
           call SMGSW_08(xa(ic),ya(ic),za(ic),xgsm,ygsm,zgsm,i_one)    ! sm to gsm
           call IGRF_GSW_08(xgsm,ygsm,zgsm,bxint,byint,bzint)
           Bfactor=sqrt(bxint*bxint+byint*byint+bzint*bzint)/bs(1)
        endif
        n5=npf-1
        do m=ic+1,ic+n5
           m1=m-ic+1 
           bm1(m)=Bfactor*bs(m1)*1.e-9   ! B in Tesla, scaled with Bi
           xa(m)=xa1(m1)
           ya(m)=ya1(m1)
           za(m)=za1(m1)
           rm(m)=rs(m1)
        enddo 

        n=ic+n5-1  ! n = no. of grid intervals from north end to south end
        do m=1,n
           xs(m)=0.5*(xa(m+1)+xa(m))
           ys(m)=0.5*(ya(m+1)+ya(m))
           zs(m)=0.5*(za(m+1)+za(m))
           bs(m)=0.5*(bm1(m+1)+bm1(m))
           d2=(xa(m)-xa(m+1))**2+(ya(m)-ya(m+1))**2+(za(m)-za(m+1))**2
           dss(m)=sqrt(d2)     ! new dss with unwanted points removed
        enddo

!  N ionosphere                       <--dss(m)-->                S ionosphere
!      xs,ys,zs
!           bs(1)                         bs(m)                     bs(n)
!         |-------|-----|......|------|----------|----|.......|---|-------|
!       bm1(1)                      bm1(m)                              bm1(n+1)
!       rm(1)                       rm(m)                               rm(n+1)
!      xa,ya,za

        ! Set up arrarys at trace grids
        xk3(im2)=0.               ! mirroring around Bmin
        !call geocorona(igeo,xa(im2),ya(im2),za(im2),h3(im2))
        call geocorona(igeo,t,xa(im2),ya(im2),za(im2),h3(im2))
        do ngrid=1,im2-1   
           n8=n+1       !  initial value
           line_integral: do ii=ngrid,n
              !call geocorona(igeo,xs(ii),ys(ii),zs(ii),hden)    
              call geocorona(igeo,t,xs(ii),ys(ii),zs(ii),hden)    
              bsi=bs(ii)
              if (bm1(ii+1).ge.bm1(ngrid)) bsi=0.5*(bm1(ngrid)+bm1(ii))
              yint(ii)=sqrt(bm1(ngrid)-bsi)
              bsbm=sqrt(1.-bsi/bm1(ngrid))
              yint1(ii)=1./bsbm
              yinth(ii)=hden/bsbm
              if (bm1(ii+1).ge.bm1(ngrid)) then
                 n8=ii+1
                 exit line_integral
              endif
           enddo line_integral
           n7=n8-1
           n6=n7-1
           dssp=dss(n7)*(bm1(ngrid)-bm1(n7))/(bm1(n8)-bm1(n7)) ! partial ds at S
           call closed(ngrid,n6,yint,dss,sso)  ! use closed form
           sso=sso+yint(n7)*dssp
           xk3(ngrid)=sso*re_m
           call closed(ngrid,n6,yint1,dss,ss1)
           ss1=ss1+yint1(n7)*dssp
           tya3(ngrid)=ss1/ro(i,j)/2.
           call closed(ngrid,n6,yinth,dss,ssh)
           ssh=ssh+yinth(n7)*dssp
           h3(ngrid)=ssh/ss1
        enddo
        tya3(im2)=tya3(im2-1)

        ! find xkcN and xkcS, the K values with mirror points at rc
        call lintp(bm1,xk3,im2,brcN,xkc1)
        call locate1(xkb,ik,xkc1,mc)
        mcN(i,j)=mc
        xkcN(i,j)=xkc1
        call lintp(bm1,xk3,im2,brcS,xkc1)
        call locate1(xkb,ik,xkc1,mc)
        mcS(i,j)=mc
        xkcS(i,j)=xkc1
   
        ! Calculate y, dmu, T(y), bounced average H density
        do m=0,ik+1
           xkm=xk(m)                 ! get Bm @ given K & location
           call lintp(xk3,bm1,im2,xkm,bmmx)
           bm(i,j,m)=bmmx 
           y(i,j,m)=sqrt(bo(i,j)/bmmx)
           if (y(i,j,m).gt.1.) y(i,j,m)=1.
           cosPit(m)=sqrt(1.-y(i,j,m)*y(i,j,m))
           call lintp(xk3,tya3,im2,xkm,tya33)
           tya(i,j,m)=tya33
           call lintp(xk3,h3,im2,xkm,h33)
           if (m.ge.1.and.m.le.ik) then
              lnbm(i,j,m)=log(bmmx)
              Hdens(i,j,m)=h33
           endif
        enddo
        do m=1,ik
           dmu(i,j,m)=0.5*(cosPit(m+1)-cosPit(m-1))
        enddo

        ! field line curvature (FLC) calculation
        if (t.ge.tflc) then
           imn=minloc(bm1(1:n),dim=1)
           call calc_flc_para(i,j,7,4,rm(imn),dss(imn-3:imn+3),xa(imn-3:imn+3),&
                              ya(imn-3:imn+3),za(imn-3:imn+3),bm1(imn-3:imn+3))
        endif

     enddo iloop
888  continue
  enddo jloop

! Make sure iba doesn't increase too much from previous j
  do j=1,ip
     j0=j-1
     if (j0.lt.1) j0=j0+ip
     if ((iba(j)-iba(j0)).gt.4) then
        iba(j)=iba(j0)+4
        rsb(j)=ro(iba(j)+1,j)
        xmltb(j)=xmlto(iba(j)+1,j)
        bob(j)=bo(iba(j)+1,j)
     endif
  enddo

! Make sure iba doesn't increase too much from next j
  do j=ip,1,-1
     j1=j+1
     if (j1.gt.ip) j1=j1-ip
     if ((iba(j)-iba(j1)).gt.4) then
        iba(j)=iba(j1)+4
        rsb(j)=ro(iba(j)+1,j)
        xmltb(j)=xmlto(iba(j)+1,j)
        bob(j)=bo(iba(j)+1,j)
     endif
  enddo

! Fill values beyond iba
  do j=1,ip
     xo_b=-rsb(j)*cos(xmltb(j)*pi/12.)
     yo_b=-rsb(j)*sin(xmltb(j)*pi/12.)
     ib=iba(j)
     ! fill values beyond iba
     do i=ib+1,ir
        ro(i,j)=rsb(j)
        xmlto(i,j)=xmltb(j)
        bo(i,j)=bob(j)
        xo(i,j)=xo_b      
        yo(i,j)=yo_b      
        xlatiS(i,j)=xlatiS(ib,j)
        phiS(i,j)=phiS(ib,j)
        BriS(i,j)=BriS(ib,j)
        BiS(i,j)=BiS(ib,j)
        mlonS(i,j)=mlonS(ib,j)
        volume(i,j)=volume(ib,j)
        y(i,j,1:ik)=y(ib,j,1:ik)
        dmu(i,j,1:ik)=dmu(ib,j,1:ik)
        tya(i,j,1:ik)=tya(ib,j,1:ik)
        bm(i,j,1:ik)=bm(ib,j,1:ik)
        lnbm(i,j,1:ik)=lnbm(ib,j,1:ik)
     enddo
  enddo

! Make sure equatorial crossing and south foot points are smooth in MLT
  do i=1,ir
     imn=minloc(xmlto(i,1:ip),dim=1)
     do j=imn+1,imn+ip-1
        j0=j-1
        if (j0.gt.ip) j0=j0-ip
        j1=j
        if (j1.gt.ip) j1=j1-ip
        j2=j+1
        if (j2.gt.ip) j2=j2-ip
        if (xmlto(i,j1).ge.xmlto(i,j2)) then
           xj1=0.5*(xo(i,j0)+xo(i,j2))
           yj1=0.5*(yo(i,j0)+yo(i,j2))
           xmlto(i,j1)=atan2(yj1,xj1)*12./pi+12.   ! mlt in hr
        endif
     enddo
     imn=minloc(phiS(i,1:ip),dim=1)
     do j=imn+1,imn+ip-2
        j0=j-1
        if (j0.gt.ip) j0=j0-ip
        j1=j
        if (j1.gt.ip) j1=j1-ip
        j2=j+1
        if (j2.gt.ip) j2=j2-ip
        if (phiS(i,j1).ge.phiS(i,j2)) then
           xj1=0.5*(cos(phiS(i,j0))+cos(phiS(i,j2)))
           yj1=0.5*(sin(phiS(i,j0))+sin(phiS(i,j2)))
           phiS(i,j1)=atan2(yj1,xj1)
        endif
     enddo
  enddo

! Calculate xjac and lifetime for loss cone particles
  do j=1,ip
     do n=1,ijs
        ! Calculate xjac
        do m=0,ik+1
           do i=1,ir
              xjac1=4.*pi*ksai(i,j)
              xjac2=xjac1*xk(m)/bm(i,j,m)**1.5
              xjac(n,i,j,0:iw+1,m)=xjac2*gridp(n,0:iw+1)**3
           enddo
           rksai=ksai(ir+1,j)/ksai(ir,j)
           xjac(n,ir+1,j,0:iw+1,m)=xjac(n,ir,j,0:iw+1,m)*rksai
        enddo
        ! Calculate Tbounce and fcone
        do i=1,ir
           mcm=min(mcN(i,j),mcS(i,j))
           ro4=4.*ro(i,j)*re_m
           do m=1,ik
              tcone1=ro4*tya(i,j,m)
              do k=1,iw
                 tcone2=tcone1/vel(n,k)                 ! Tbounce
                 Tbounce(n,i,j,k,m)=tcone2              ! bounce time
                 fcone(n,i,j,k,m)=1.                    ! outside losscone
                 if (m.gt.mcm) fcone(n,i,j,k,m)=1./exp(dt/tcone2) 
              enddo
           enddo
        enddo  
        ! field line curvature (FLC) scattering rate calculation
        if (t.ge.tflc) then
           iflc=2
           do i=1,iba(j)
              call calc_Daa_flc(n,i,j,gridp(n,1:iw),&
                                Tbounce(n,i,j,1:iw,1:ik),y(i,j,1:ik),&
                                bo(i,j),BiN(i,j),BiS(i,j),mcN(i,j),mcS(i,j))
           enddo
           Daa_flc(n,iba(j)+1:ir,j,:,:)=0.
        else
           Daa_flc(n,:,j,:,:)=0.
        endif
     enddo  
  enddo

! Calculate dlnbm/dphi and dlnbm/dvarL at i,j
  dphi2=dphi*2.
  do j=1,ip
     j0=j-1
     if (j0.lt.1) j0=j0+ip
     j1=j+1
     if (j1.gt.ip) j1=j1-ip
     do i=1,ir
        i0=i-1
        if (i0.lt.1) i0=1
        i1=i+1
        if (i1.gt.ir) i1=ir
        dvarL1=varL(i1)-varL(i0)
        do m=1,ik
           dlnBmdL(i,j,m)=log(bm(i1,j,m)/bm(i0,j,m))/dvarL1
           dlnBmdp(i,j,m)=log(bm(i,j1,m)/bm(i,j0,m))/dphi2
        enddo
     enddo
  enddo

! Calculate dlnbm/dphi at i+0.5,j and dlnbm/dvarL at i,j+0.5
  do j=1,ip
     j0=j-1
     if (j0.lt.1) j0=j0+ip
     j1=j+1
     if (j1.gt.ip) j1=j1-ip
     do i=1,ir
        i0=i-1
        if (i0.lt.1) i0=1
        i1=i+1
        if (i1.gt.ir) i1=ir
        dvarL1=varL(i1)-varL(i0)
        do m=1,ik
           Bmi0j5=0.5*(bm(i0,j1,m)+bm(i0,j,m))         ! Bm at i0,j+0.5
           Bmi1j5=0.5*(bm(i1,j1,m)+bm(i1,j,m))         ! Bm at i1,j+0.5
           dlnBmdL1(i,j,m)=log(Bmi1j5/Bmi0j5)/dvarL1   ! dlnBm/dvarL at i,j+0.5
           Bmi5j0=0.5*(bm(i1,j0,m)+bm(i,j0,m))         ! Bm at i+0.5,j0
           Bmi5j1=0.5*(bm(i1,j1,m)+bm(i,j1,m))         ! Bm at i+0.5,j1
           dlnBmdp1(i,j,m)=log(Bmi5j1/Bmi5j0)/dphi2    ! dlnBm/dphi at i+0.5,j
        enddo
     enddo
  enddo
 
! Calculate change of f2 due to dBm/dt
  if (t.gt.tstart) then
     do j=1,ip
        do n=1,ijs
           pp3=(gridp(n,1)/gridp(n,0))**3
           ib1=min(iba(j),iba0(j))
           do i=1,ib1
           do m=1,ik
              ! find cww and nrun
              deltaW=log(bm(i,j,m)/bm0(i,j,m))/2.
              cww1=deltaW/dlnp(n)
              nrun=ifix(abs(cww1)/0.5)+1
              cww(:)=cww1/nrun

              ! advection in W (lnp)
              fww(1:iw)=f2(n,i,j,1:iw,m)
              do nn=1,nrun
                 fw0(1:iw)=fww(1:iw)
                 if (js(n).lt.5) fww(0)=fww(1)         ! f2(0)=f2(1) for ions
                 if (js(n).ge.5) fww(0)=fww(1)/pp3     ! psd(0)=psd(1) for e-
                 fww(iw+1)=fww(iw)    ! f2(iw+1)=f2(iw)
                 call FLS_2O_1D(cww,fww,faw)
                 do k=1,iw
                    fww(k)=fw0(k)+cww(k-1)*faw(k-1)-cww(k)*faw(k)
                 enddo
              enddo

              ! update f2
              f2(n,i,j,1:iw,m)=fww(1:iw)
           enddo
           enddo
        enddo        ! end of j loop
     enddo           ! end of n loop
  endif              ! end of if (t.gt.tstart)

  end subroutine fieldpara


!***********************************************************************
  subroutine EquatorialXing(i,j,xa1,ya1,za1,rs,bs,npf,imn,im0,dr0,iout)
!***********************************************************************
  use cread2,only: rb
  implicit none
  integer,parameter :: mmax=99
  integer npf,im0,m,m1,imx,imn,i,j,ic,iout
  integer m4(mmax),indx(npf),ibmin(4)
  real xa1(npf),ya1(npf),za1(npf),rs(npf),bs(npf),r_xy(npf),rslope0,rslope, &
       rslope2,rslope3,zabs(npf),rbmin(3),dr0,rbdr

  ! Locate multiple crossing of equatorial region or magnetic island
  do m=1,npf
     r_xy(m)=sqrt(xa1(m)*xa1(m)+ya1(m)*ya1(m))
  enddo
  rslope0=(r_xy(2)-r_xy(1))/(za1(2)-za1(1))
  imx=0.
  do m=2,npf-1
     rslope=(r_xy(m+1)-r_xy(m))/(za1(m+1)-za1(m))
     rslope2=rslope0*rslope
     if (rslope2.lt.0..and.za1(m+1).lt.za1(m)) then
        rslope3=abs(rslope0)+abs(rslope)
        if (rslope3.lt.2.) then
           imx=imx+1
           if (imx.gt.mmax) then
              write(*,*) 'Error: imx.gt.mmax, i,j ',i,j
              stop
           endif
           m4(imx)=m                        ! equatorial crossing
        endif
     endif
     rslope0=rslope
  enddo 

! cut short the field line if multiple crossing
  if (imx.gt.1) then    
     ! add a point between m4(1) and m4(imx)+1
     m=m4(1)+1
     m1=m4(imx)+1
     xa1(m)=0.5*(xa1(m-1)+xa1(m1))
     ya1(m)=0.5*(ya1(m-1)+ya1(m1))
     za1(m)=0.5*(za1(m-1)+za1(m1))
     rs(m)=sqrt(xa1(m)*xa1(m)+ya1(m)*ya1(m)+za1(m)*za1(m))
     bs(m)=0.5*(bs(m-1)+bs(m1)) 
     r_xy(m)=sqrt(xa1(m)*xa1(m)+ya1(m)*ya1(m))
     ! renumber the S part
     ic=m
     do m=m1,npf
        ic=ic+1
        xa1(ic)=xa1(m)
        ya1(ic)=ya1(m)
        za1(ic)=za1(m)
        rs(ic)=rs(m)
        bs(ic)=bs(m)
        r_xy(ic)=r_xy(m)
     enddo
     write(*,'(a,4i8,f7.3)')'multi crossing, i,j,npfOld,npfNew ',i,j,npf,ic
     npf=ic               ! new npf
  endif                   ! end of cut short the field line
 
! find local Bmin with largest r_xy
  ic=0
  do m=2,npf-1
     if (bs(m).lt.bs(m-1).and.bs(m).lt.bs(m+1)) then
        ic=ic+1
        if (ic.gt.3) then
           iout=1
           return
        endif
        ibmin(ic)=m        ! local Bmin
        rbmin(ic)=r_xy(m) 
     endif
  enddo
  imn=ibmin(1)
  if (ic.gt.1) then
     call indexx(ic,rbmin,indx)
     imn=ibmin(indx(ic))
  endif

! find crossing magnetic equator
  zabs(1:npf)=abs(za1(1:npf))
  im0=minloc(zabs(1:npf),dim=1)         ! the point near Zsm=0

! check rs(im0)
  rbdr=rb+dr0
  if (rs(im0).gt.rbdr) iout=1

  end subroutine EquatorialXing
 

!***********************************************************************
  subroutine find_ceSigma
!***********************************************************************
! Routine calculates the cross sections of charge exchange

  use ccepara 
  use cgrid, only: ekev
  use cread2, only: ijs,js 
  implicit none
  integer n,k
  real Ekev1,flog_keV,flog_eV,a(10),nene(10),nsig(10),t1,t2,t3,d,x

  data a/1.2e-19,32.0,4.15e-3,2.35e-9,0.792,    &
         2.94e-2,7.26e-3,8.45e-4,8.35e-4,0.563/

  !N+ cross sections from Phaneuf 1978, Table III (in keV, m^2)
  data nene/10.0,20.0,50.0,100.0,150.0, &
            200.0,250.0,300.0,450.0,550.0/

  data nsig/6.2e-20,8.2e-20,7.3e-20,4.6e-20,4.6e-20, &
            3.2e-20,2.8e-20,2.5e-20,1.8e-20,1.4e-20/

  do n=1,ijs
     do k=1,iw
        Ekev1=ekev(n,k)
        flog_keV=log10(Ekev1)
        flog_eV=flog_keV+3.

        ! cross section of ion-H 
        x=flog_keV
        if (x.lt.-2.) x=-2.
        if (js(n).eq.1) &      ! cross section of H+-H
           d=-18.767-0.11017*x-3.8173e-2*x**2-0.1232*x**3-5.0488e-2*x**4
        if (js(n).eq.2) &      ! cross section of O+-H
           d=-18.987-0.10613*x-5.4841e-3*x**2-1.6262e-2*x**3-7.0554e-3*x**4
        if (js(n).eq.3) &      ! cross section of He+-H
           d=-20.789+0.92316*x-0.68017*x**2+0.66153*x**3-0.20998*x**4
        if (js(n).eq.4) &      ! cross section of N+-H
          call lintp(nene,log10(nsig),10,Ekev1,d)
        ceSigma(n,k,1)=10.**d      ! cross section with [H] in m^2

        ! cross section of ion-O
        ceSigma(n,k,2)=0. 
        if (js(n).eq.1) then   ! cross section of H+-O
           t1=exp(-a(2)/Ekev1)/(1.+a(3)*Ekev1**2+a(4)*Ekev1**4.5)
           t2=a(5)*exp(-a(6)*Ekev1)/Ekev1**a(7)
           t3=a(8)*exp(-a(9)*Ekev1)/Ekev1**a(10)
           ceSigma(n,k,2)=a(1)*(t1+t2+t3)    ! in m^2, Perez   
        endif
        if (js(n).eq.2) &      ! cross section of O+-O
           ceSigma(n,k,2)=1.e-20*(32.078-4.911*flog_eV)  ! in m^2, Lindsay
     enddo
  enddo

  end subroutine find_ceSigma


!***********************************************************************
  subroutine cepara
!***********************************************************************
! routine calculates the charge exchange decay rate of ring
! current ion species with geocorona

  use ccepara      
  use cfield,only: Hdens,iba
  use cgrid,only: vel   
  use cread2, only: ijs,js,dt 
  implicit none
  integer n,k,m,i,j
  real vsigmat1,alpha1

  do n=1,ijs
     if (js(n).le.4) then        ! ions only
        do k=1,iw
           vsigmat1=vel(n,k)*ceSigma(n,k,1)*dt
           do j=1,ip
           do i=1,iba(j)
           do m=1,ik
              alpha1=vsigmat1*Hdens(i,j,m)
              achar(n,i,j,k,m)=exp(-alpha1) ! charge exchange decay rate
           enddo
           enddo
           enddo
        enddo
     endif   
  enddo                ! end of do n=1,ijs

  end subroutine cepara


!*******************************************************************************
  subroutine coulpara
!*******************************************************************************
! Routine calcalates the Coulomb drag coefficients (coulii and coulee)
! 1/p*dp/dt=1/pv*dE/dt, where dE/dt is given by Eq 21 in Fok et al, 1993 JGR.

  use constants
  use cCoulpara
  use cread2,only: ijs
  use cgrid
  use cPlasmasphere,only: Tcold,Rcold
  implicit none
  integer n,k,jb
  real e2e,ga1,rfactor,xmasss,gammas,xmm1,vth(4),vvth,coul1,coul2,vup,pup,gfun

  e2e=echarge*echarge/epsilon0
  ga1=e2e*e2e*coullog/4./pi
  coulii(:,:)=0.
  coulee(:,:)=0.

! Calculate plasmasphere thermal velocity
  xmm1=echarge/xmp
  do jb=1,4         ! loop over plasmaspheric species
     vth(jb)=sqrt(2.*Tcold*xmm1/xmass1(jb))
  enddo

  do n=1,ijs
     do k=0,iw
        vup=sqrt(vel(n,k)*vel(n,k+1))      ! vel at upper grid boundary
        pup=sqrt(gridp(n,k)*gridp(n,k+1))  ! p at upper grid boundary
        rfactor=sqrt(1.-(vup/EM_speed)**2)
        xmasss=xmass(n)/rfactor
        gammas=ga1/xmasss/xmasss
        do jb=1,4         ! loop over plasmaspheric species
           xmm1=xmasss/xmass1(jb)/xmp+1.
           vvth=vup/vth(jb)
           coul1=Rcold(jb)*(2.*vvth*vvth*xmm1*gfun(vvth)-erf(vvth))
           if (jb.le.3.and.coul1.gt.0.)coulii(n,k)=coulii(n,k)+coul1
           if (jb.gt.3.and.coul1.gt.0.)coulee(n,k)=coulee(n,k)+coul1
        enddo
        coul2=-xmass(n)*gammas/pup/vup**2      ! at upper grid
        coulii(n,k)=coul2*coulii(n,k)
        coulee(n,k)=coul2*coulee(n,k)
     enddo
  enddo

  end subroutine coulpara


!*******************************************************************************
!                                convection
!  Routine calculates the convection potential
!*******************************************************************************
  subroutine convection(t)
      use constants
      use cread1
      use cread2
      use convect
      use rcmgrid
      use potential
      use conductance
      use plasma
      use coefficients
      use current
      use closs        
      use cfield
      use cgrid
      use cinitial,only: f2
      implicit none

    logical UseAL
    integer kp1,imlt,i,j,n,i0,ic,jc,k,jm,ibmax,j0,m
    real t,SWDen,bx,by,bz,angle,Bt,theta,sina1,colatm1, &
         colatp1,xlat1,dPdL,dPdp,xmlt1,Brik, &
         delmlt,c0,Tilt,fac0,gLat,Ediff1,Ediff2,dvarL2,birk1, &
         gMLT,BnLat,vdrop,Psw,Esw,phim,sigmaHi,phisa,pcp_Hill,phi1,thetab, &
         xlathd,sigmap,sigmah,BoundaryLat,EpotVal,theta1,hallB(idim,jdim), &
         sigmap2I,sigmah2I,pedpsiB(idim,jdim),pedlamB(idim,jdim), &
         Etot,Ptot,Eflux,meanE,EEF,EmeanE,dphi2,Ecri, &
         facm1,facp1,birk1d(idim),Peta(ir,ip,ik),sigmapi, &
         Peta1,colat1d(idim),sigmapa(idim),sigmaha(idim),colata(idim)

  dphi2=dphi*2.       ! parameter for E field calculation
  Ecri=35.            ! critical E field value in mV/m
  AL=100.             ! arbitrary initial value
  if (nAE.gt.0) then
     call lintp(tAE,AEa,nAE,t,AE)
     call lintp(tAE,ALa,nAE,t,AL)
  endif

! Setup for Hill's and Weimer's models 
      delmlt=0.
      c0=0.77           ! a magic no. see Ober et al. (2003, JGR, 108(A12))
      UseAL=.false. 
      if (nAE.gt.0) UseAL=.true.  
      Tilt=psi1*180./pi
      call lintp(tsw,vswa,nsw,t,SWVel)
      call lintp(tsw,xnswa,nsw,t,SWDen)
      call lintp(timf,bxw,nimf,t,bx)
      call lintp(timf,byw,nimf,t,by)
      call lintp(timf,bzw,nimf,t,bz)
      SWBz=bz
      angle=atan2(by,bz)*180./pi
      Bt=sqrt(bx*bx+by*by+bz*bz)
      theta=acos(bz/Bt)
      sina1=sin(0.5*theta)         
      call SetModel(angle,Bt,Tilt,SWVel,SWDen,AL,UseAL)    

! Setup the baseline conductance when t=tstart if ipot.ge.3
  if (t.eq.tstart.and.ipot.ge.3) then
     kp1=2
     call condtot(kp1,ExAC)   ! quiet time conductance
     pedpsiB(:,:)=pedpsi(:,:)
     pedlamB(:,:)=pedlam(:,:)
     hallB(:,:)=hall(:,:)
  endif

! Calculate conductances when ipot eq 2
  if (ipot.eq.2) then   
     kp1=nint(zkp)
     if (kp1.lt.1) kp1=1
     if (kp1.gt.n_kp) kp1=n_kp
     call condtot(kp1,ExAC)  
  endif

! find conductances if ipot ge 3
  if (ipot.ge.3) then
     do j=1,jdim
        jc=j-jdif      ! j index in cimi
        if (jc.lt.1) jc=jc+ip
        if (jc.gt.ip) jc=jc-ip
        jm=j
        if (j.eq.1) jm=jdim-2
        if (j.eq.jdim) jm=3
        birk1d(:)=-1.*birk(:,j)  ! -1 for conversion between rcm and Robinson
        colat1d(:)=colat(:,j)
        colata(:)=colat1d(1)+(colat1d(:)-colat1d(1))/ExAC
        ! find sigmap and sigmah by R2 current or e- precipitation
        do i=1,idim
           ic=idim+1-i    ! i index in cimi
           if (ipot.eq.3) then   ! use J|| driven conductance
              imlt=nint(12.+aloct(i,j)*12./pi)  ! imlt in hour
              if (imlt.gt.24) imlt=imlt-24
              colatp1=colat(i,j)+pi/180.         
              facm1=0.  
              if (colatp1.le.colat(idim,j)) &
                  call lintp(colat1d,birk1d,idim,colatp1,facm1) 
              colatm1=colat(i,j)-pi/180.         
              facp1=0.  
              if (colatm1.ge.colat(1,j)) & 
                  call lintp(colat1d,birk1d,idim,colatm1,facp1) 
              fac0=birk1d(i)
              call conductance_model_v2(fac0,facm1,facp1,imlt,sigmap,sigmah)
           else              ! use e- and H+ precipitation driven conductance
              ! contribution from e-
              Etot=0.
              Ptot=0.
              Eflux=0.
              do k=1,je
                 if (gride_e(k).lt.40.) then
                    Etot=Etot+prePe(ic,jc,k)*gride_e(k)
                    Ptot=Ptot+prePe(ic,jc,k)
                    Eflux=Eflux+preFe(ic,jc,k)
                 endif
              enddo
              meanE=0.
              if (Ptot.gt.0.) meanE=Etot/Ptot
              EEF=40.*meanE*sqrt(0.5*Eflux)
              sigmap=2.*EEF/(16.+meanE*meanE)
              sigmah=0.45*sigmap*meanE**0.85
              ! add contribution from H+ (Eqs 6-7, Galand & Richmond, 2001)
              Etot=0.
              Ptot=0.
              Eflux=0.
              do k=1,je
                 if (gride_i(k).lt.100.) then
                    Etot=Etot+prePi(ic,jc,k)*gride_i(k)
                    Ptot=Ptot+prePi(ic,jc,k)
                    Eflux=Eflux+preFi(ic,jc,k)
                 endif
              enddo
              meanE=0.
              if (Ptot.gt.0.) meanE=Etot/Ptot
              sigmapi=5.7*Eflux**0.5
              sigmap=sigmap+sigmapi
              sigmah=sigmah+0.45*sigmapi*meanE**0.3
           endif
           sigmapa(i)=sigmap
           sigmaha(i)=sigmah
        enddo
        do i=1,idim
           ! map sigmap and sigmah form colata to colat1d grid
           call lintp(colata,sigmapa,idim,colat1d(i),sigmap)
           call lintp(colata,sigmaha,idim,colat1d(i),sigmah)
           ! Find sigmas
           sigmap2I=sigmap/(sini(i,j)**2)
           sigmah2I=sigmah/sini(i,j)
           pedpsi(i,j)=pedpsiB(i,j)+sigmap
           pedlam(i,j)=pedlamB(i,j)+sigmap2I
           hall(i,j)=hallB(i,j)+sigmah2I
        enddo   
     enddo      ! end of do j=1,jdim
  endif         ! end of if (ipot.ge.3)

! Find ionospheric potentials at the polarward boundary if ipot.ge.2
   if (ipot.ge.2) then       ! use self-consistent E field (SCE)
      ibmax=maxval(iba)
      do j=1,jdim
         ! determine ain
         if (iain.eq.1) ain(j)=1.               ! high lat boundary at xlati(ir)
         if (iain.eq.2) ain(j)=1.*(ir-ibmax+1)  ! high lat boundary@xlati(ibmax)
         i=ifix(ain(j))
         xlat1=pi/2.-colat(i,j)
         xmlt1=12.+aloct(i,j)*12./pi
         if (xmlt1.gt.24.) xmlt1=xmlt1-24.
         ! find potentials at the polarward boundary, vpob
         if (ihigh.eq.0) then                       ! Weimer's
            gLAT=acos(cos(xlat1)/sqrt(rc))*180./pi  ! invariant lat in deg
            gMLT=xmlt1                              ! mlt in hour 
            BnLat=BoundaryLat(gMLT)
            if (gLAT.LE.BnLat) vpob(j)=0.0
            if (gLAT.gt.BnLat) vpob(j)=EpotVal(gLAT,gMLT)*1000.  ! in Volt
         endif                                 
         if (ihigh.eq.1.and.j.eq.1) then     ! Boyle's
            vdrop=1000.*0.5*(1.1e-4*SWVel*SWVel+11.1*Bt*sina1**3)
            call hilat_bndy(vdrop,delmlt)    ! vpob is calculated
         endif                           
         if (ihigh.eq.2.and.j.eq.1) then     ! Hill's
            Psw=1.67e-27*SWden*1.e6*SWvel*SWvel*1.e6*1.e9  ! Psw(nPa)
            Esw=SWvel*1.e3*sqrt(by*by+bz*bz)*1.e-9*1.e3    ! Esw (mV/m)
            phim=30.+57.6*Esw/(Psw**(1./6.))*sin(theta/2.)
            sigmaHi=c0*sqrt(F107)
            phisa=(1600.*Psw**(1./3.))/(sigmaHi)
            pcp_Hill=(phim*phisa)/(phim+phisa)  ! PCP in kV
            vdrop=0.5*pcp_Hill*1000.            ! half of pcp in Volt
            call hilat_bndy(vdrop,delmlt)       ! vpob is calculated
         endif
         if (ihigh.eq.3) then                 ! use data from hybrid code
         endif
      enddo               ! end of do j=1,jdim 
   endif                  ! end of if (ipot.ge.2) 

! Compute ionospheric potential
  if (ipot.ge.2) then
     ! Calculate Peta
     do i=1,ir
        do j=1,ip
           do m=1,ik
              Peta(i,j,m)=0.
              do n=1,ijs
     !           if (js(n).ne.4) then      ! only count ions
                 do k=1,iw
                    Peta1=pcEo(n,k)*f2(n,i,j,k,m)*dlnp(n)*dlnK/ksai(i,j)
                    Peta(i,j,m)=Peta(i,j,m)+Peta1
                 enddo
      !          endif
              enddo
              Peta(i,j,m)=Peta(i,j,m)*1.6e-16    ! *1.6e-16 because pcEo in keV
           enddo
        enddo
     enddo
     ! Calculate birk_f (in microAmp/m^2)
     dvarL2=2.*dvarL
     birk_f(:,:)=0.
     do j=1,ip
        j0=j-1
        if (j0.lt.1) j0=j0+ip
        jc=j+1
        if (jc.gt.ip) jc=jc-ip
        do i=2,iba(j)-1          ! Set birk_f(1,j)=birk_f(iba-1,j)=0
           Brik=1.e6*BriN(i,j)/ksai(i,j)  ! *1.e6 because birk_f in microAmp/m^2
           if (i.le.iba(j0).and.i.le.iba(jc)) then
              do m=1,ik
                 dPdL=(Peta(i+1,j,m)-Peta(i-1,j,m))/dvarL2
                 dPdp=(Peta(i,jc,m)-Peta(i,j0,m))/dphi2
                 birk1=dlnBmdp(i,j,m)*dPdL-dlnBmdL(i,j,m)*dPdp
                 birk_f(i,j)=birk_f(i,j)+birk1*Brik   
              enddo
           endif
        enddo
     enddo
     ! Setup birk (in RCM grid)
     do j=1,jdim
        jc=j-jdif      ! j index in cimi
        if (jc.lt.1) jc=jc+ip
        if (jc.gt.ip) jc=jc-ip
        do i=1,idim
           ic=idim+1-i
           birk(i,j)=birk_f(ic,jc)        ! birk in Amp/km^2=microAmp/m^2
        enddo
     enddo
     call tube_volume(ro,xmlto,volume)
     call computev(potent)
 !   if (t.eq.tstart) then
 !      do n=1,ijs
 !         if (js(n).eq.4) write(*,*) 'Electrons are not included in SCE'
 !      enddo
 !   endif
  endif

! fill potential beyond ain 
      if (ipot.eq.1) ain(1:jdim)=float(ir+1)   ! use Weimer potential everywhere
      do j=3,jdim-1
         jc=j-jdif      ! j index in cimi
         if (jc.lt.1) jc=jc+ip
         if (jc.gt.ip) jc=jc-ip
         i0=ir+2-ifix(ain(j))
         gMLT=xmlt(jc)                      ! mlt in hour 
         do i=i0,ir
            if (ihigh.eq.0) then            ! Weimer model
               gLAT=acos(cos(xlati(i,jc))/sqrt(rc))*180./pi  ! invariant lat in deg
               BnLat=BoundaryLat(gMLT)
               if (gLAT.LE.BnLat) potent(i,jc)=0.0
               if (gLAT.gt.BnLat) potent(i,jc)=EpotVal(gLAT,gMLT)*1000.  ! Volt
            else                            ! Hill or Boyle's model
               phi1=(gMLT-delmlt)*pi/12     ! mlt in radian
               theta1=pi/2.-xlati(i,jc)
               thetab=pi/2.-xlati(i0-1,jc)
               potent(i,jc)=vdrop*sin(phi1)*sin(theta1)/sin(thetab)
            endif
         enddo
      enddo

! find potenth at xlath
  do i=0,irh
     do j=1,ip
        if (i.eq.0) xlathd=xlatd(ir+1,j)*180./pi
        if (i.gt.0) xlathd=xlath(i)
        gMLT=xmlt(j)
        gLAT=acos(cosd(xlathd)/sqrt(rc))*180./pi  ! invariant lat in deg
        BnLat=BoundaryLat(gMLT)
        if (gLAT.LE.BnLat) potenth(i,j)=0.0
        if (gLAT.gt.BnLat) potenth(i,j)=EpotVal(gLAT,gMLT)*1000.  ! Volt
     enddo
  enddo

  end subroutine convection


!***************************************************************************
!                               VdriftB
!  Routine calculates magnetic drifts 
!***************************************************************************
      subroutine VdriftB(ijs,js)
                     
      use constants
      use cread1
      use cVdrift     
      use cgrid
      use cfield
      implicit none
      real ksai1,Eo,pcEo1
      integer js(ns),n,i,j,k,m,i0,i1,j0,j1,icharge,ijs

      do n=1,ijs
         if (js(n).ge.5) icharge=-1          ! electrons
         if (js(n).lt.5) icharge=1           ! ions       
         Eo=xmass(n)*Em_speed*Em_speed       ! rest energy
         do k=1,iw
            pcEo1=pcEo(n,k)*1000./icharge    ! *1000 because pcEo in keV
            do m=1,ik
               do i=1,ir
                  i0=i-1
                  if (i0.lt.1) i0=1
                  i1=i+1
                  if (i1.gt.ir) i1=ir
                  do j=1,ip
                     j0=j-1
                     if (j0.lt.1) j0=j0+ip
                     j1=j+1
                     if (j1.gt.ip) j1=j1-ip
                     ! calculate vlB at (i+0.5,j) and vpB at (i,j+0.5)
                     ksai1=0.5*(ksai(i1,j)+ksai(i,j)) 
                     if (i1.le.iba(j0).and.i1.le.iba(j1)) then
                        vlB(n,i,j,k,m)=-pcEo1*dlnBmdp1(i,j,m)/ksai1
                     else
                        vlB(n,i,j,k,m)=vlB(n,i-1,j,k,m)
                     endif
                     if (i1.le.iba(j1).and.i1.le.iba(j)) then
                        vpB(n,i,j,k,m)=pcEo1*dlnBmdL1(i,j,m)/ksai(i,j)
                     else
                        vpB(n,i,j,k,m)=vpB(n,i-1,j,k,m)  
                     endif
                  enddo
               enddo          ! end of i loop
            enddo             ! end of j loop
         enddo                ! end of k loop
      enddo                   ! end of do n=1,ijs

      end subroutine VdriftB


!***************************************************************************
!                               VdriftE
!  Routine calculates electric drift and the corresponding drift in lnp.
!***************************************************************************
      subroutine VdriftE(potent,ijs,js)

      use constants
      use cread1
      use cVdrift
      use cgrid
      use cfield
      use convect,only : potenth
      implicit none
      real potentr(ir+1,ip),potent(ir,ip),ksai1,dphi2, &
           pt0,pt1,dPdL,dPdp,vdBm,vdL,vdp,dvarL1
      integer js(ns),n,i,j,k,m,i0,i1,j0,j1,ijs

      dphi2=dphi*2.

! Setup potentr
      potentr(1:ir,1:ip)=potent(1:ir,1:ip)
      potentr(ir+1,1:ip)=potenth(0,1:ip)

      do i=1,ir
         i0=i-1
         if (i0.lt.1) i0=1
         i1=i+1
         dvarL1=(i1-i0)*dvarL
         do j=1,ip
            j0=j-1
            if (j0.lt.1) j0=j0+ip
            j1=j+1
            if (j1.gt.ip) j1=j1-ip

            ! calculate vlE at (i+0.5,j)
            ksai1=0.5*(ksai(i1,j)+ksai(i,j))
            pt0=0.5*(potentr(i1,j0)+potentr(i,j0))
            pt1=0.5*(potentr(i1,j1)+potentr(i,j1))
            vlE(i,j)=-(pt1-pt0)/dphi2/ksai1

            ! calculate vpE at (i,j+0.5)
            pt0=0.5*(potentr(i0,j1)+potentr(i0,j))
            pt1=0.5*(potentr(i1,j1)+potentr(i1,j))
            vpE(i,j)=(pt1-pt0)/dvarL1/ksai(i,j)

            ! calculate vlnp at (i,j)
            dPdL=(potentr(i1,j)-potentr(i0,j))/dvarL1
            dPdp=(potentr(i,j1)-potentr(i,j0))/dphi2
            vdL=-dPdp/ksai(i,j)
            vdp=dPdL/ksai(i,j)
            do m=1,ik
               if (i.le.iba(j0).and.i.le.iba(j1)) then
                  vdBm=vdL*dlnBmdL(i,j,m)+vdp*dlnBmdp(i,j,m)
                  vlnp(i,j,m)=vdBm/2.
               else
                  vlnp(i,j,m)=vlnp(i-1,j,m)      
               endif
            enddo     
         enddo          ! end of i loop
      enddo             ! end of j loop

      end subroutine VdriftE

!******************************************************************************
!                           initial
! Initially set up the distribution function (f2) 
!******************************************************************************
      subroutine initial(itype,t,outname,st2,rb,gride,xmass,js,ijs,init)
      use constants
      use cinitial     
      use cgrid,only: d4,gridp,ekev,ebound
      use cfield
      use closs        
      use conductance
      real gride1(0:je+1),f(ir,ip,iw,ik),xmass(ns),gride(ns,je)
      real,allocatable,dimension(:) :: roi,ei,eilog
      real,allocatable,dimension(:,:) :: fi,psdi 
      integer js(ns),init(ns)
      character outname*13,st2(ns)*2,header*80

      eout(:,:,:,:)=0.   
      esum3(:,:,:,:)=0.
      psum3(:,:,:,:)=0.
      esum(:,:,:,:)=0.
      psum(:,:,:,:)=0.
      xled(:,:,:,:)=0.
      xlee(:,:,:,:)=0.
      xlec(:,:,:,:)=0.
      xlel(:,:,:,:)=0.
      xPreN(:,:,:,:)=0.
      xPreS(:,:,:,:)=0.
      xleb(:,:,:,:)=0.
      pled(:,:,:,:)=0.
      plee(:,:,:,:)=0.
      plec(:,:,:,:)=0.
      plel(:,:,:,:)=0.
      pPreN(:,:,:,:)=0.
      pPreS(:,:,:,:)=0.
      pleb(:,:,:,:)=0.
      HRPe(:,:,:)=0.
      HRPi(:,:,:)=0.
      f2(:,:,:,:,:)=0.
      pedpsi(:,:)=0.   
      pedlam(:,:)=0.   
      hall(:,:)=0.   
      preFe(:,:,:)=0.
      prePe(:,:,:)=0. 
      preFi(:,:,:)=0.
      prePi(:,:,:)=0. 
      ib0(1:ip)=iba(1:ip)

! Initial setup
  If (itype.eq.1) then
     call calc_Lstar2(jnoon,Lstar,Lstar_max)
     do n=1,ijs
        open(unit=3,file=trim(outname)//st2(n)//'.ece',status='unknown')
        write(3,*) je,'           ! je, ebound(je)'
        write(3,'(8f9.3)') (ebound(n,k),k=1,je)
        close(3)
     enddo
  else
     ! Read energy change from drift, charge exchange, wave diffusion,
     ! losscone loss and so on
     open(unit=11,file=trim(outname)//'.le',status='old')
     read(11,'(a80)') header
     read(11,*) parmod0
     read(11,*) ib0 
     read(11,*) eout
     read(11,*) xled
     read(11,*) xlee
     read(11,*) xlec
     read(11,*) xlel
     read(11,*) xPreN
     read(11,*) xPreS
     read(11,*) xleb
     read(11,*) pled
     read(11,*) plee
     read(11,*) plec
     read(11,*) plel
     read(11,*) pPreN
     read(11,*) pPreS
     read(11,*) pleb
     read(11,*) HRPe
     read(11,*) HRPi
     read(11,*) pedpsi
     read(11,*) pedlam
     read(11,*) hall
     read(11,*) preFe
     read(11,*) prePe
     read(11,*) preFi
     read(11,*) prePi
     close(11)
  endif

! Setup conservative psd at t=tstart
  if (itype.eq.1) then     ! initial run
     do n=1,ijs
        open(unit=8,file='quiet'//st2(n)//'.fin',status='old')
        read(8,*) il,ie
        allocate (roi(il),ei(ie),eilog(ie),fi(il,ie),psdi(il,ie))
        read(8,*) iunit   ! 1=flux in (cm2 s sr keV)^-1, 2=in (cm2 s MeV)^-1
        read(8,*) roi
        read(8,*) ei      ! ei in keV
        nskip=0
        if (js(n).ge.5) nskip=1+abs(init(n))+abs(init(n))*ie    ! for e-
        do i=1,nskip
           read(8,'(a)') header
        enddo
        read(8,*) fi
        close(8)
        if(iunit.eq.2) fi(:,:)=fi(:,:)/4./pi/1000. !con.To(cm^2 s sr keV)^-1
        if (js(n).lt.5.and.init(n).eq.0) then       ! for ion only
           fi(:,:)=1.e-30
           write(*,*) ' initial distribution = 0, species = ',st2(n)
        endif
        e_rest=xmass(n)*EM_speed*EM_speed/1.6e-16   ! rest energy in keV
        do k=1,ie
           eilog(k)=log10(ei(k))
           eir=ei(k)/e_rest             ! normalized energy in rest energy
           pp=xmass(n)*EM_speed*sqrt(eir*(eir+2.))
           do i=1,il
              psdi(i,k)=log10(fi(i,k)*6.25e19/pp/pp)      ! log(psd)
           enddo
        enddo

        do j=1,ip
           do i=1,iba(j)
              do m=1,ik
                 roii=Lstar(i,j,m)
	         if (roii.lt.roi(1)) roii=roi(1) ! flat dist @ lowL
	         if (roii.gt.roi(il)) roii=roi(il) ! flat @ high L

                 do k=1,iw
                    psd=0.
                    e1=log10(ekev(n,k)) 
                    if (e1.lt.eilog(1)) e1=eilog(1)   ! flat psd @ low E
                    if (e1.gt.eilog(ie)) e1=eilog(ie)  ! flat psd @ high E
                    call lintp2(roi,eilog,psdi,il,ie,roii,e1,x)
                    psd=10.**x         
                    f2(n,i,j,k,m)=psd*xjac(n,i,j,k,m) ! conser f2
                 enddo                                ! end of k loop
              enddo                                   ! end of m loop
  
           enddo                                      ! end of i loop
        enddo                                         ! end of j loop
        deallocate (roi,ei,eilog,fi,psdi)
        if (init(n).lt.0) then     ! add local injection
           call local_inj(f2,Tbounce,ekev,gridp,ro,xmlto,js,ijs,iba)
        endif
     enddo                                            ! end of n loop

  else                     ! read in flux from previous run

     open(unit=32,file=trim(outname)//'_c.f2',status='old',form='unformatted')
     do n=1,ijs
        read(32) f 
        f2(n,:,:,:,:)=f(:,:,:,:)           
     enddo         
     close(32)               

  endif            ! end of if (itype.eq.1)

! Find rbsum and rcsum           
      gride1(0)=0.
      gride1(je+1)=1.e10     ! arbitrary large number
      do n=1,ijs
         rbsum(n)=0.
         rcsum(n)=0.
         gride1(1:je)=ebound(n,1:je)
         do j=1,ip
            do i=1,iba(j)
               ii=2
               if (ro(i,j).le.xL1) ii=1
               if (ro(i,j).ge.xL2) ii=3
               do m=1,ik
                  do k=1,iw
                     ekev1=ekev(n,k)
                     weight=f2(n,i,j,k,m)*d4(n)
                     weighte=weight*ekev1
                     esum(n,i,j,je+2)=esum(n,i,j,je+2)+weighte
                     psum(n,i,j,je+2)=psum(n,i,j,je+2)+weight
                     psum3(n,ii,j,je+2)=psum3(n,ii,j,je+2)+weight
                     esum3(n,ii,j,je+2)=esum3(n,ii,j,je+2)+weighte
                     rbsum(n)=rbsum(n)+weighte
                     if (ro(i,j).le.6.6) rcsum(n)=rcsum(n)+weighte
                     kkloop: do kk=1,je+1
                        if (ekev1.gt.gride1(kk-1).and.ekev1.le.gride1(kk)) then
                           psum(n,i,j,kk)=psum(n,i,j,kk)+weight
                           esum(n,i,j,kk)=esum(n,i,j,kk)+weighte
                           psum3(n,ii,j,kk)=psum3(n,ii,j,kk)+weight
                           esum3(n,ii,j,kk)=esum3(n,ii,j,kk)+weighte
                           exit kkloop
                        endif
                     enddo kkloop
                  enddo
               enddo
            enddo
         enddo
      enddo

! Setup xPre0, pPre0, HRPe0 and HRPi0 
  xPre0N=xPreN
  xPre0S=xPreS
  pPre0N=pPreN
  pPre0S=pPreS
  HRPe0=HRPe
  HRPi0=HRPi

  end subroutine initial


!*******************************************************************************
  subroutine local_inj(f2,Tbounce,ekev,gridp,ro,xmlto,js,ijs,iba)
!*******************************************************************************
! Routine initializes a local injection of electrons.
! input: Tbounce,ekev,gridp,ro,xmlto,js,ijs,iba
! input/output: f2

  use constants
  use cimigrid_dim
  use cfield,only: xjac
  implicit none
  integer n,i,j,k,m,ijs,js(ns),iba(ip)
  real f2(ns,ir,ip,iw,ik),ekev(ns,0:iw+1),ro(ir,ip),xmlto(ir,ip),E1000,E2000, &
      gridp(ns,0:iw+1),Tbounce(ns,ir,ip,iw,ik), &
      Rp,LTp,dro,dmlt,dro2,dmlt2,rrp,rrpe,LLp,LLpe,rLpe,rLpee,Jpo,Jp,df2

! parameters for local injection 
  Rp=2.1
  LTp=8.0
  dro=0.18
  dmlt=0.5

! Add localized injection
  dro2=dro*dro
  dmlt2=dmlt*dmlt
  do n=1,ijs
     if (js(n).eq.5) then      ! high-energy electrons only
        do j=1,ip
           do i=1,iba(j)
              rrp=ro(i,j)-Rp
              rrpe=rrp*rrp/dro2/2.
              LLp=abs(xmlto(i,j)-LTp)
              if (LLp.gt.12.) LLp=24.-LLp
              LLpe=LLp*LLp/dmlt2/2.
              rLpe=rrpe+LLpe
              rLpee=0.
              if (rLpe.lt.80.) rLpee=exp(-rLpe)
              mloop: do m=1,ik
                 if (rLpee.eq.0.) exit mloop
                 do k=1,iw
                    E1000=ekev(n,k)/1000.
                    E2000=E1000*(0.575+0.055*E1000)
                    Jpo=2.2e5/exp(E2000)/Tbounce(n,i,j,k,m)
                    Jp=Jpo*rLpee/4./pi      ! flux in cm-2 sr-1 s-1 keV-1
                    df2=6.25e19*Jp*xjac(n,i,j,k,m)/gridp(n,k)/gridp(n,k)
                    f2(n,i,j,k,m)=f2(n,i,j,k,m)+df2
                 enddo
              enddo mloop
           enddo            ! i loop
        enddo               ! j loop
     endif
  enddo                     ! n loop

  end subroutine local_inj


!*******************************************************************************
      subroutine boundary(t,f2)
!*******************************************************************************
! Rountine find the instantaneous boundary conditions on the model boundary

  use constants
  use cbound
  use cfield
  use cread2
  use cgrid
  use dub_plasma
  implicit none
  integer iHi,iOi,ib,ib1,n,i,j,m,k
  real f2(ns,ir,ip,iw,ik),zkappa,zk1,zk2,gamma_dif,v2n,Psw,xmltj, &
       gamma_ra,t1,t,xnsw,vsw,Bz,Rtsy,phit,O_H_ratio,O_H_ratio0,N_H_ratio,chmass, &
       parE1,perE1,den1,temp32,fbb,fbb1,y2,x2,Vsq,Vpar2,Vper2,erpp,Psq,Ppar2, &
       Pper2,xn1(ip),zkt1(ip),xnn(ijs,ip),zktn(ijs,ip),logHe,HeFactor, &
       Rdub,zkte,xnne,BsAveN,BsAveT,BnAveT,NswAveN,VswAveT
  real gammln

  iHi=0
  iOi=0

! Use Ebihara-Ejiri-Borovsky or Tsyganenko-Mukai model for boundary n and kT
  t1=t
  if (iplsh.eq.1) t1=t-2.*3600.         ! 2 hours lag from SW to PS in Ebihara
  call lintp(tsw,xnswa,nsw,t1,xnsw)
  call lintp(tsw,vswa,nsw,t1,vsw)       ! vsw in km/s
  call lintp(timf,bzw,nimf,t,Bz)
  v2n=xnsw*vsw*vsw
  Psw=xmp*v2n*1.e12/1.e-9    ! Pdyn in nPa
  if (iplsh.eq.1) then
     xn1(1:ip)=0.025*xnsw+0.395         ! xnn in cm^-3, Ebihara & Ejiri 2000
     zkt1(1:ip)=0.019*vsw-3.65          ! kT_ion in keV, Borovsky etal 1998
  else
     do j=1,ip                          ! Tsyganenko-Mukai PS model
        ib=min(ir,iba(j)+1)
        Rtsy=ro(ib,j)
        if (Rtsy.lt.10.) Rtsy=10.   ! model only for r > 10
        phit=xmlto(ib,j)*pi/12.
        call get_tsy_plasma(Bz,vsw,xnsw,Rtsy,phit,zkt1(j),xn1(j))
     enddo
  endif

! Find average Bs, Bn, Nsw and Vsw if eplsh=3 (Dubyagin's e- plasmasheet model)
  if (eplsh.eq.3) then
     t1=t
     call lintp(timf,BsAveNa,nimf,t1,BsAveN)
     call lintp(timf,BsAveTa,nimf,t1,BsAveT)
     call lintp(timf,BnAvea,nimf,t1,BnAveT)
     call lintp(tsw,NswAvea,nsw,t1,NswAveN)
     call lintp(tsw,VswAvea,nsw,t1,VswAveT)
  endif

! Find the density and mean energy (temperature) at the boundary         
  if (icom.eq.2) then   ! calculate Psw with no time delay
     call lintp(tsw,xnswa,nsw,t,xnsw)
     call lintp(tsw,vswa,nsw,t,vsw)       ! vsw in km/s
     v2n=xnsw*vsw*vsw
     Psw=xmp*v2n*1.e12/1.e-9    ! Pdyn in nPa
  endif
  do n=1,ijs
     xnn(n,1:ip)=xn1(1:ip)
     zktn(n,1:ip)=zkt1(1:ip)
     if (js(n).eq.3) then               ! He+
        logHe=4.5e-3*(F107-50.)-2.4  ! log[He+] (Figure 12b)
        xnn(n,1:ip)=10.**logHe       ! [He+] in cm^-3
        if (icom.eq.2) then          ! add Pandya et al 2018
           HeFactor=exp(0.085598*Psw-4.365)    ! Fig 7c
           do j=1,ip
              xnn(n,j)=xnn(n,j)+xn1(j)*HeFactor
           enddo
        endif
     endif
     if (js(n).ge.5) then               ! electrons
        if (eplsh.eq.1) &   
           zktn(n,1:ip)=zktn(n,1:ip)/TiTe  ! Te, Ti/Te=TiTe (Baumjohann etal 89)
        if (eplsh.eq.2) then      ! RBE boundary condition for e-
           xnn(n,1:ip)=(0.02*xnsw+0.316)*sqrt(xmass(n)/xmp)   ! xnn in cm^-3
           zktn(n,1:ip)=0.016*vsw-2.4                         ! e- Eo in keV
        endif
        if (eplsh.eq.3) then      ! use Dubyagin's electron plasmasheet model
           do j=1,ip
              ib=min(ir,iba(j)+1)
              Rdub=ro(ib,j)
              if (Rdub.lt.6.) Rdub=6.     ! model only for 6 < r < 11
              if (Rdub.gt.11.) Rdub=11.   ! 
              phit=xmlto(ib,j)*15.        ! MLT in degree
              call get_dub_plasma(BsAveN,BsAveT,BnAveT,NswAveN,VswAveT, &
                                  Rdub,phit,zktn(n,j),xnn(n,j))
           enddo
        endif
     endif
     if (js(n).eq.1) iHi=1
     if (js(n).eq.2) iOi=1
     if (zktn(n,1).le.0.) then
        write(*,*) ' zktn(n,1).le.0. '
        stop
     endif
  enddo

! Redo xnn if consider both H+ and O+
   if (iHi.eq.1.and.iOi.eq.1) then
      O_H_ratio0=4.5e-2*exp(0.17*zkp+0.01*F107)                 ! Eq5,Young etal
      if (icom.eq.2) O_H_ratio0=O_H_ratio0+exp(0.23425*Psw-3.017)  ! Fig7a,Pandya
      O_H_ratio=0.67*O_H_ratio0
      N_H_ratio=0.33*O_H_ratio0
      do n=1,ijs
         if (js(n).eq.1) xnn(n,1:ip)=xnn(n,1:ip)/(O_H_ratio0+1.)             ! H+
         if (js(n).eq.2) xnn(n,1:ip)=xnn(n,1:ip)*O_H_ratio/(O_H_ratio0+1.)   ! O+
         if (js(n).eq.4) xnn(n,1:ip)=xnn(n,1:ip)*N_H_ratio/(O_H_ratio0+1.)   ! N+
      enddo
   endif

!  Calculate psd at boundary
   do n=1,ijs
      zkappa=ibset(n)*1.
      zk1=zkappa+1.
      zk2=zkappa-0.5
      gamma_dif=gammln(zk1)-gammln(zk2)
      gamma_ra=exp(gamma_dif)            ! gamma(kappa+1)/gamma(kappa-0.5)
      chmass=1.6e-19/xmass(n)
      do j=1,ip
         ib=iba(j)
         ib1=ib+1
         if (ib1.gt.ir) ib1=ir
         den1=0.
         parE1=0.
         perE1=0.
         if (xmlto(ib1,j).lt.6..or.xmlto(ib1,j).gt.18.) then  ! fb ne 0 at nite
            den1=xnn(n,j)*1.e6                       ! density in m-3
            parE1=zktn(n,j)                          ! characteristic in keV
            perE1=zktn(n,j)                          ! characteristic in keV
            temp32=sqrt(parE1)*perE1*1.6e-16**1.5
            if (ibset(n).eq.99) fbb1=den1/temp32/(2.*pi*xmass(n))**1.5 ! Mxwelin
            if (ibset(n).lt.99) fbb1=den1*gamma_ra/temp32/&
                                    (2.*pi*zkappa*xmass(n))**1.5       ! Kappa
         endif
         dena(n,j)=den1/1.e6                      ! density in cm-3
         parE(n,j)=parE1
         perE(n,j)=perE1
         do m=1,ik
            y2=y(ib,j,m)*y(ib,j,m)
            x2=1.-y2
            do k=1,iw
               fbb=0.
               if (den1.gt.0.) then
                  if (ibset(n).eq.99) then        ! Maxwellian
                     Vsq=vel(n,k)*vel(n,k)
                     Vpar2=Vsq*x2
                     Vper2=Vsq*y2
                     erpp=(Vpar2/parE1+Vper2/perE1)/2./1000./chmass
                     fbb=0.
                     if (erpp.lt.500.) fbb=fbb1/exp(erpp)
                  endif
                  if (ibset(n).lt.99) then
                     Psq=gridp(n,k)*gridp(n,k)
                     Ppar2=Psq*x2
                     Pper2=Psq*y2
                     erpp=(Ppar2/parE1+Pper2/perE1)/2./xmass(n)/1.6e-16
                     fbb=fbb1/(1.+erpp/zkappa)**zk1
                  endif
               endif      ! end of den1.gt.0
               fb(n,j,k,m)=fbb
               do i=ib+1,ir
                  f2(n,i,j,k,m)=fbb*xjac(n,i,j,k,m) ! f2 beyond rb
               enddo
            enddo                                   ! end of k loop
         enddo                                      ! end of m loop
      enddo                                         ! end of j loop

   enddo                                       ! end of n loop

   end subroutine boundary


!*******************************************************************************
  subroutine geoboundary(t,f2)
!*******************************************************************************
! Routine setup boundary condition at the geosyn orbit.

  use constants
  use cread2
  use cbound
  use cgrid,only: dlnp,gridp,ekev
  use cfield,only: iba,xmlto,dmu,xjac,F107,zkp
  implicit none
  
  integer iOi,iele,n,i,j,k,m,ntg1,ib
  real f2(ns,ir,ip,iw,ik),H_factor,O_factor,O_H_ratio,f_factor,MLT1D(nMLTg+1),&
       t,iflux1D(nMLTg+1),eflux1D(nMLTg+1),iflux2D(ip,iw),eflux2D(ip,iw), &
       flux2D(ip,iw),xmlto1,pre1(ijs,ip),dlnp4,fbb,weight,p2dp(ijs,iw),Bz, &
       xnsw,vsw,v2n,Psw

! Determine whether oxygen and electrons are included
  if (t.eq.tstart) then
     iOi=0
     iele=0
     do n=1,ijs
        if (js(n).eq.2) iOi=1
        if (js(n).ge.5) iele=1
     enddo 
  endif

! Determine H_factor and O_factor
  if (icom.eq.2) then
     call lintp(timf,bzw,nimf,t,Bz)
     call lintp(tsw,xnswa,nsw,t,xnsw)
     call lintp(tsw,vswa,nsw,t,vsw)       ! vsw in km/s
     v2n=xnsw*vsw*vsw
     Psw=xmp*v2n*1.e12/1.e-9    ! Pdyn in nPa
  endif
  H_factor=1.
  if (iOi.eq.1) then
     O_H_ratio=4.5e-2*exp(0.17*zkp+0.01*F107)                  ! Young etal 1982
     if (icom.eq.2) O_H_ratio=O_H_ratio+exp(0.23425*Psw-3.017)    ! Fig7a,Pandya
     H_factor=1./(O_H_ratio+1.)
     O_factor=O_H_ratio/(O_H_ratio+1.)
  endif

! Determine the time index of geosyn data
  if (t.eq.tstart) ntg0=-99.    ! arbitrary -ve number
  call locate1(tgeo,ntg,t,ntg1)
  if (ntg1.eq.0) ntg1=1
  if (ntg1.lt.ntg) then
     if ((tgeo(ntg1+1)-t).lt.(t-tgeo(ntg1))) ntg1=ntg1+1
  endif
  if (ntg1.eq.ntg0) return

! Setup 1D array of MLT
  MLT1D(1:nMLTg)=MLTgeo(ntg1,1:nMLTg)
  MLT1D(nMLTg+1)=MLT1D(1)+24.

! find ion and electron fluxes at geosyn
  do k=1,iw
     iflux1D(1:nMLTg)=iflux66(ntg1,1:nMLTg,k)
     iflux1D(nMLTg+1)=iflux1D(1)
     eflux1D(1:nMLTg)=eflux66(ntg1,1:nMLTg,k)
     eflux1D(nMLTg+1)=eflux1D(1)
     do j=1,ip
        ib=min(ir,iba(j)+1)
        xmlto1=xmlto(ib,j)
        if (xmlto1.lt.MLT1D(1)) xmlto1=xmlto1+24.
        call lintp(MLT1D,iflux1D,nMLTg+1,xmlto1,iflux2D(j,k))
        if (iele.eq.1) call lintp(MLT1D,eflux1D,nMLTg+1,xmlto1,eflux2D(j,k))
     enddo
  enddo

! Calculate psd, density and characteristic energy at rb
  dena(:,:)=0.
  pre1(:,:)=0.
  do n=1,ijs
     dlnp4=8.*pi*dlnp(n)/3.
     do k=1,iw
        p2dp(n,k)=(gridp(n,k)**3)*dlnp4
     enddo
     if (js(n).lt.5) then
        flux2D(:,:)=iflux2D(:,:)
        if (js(n).eq.1) f_factor=H_factor
        if (js(n).eq.2) f_factor=O_factor
        if (js(n).eq.3) then
           write(*,*) 'Error: He+ is not ready for geosynchronous boundary'
           stop
        endif
     else
        flux2D(:,:)=eflux2D(:,:)
        f_factor=1.
     endif
     do j=1,ip
        ib=min(ir,iba(j)+1)
        do k=1,iw
           fbb=f_factor*flux2D(j,k)/gridp(n,k)/gridp(n,k)/1.6e-20  ! psd
           do m=1,ik
              fb(n,j,k,m)=fbb
              weight=fbb*p2dp(n,k)*dmu(ib,j,m)
              dena(n,j)=dena(n,j)+weight
              pre1(n,j)=pre1(n,j)+ekev(n,k)*weight ! pressure
              do i=ib+1,ir
                 f2(n,i,j,k,m)=fbb*xjac(n,i,j,k,m) ! f2 beyond rb
              enddo
           enddo
        enddo
        parE(n,j)=0.
        if (dena(n,j).gt.0.) parE(n,j)=pre1(n,j)/dena(n,j)
        perE(n,j)=parE(n,j)             ! isotropic
        dena(n,j)=dena(n,j)*1.e-6       ! convert dena to cm-3
     enddo
  enddo

  ntg0=ntg1
 
  end subroutine geoboundary

!*******************************************************************************
!                           p_result
!          Routine prints all the results at selective times
!*******************************************************************************
   subroutine p_result(t,tstart,parE,perE,dena,ompe,CHpower, &
                       HIpower,nprint,js,ijs,ibset,tpls,ipot)

      use constants
      use cread1
      use cread2,only: tint,icoul,iEMIC,elon,ctp,stp
      use cgrid
      use cfield
      use cinitial
      use WaveGrowth
      use closs
      use rcmgrid
      use current
      use conductance
      use convect
      use cPlasmasphere_new, only: nlp,npp,Nion
      use ModCurvScatt, only: write_flc
      integer ibset(ns),js(ns)
      real preF(je+2),preP(je+2),HRPee(ijs,ir,ip),p2dp(ijs,iw), &
           pedpsi_f(ir,ip),HRPii(ijs,ir,ip), &
           f(ijs,ir,ip,iw,ik),dena(ns,ip),xlatiS1(ir,ip),xmltS(ir,ip), &
           xlati1(ir,ip),parE(ns,ip),perE(ns,ip),f5d(ijs,ir,ip,iw,ik), &
           ompe(ir,ip),CHpower(ir,ip),HIpower(ir,ip)

      hour=t/3600.
      tintf=tint

! Convert conjugate grid in degree and hour
  do i=1,ir
     do j=1,ip
        xlatiS1(i,j)=xlatiS(i,j)*180./pi   ! lat. at S ionosphere in degree
        xmltS(i,j)=phiS(i,j)*12./pi        ! mlt at S ionosphere in hour
     enddo
  enddo

! Calculate 4*pi*p^2*dp (or 4*pi*p^3*dlnp)
  do n=1,ijs
     dlnp4=4.*pi*dlnp(n)
     do k=1,iw
        p2dp(n,k)=(gridp(n,k)**3)*dlnp4
     enddo
  enddo

! Calculate fluxes
      f(:,:,:,:,:)=0.
      f5d(:,:,:,:,:)=0.
      do n=1,ijs
         do i=1,ir
            do j=1,ip

               do m=0,ik+1
               enddo

               denWP(n,i,j)=0.       
               TparaWP(n,i,j)=0.     
               TperpWP(n,i,j)=0.     
               PparaWP=0.    
               PperpWP=0.   
               do m=1,ik
                  sin_Sq=y(i,j,m)*y(i,j,m)
                  cos_Sq=1.-sin_Sq
                  do k=1,iw
                     xm1=gridp(n,k)/vel(n,k)
                     xm3=xm1**3
                     psd=f2(n,i,j,k,m)/xjac(n,i,j,k,m)  ! kg^-3m^-6s^3
                     f5d(n,i,j,k,m)=psd*xm3           ! psd in m^-6s^3
                     f(n,i,j,k,m)=psd*1.6e-20*gridp(n,k)*gridp(n,k) ! flux
                     weight=psd*p2dp(n,k)*dmu(i,j,m)
                     denWP(n,i,j)=denWP(n,i,j)+weight
                     PparaWP=PparaWP+2.*ekeV(n,k)*cos_Sq*weight
                     PperpWP=PperpWP+ekeV(n,k)*sin_Sq*weight
                  enddo
               enddo
               if (denWP(n,i,j).gt.0.) then
                  TparaWP(n,i,j)=PparaWP/denWP(n,i,j)
                  TperpWP(n,i,j)=PperpWP/denWP(n,i,j)
               endif

            enddo
         enddo
      enddo

! Get ain, conductance at cimi grid
  pedpsi_f=0.       !
  if (ipot.ge.2) then
      do j=3,jdim-1
         jc=j-jdif      ! j index in cimi
         if (jc.lt.1) jc=jc+ip
         if (jc.gt.ip) jc=jc-ip
         do i=ifix(ain(j)),idim
            ic=idim+1-i
            pedpsi_f(ic,jc)=pedpsi(i,j)
         enddo
      enddo
  endif

! Calculate heating rate per unit volume to plasmasphere if icoul eq 1
  HRPee(:,:,:)=0.
  HRPii(:,:,:)=0.
  if (icoul.eq.1) then
     do j=1,ip
        do i=1,iba(j)
           Vtube=volume(i,j)*ksai(i,j)*dphi*dvarL    ! flux tube volume in m^3
           Vtubet=Vtube*tint
           do n=1,ijs
              HRPee(n,i,j)=(HRPe(n,i,j)-HRPe0(n,i,j))/Vtubet   ! HRPee,HRPii in
              HRPii(n,i,j)=(HRPi(n,i,j)-HRPi0(n,i,j))/Vtubet   ! keV/m^3/sec
              if (HRPee(n,i,j).lt.1.e-30) HRPee(n,i,j)=0.
              if (HRPii(n,i,j).lt.1.e-30) HRPii(n,i,j)=0.
           enddo
        enddo
     enddo
     HRPe0=HRPe
     HRPi0=HRPi
  endif

! Calculate fluxes at fixed E and y grides and output to .fls
      do i=1,ir
         do j=1,ip
            xlati1(i,j)=xlati(i,j)*180./pi   ! lat. at ionosphere in degree
         enddo
      enddo
      do n=1,ijs
         if (t.eq.tstart) then
            open(unit=26,file=trim(outname)//st2(n)//'.fls',status='unknown')
            write(26,'(f10.6,4i6,a)') rc,ir,ip,je,ig, &
                  '  ! rc_Re,ir,ip,je,ig, no. of grid in Li,mlong,E,sinPAo'
            write(26,'(7f11.5)') (varL(i),i=1,ir)
            write(26,'(7f11.6)') (mphi(j),j=1,ip)
            write(26,'(7f11.5)') (gride(n,k),k=1,je)
            write(26,'(7f11.5)') (gridy(m),m=1,ig)
         else
            open(unit=26,file=trim(outname)//st2(n)//'.fls',status='old',&
                 position='append')
         endif
         if (mod(t,tintf).eq.0) then
            if (n.eq.1.and.t.gt.tstart) call calc_Lstar2(jnoon,Lstar,Lstar_max)
            call fluxes(f,hour,xlati1,xlatiS1,xmltS,parmod,HRPee,HRPii,Lstar, &
                        Lstar_max,ompe,CHpower,HIpower,js,ijs,n)
         endif
         close(26)
         if (n.eq.ijs) write(*,*) 'finish writing flux to .fls'
      enddo

! Write boundary Eo (keV) and density (cm-3)
      do n=1,ijs
         if (t.eq.tstart) then
            open(unit=12,file=trim(outname)//st2(n)//'.bc',status='unknown')
            if (ibset(n).eq.99)write(12,*) 'Maxwellian distribution at boundary'
            if (ibset(n).lt.99)write(12,*) 'Kappa distribution at boundary'
         else
            open(unit=12,file=trim(outname)//st2(n)//'.bc',status='old',&
                 position='append')
         endif
         write(12,*) hour,'       ! hour'
         write(12,*) &
       '        ro          mlto     Epar(keV)    Eper(keV)   density(cm-3)'
         do j=1,ip
            ib=min(ir,iba(j)+1)
            write(12,'(5f13.5)') rsb(j),xmltb(j),parE(n,j),perE(n,j),&
                                 dena(n,j)
         enddo
         close(12)
      enddo

! Write the energy changes in .ece fiels
   if (t.gt.tstart) then
      totbsumi=0.
      totbsum=0.
      do n=1,ijs
        open(unit=3,file=trim(outname)//st2(n)//'.ece',status='old', &
             position='append')
        call echannel(hour,Dst,eout,esum3,bsum,xL1,xL2,n)
        if (js(n).lt.5) totbsumi=totbsumi+bsum  ! don't count e-
        if (js(n).lt.6) totbsum=totbsum+bsum    ! don't count low-E e-
        close(3)
      enddo
   endif

!  Write potential, conductance, region 2 current into file *.pot
      if (t.eq.tstart) then
         open(unit=13,file=trim(outname)//'.pot',status='unknown')
         write(13,*) rc,'               ! radius of ionosphere in Re'
         write(13,*) ir,'               ! no. of latitude (deg) '
         write(13,*) irh,'              ! no. of latitude beyond xlati(ir)'
         write(13,*) ip,'               ! no. of local time (hour)'
         write(13,*) ik,'               ! no. of longitudial invariant'
         write(13,'(10f8.3)') xlati1    ! xlati1(ir,ip) in degree
         write(13,'(10f8.3)') xlath
         write(13,'(8f10.1)') potentc   ! Corotation potentials in Volt
         write(13,'(1p,7e11.3)') BriN   ! radial component of B at N ionosphere
      else
         open(unit=13,file=trim(outname)//'.pot',status='old',position='append')
      endif
      write(13,*) hour,'            ! time in hour'
      write(13,'(16i5)') iba
      write(13,'(10f8.2)') (xmlt(j),j=1,ip)
      write(13,'(10f8.3)') ((xlatiS1(i,j),i=1,ir),j=1,ip)
      write(13,'(10f8.2)') ((xmltS(i,j),i=1,ir),j=1,ip)
      write(13,'(8f10.2)') ((ro(i,j),i=1,ir),j=1,ip)
      write(13,'(8f10.2)') ((xmlto(i,j),i=1,ir),j=1,ip)
      write(13,'(1p,7e11.3)') ((bo(i,j),i=1,ir),j=1,ip)
      write(13,'(8f10.5)') (((y(i,j,m),i=1,ir),j=1,ip),m=1,ik)
      write(13,'(8f10.1)') ((potent(i,j),i=1,ir),j=1,ip)       ! in Volt
      write(13,'(8f10.1)') ((potenth(i,j),i=1,irh),j=1,ip)     ! in Volt
      write(13,'(8f10.2)') ((pedpsi_f(i,j),i=1,ir),j=1,ip)  ! both hemispheres
      write(13,'(1p,6e13.4)') ((birk_f(i,j),i=1,ir),j=1,ip)
      close(13)

! Write total energy in .db files
  totRCi=0.
  totEi=0.
  totRC=0.
  totE=0.
  do n=1,ijs
     if (js(n).lt.5) totRCi=totRCi+rcsum(n)     ! don't count e-
     if (js(n).lt.5) totEi=totEi+rbsum(n)       ! don't count e-
     if (js(n).lt.6) totRC=totRC+rcsum(n)       ! don't count low-E e-
     if (js(n).lt.6) totE=totE+rbsum(n)         ! don't count low-E e-
  enddo
  totE_Bi=totEi-totbsumi
  totE_B=totE-totbsum
  totRCe=totRC-totRCi
  totEe=totE-totEi
  totE_Be=totE_B-totE_Bi
  if (t.eq.tstart) then
    open(unit=2,file=trim(outname)//'.dbi',status='unknown')
    write(2,'(a)')'     hour    Kp      AE       AL       Dst     DstRC   '//&
               'rcsum       totE     totE-Bchange'
    open(unit=23,file=trim(outname)//'.dbe',status='unknown')
    write(23,'(a)')'     hour    Kp      AE       AL       Dst     DstRC   '//&
               'rcsum       totE     totE-Bchange'
    open(unit=7,file=trim(outname)//'.db',status='unknown')
    write(7,'(a)')'     hour    Kp      AE       AL       Dst     DstRC   '//&
               'rcsum       totE     totE-Bchange'
  else
    open(unit=2,file=trim(outname)//'.dbi',status='old',position='append')
    open(unit=23,file=trim(outname)//'.dbe',status='old',position='append')
    open(unit=7,file=trim(outname)//'.db',status='old',position='append')
  endif
  write(2,'(f9.3,f6.2,4f9.1,1p,3e12.4)') hour,zkp,AE,AL,Dst,DstRC,totRCi, &
                                         totEi,totE_Bi
  write(23,'(f9.3,f6.2,4f9.1,1p,3e12.4)') hour,zkp,AE,AL,Dst,DstRC,totRCe, &
                                         totEe,totE_Be
  write(7,'(f9.3,f6.2,4f9.1,1p,3e12.4)') hour,zkp,AE,AL,Dst,DstRC,totRC, &
                                         totE,totE_B
  close(2)
  close(23)
  close(7)

! Output precipitating energy and particle flux on both hemispheres
  preFe(:,:,:)=0.       ! energy flux of precipitating electrons in mW/m2
  prePe(:,:,:)=0.       ! number of precipitating electrons 
  preFi(:,:,:)=0.       ! energy flux of precipitating protons in mW/m2
  prePi(:,:,:)=0.       ! number of precipitating protons   
  rc2=rc*rc*re_m*re_m
  if (t.eq.tstart) then
     nprint1=nprint-1
     do n=1,ijs
        do m=1,2
           if (m.eq.1) open (unit=m,file=trim(outname)//st2(n)//'_N.preci')
           if (m.eq.2) open (unit=m,file=trim(outname)//st2(n)//'_S.preci')
           write(m,"(f10.6,5i6,6x,'! rc in Re,ir,ip,je,ik,ntime')") &
                                  rc,ir,ip,je,ik,nprint1
           write(m,"(f11.5,2f11.6,' ! elon,ctp,stp')") elon,ctp,stp
           write(m,'(8f9.3)') (gride(n,k),k=1,je)    ! mean E in each bin
           close(m)
        enddo
     enddo
  else
     do n=1,ijs
        do m=1,2
           if (m.eq.1) open (unit=m,file=trim(outname)//st2(n)//'_N.preci', &
                             status='old',position='append')
           if (m.eq.2) open (unit=m,file=trim(outname)//st2(n)//'_S.preci', &
                             status='old',position='append')
           write(m,*) hour,Dst,'       ! hour, Dst'
           area1=rc2*dphi
           do i=1,ir
              do j=1,ip
                 area=area1*cos(xlati(i,j))*dlati(i,j)            ! area in m^2
                 if (m.eq.2) area=area*BriN(i,j)/BriS(i,j)
                 Asec=area*tint
                 do k=1,je+2
                    if (m.eq.1) dPre=xPreN(n,i,j,k)-xPre0N(n,i,j,k)
                    if (m.eq.1) dpPre=pPreN(n,i,j,k)-pPre0N(n,i,j,k)
                    if (m.eq.2) dPre=xPreS(n,i,j,k)-xPre0S(n,i,j,k)
                    if (m.eq.2) dpPre=pPreS(n,i,j,k)-pPre0S(n,i,j,k)
                    preF(k)=0.
                    preP(k)=0.
                    if (dPre.lt.0..and.dpPre.lt.0.) then
                       preF(k)=-dPre*1.6e-13/Asec        ! E flux in mW/m2
                       preP(k)=-dpPre                    ! number of particles
                       if (js(n).eq.5.and.k.le.je) then
                          preFe(i,j,k)=preFe(i,j,k)+preF(k)
                          prePe(i,j,k)=prePe(i,j,k)+preP(k)
                       endif
                       if (js(n).eq.1.and.k.le.je) then
                          preFi(i,j,k)=preFi(i,j,k)+preF(k)
                          prePi(i,j,k)=prePi(i,j,k)+preP(k)
                       endif
                    endif
                 enddo
                 if (m.eq.1) write(m,'(f8.3,f7.2,f9.2,1pE12.4,i4,1pE12.4,a)') &
                    xlati1(i,j),xmlt(j),mlon(j),area,mcN(i,j),BiN(i,j), &
                    ' ! mlat,mlt,mlon,area,mc,Bi'
                 if (m.eq.2) write(m,'(f8.3,f7.2,f9.2,1pE12.4,i4,1pE12.4,a)') &
                    xlatiS1(i,j),xmltS(i,j),mlonS(i,j),area,mcS(i,j),BiS(i,j), &
                    ' ! mlat,mlt,mlon,area,mc,Bi'
                 write(m,'(1p,6e12.4)') preF
                 write(m,'(1p,6e12.4)') preP
              enddo
           enddo
           close(m)
        enddo        ! end of do m=1,2
     enddo           ! end of do n=1,ijs
  endif
  xPre0N=xPreN
  xPre0S=xPreS
  pPre0N=pPreN
  pPre0S=pPreS

! Calculate EMIC convective growth rates
  if (t.eq.tstart) then
     ! Setup Xpar(1:nw), parameter X in George's note, normalized frequency
     dXp=1./nw
     Xpar(1)=0.5*dXp
     do m=2,nw
        Xpar(m)=Xpar(m-1)+dXp 
     enddo    
     open(unit=9,file=trim(outname)//'.gamma')
     write(9,*) ir,ip,'     ! ir,ip '
     write(9,*) nw,'     ! nw, number of normalized frequency '
     write(9,'(10f8.4)') Xpar
  else
     open(unit=9,file=trim(outname)//'.gamma',status='old',position='append')
  endif
  if (t.ge.tpls.and.iEMIC.gt.0) then
     call MapDistribution(f5d)
     call EMICGrowthRate(hour,js,ijs)
  endif
  close(9)

! Open files to write fluxes, energy loss due to different processes and
! plasmasphere Nion for continuous run
      if (t.gt.tstart) then
         open(unit=32,file=trim(outname)//'_c.f2',form='unformatted')
         do n=1,ijs
            write(32) f2(n,:,:,:,:)     ! write f2   
         enddo
         close(32)

         open(unit=11,file=trim(outname)//'.le',status='unknown')
         write(11,*) hour,'   ! hour'
         write(11,*) parmod0
         write(11,*) ib0 
         write(11,*) eout
         write(11,*) xled
         write(11,*) xlee
         write(11,*) xlec
         write(11,*) xlel
         write(11,*) xPreN
         write(11,*) xPreS
         write(11,*) xleb
         write(11,*) pled
         write(11,*) plee
         write(11,*) plec
         write(11,*) plel
         write(11,*) pPreN
         write(11,*) pPreS
         write(11,*) pleb
         write(11,*) HRPe
         write(11,*) HRPi
         write(11,*) pedpsi
         write(11,*) pedlam
         write(11,*) hall   
         write(11,*) preFe
         write(11,*) prePe
         write(11,*) preFi
         write(11,*) prePi
         close(11)

         open(unit=31,file=trim(outname)//'_c.Nion',form='unformatted')
         write(31) Nion
         close(31)
      endif

      end subroutine p_result


!*******************************************************************************
  subroutine p_rtp(t)
!*******************************************************************************
! Routine finds the B fields and the equatorial crossing point of each point in 
! a spherical grid (rtp) in SM. This information, combined with .fls and .pot, 
! the differential flux of energetic particles, cold plasma density and
! potential at any point in the CIMI domain can be determined.
  
  use cread1
  use cread2, only: tstart,rb,imod,intB
  use cgrid, only: xlati,mphi
  use constants, only: pi
  use crtp
  implicit none
  integer i,j,m,iout
  real t,dir,rlim,RsinT,xi,yi,zi,xf,yf,zf,bx,by,bz,hour,fn(3),rxy
  real BXrtp(ir,nth,ip+1),BYrtp(ir,nth,ip+1),BZrtp(ir,nth,ip+1)
  real Rortp(ir,nth,ip+1),MLTrtp(ir,nth,ip+1)

! open .rtp file and write the spherical grid, rtp, when t=tstart
  if (t.eq.tstart) then
     open(unit=22,file=trim(outname)//'.rtp')
     write(22,*) ir,nth,ip+1,'    ! nr,nth,nphi, grid dimension in SM'
     write(22,'(7f11.5)') Rrtp
     write(22,'(7f11.6)') Trtp
     write(22,'(7f11.6)') Prtp
  else
     open(unit=22,file=trim(outname)//'.rtp',status='old',position='append')
  endif

! Find equatorial crossing point and B field at each SM rtp grid pt
  rlim=1.2*rb
  do i=1,ir
     do m=1,nth
        RsinT=Rrtp(i)*sin(Trtp(m))
        zi=Rrtp(i)*cos(Trtp(m))
        do j=1,ip
           xi=RsinT*cos(Prtp(j))
           yi=RsinT*sin(Prtp(j))
           ! find B at the grid point
           call Btotal_SM(xi,yi,zi,bx,by,bz,fn)
           BXrtp(i,m,j)=bx*1.e-9    ! BXrtp in Tesla
           BYrtp(i,m,j)=by*1.e-9    ! BYrtp in Tesla
           BZrtp(i,m,j)=bz*1.e-9    ! BZrtp in Tesla
           ! find equatorial crossing of the grid point
           call trace_SMeq(xi,yi,zi,rlim,xf,yf,zf,iout)
           MLTrtp(i,m,j)=atan2(yf,xf)*12./pi+12.   ! mlt in hr
           rxy=sqrt(xf*xf+yf*yf)
           Rortp(i,m,j)=99.
           if (iout.eq.0) Rortp(i,m,j)=rxy
        enddo        ! end of do j=1,ip
        Rortp(i,m,ip+1)=Rortp(i,m,1)
        MLTrtp(i,m,ip+1)=MLTrtp(i,m,1)
        BXrtp(i,m,ip+1)=BXrtp(i,m,1)
        BYrtp(i,m,ip+1)=BYrtp(i,m,1)
        BZrtp(i,m,ip+1)=BZrtp(i,m,1)
     enddo
  enddo

! Write the results
  hour=t/3600.
  write(22,*) hour,'    ! hour'
  write(22,'(1p,7e11.3)') BXrtp
  write(22,'(1p,7e11.3)') BYrtp
  write(22,'(1p,7e11.3)') BZrtp
  write(22,'(7f11.5)') Rortp
  write(22,'(7f11.6)') MLTrtp
  close(22)
  write(*,*) 'finish writing to .rtp'
         
  end subroutine p_rtp


!*******************************************************************************
  subroutine Btotal_SM(xsm,ysm,zsm,bxsm,bysm,bzsm,fn)
!*******************************************************************************
! Routine calculate the total field from the external fieid (T96 or T04) and 
! internal field (dipole or IGRF) in SM coor.
 
  use constants, only: i_one,m_one
  use cread2, only: intB,imod
  use cfield, only: parmod,psi1
  implicit none
  real xsm,ysm,zsm,bxsm,bysm,bzsm,xgsm,ygsm,zgsm,fn(3),btot
  real bxint,byint,bzint,bxext,byext,bzext,bxgsm,bygsm,bzgsm

  call smgsw_08(xsm,ysm,zsm,xgsm,ygsm,zgsm,i_one)
  if (intB.eq.0) call dip_08(xgsm,ygsm,zgsm,bxint,byint,bzint)
  if (intB.eq.1) call IGRF_GSW_08(xgsm,ygsm,zgsm,bxint,byint,bzint)
  if (imod.eq.0) call zeroB(imod,parmod,psi1,xgsm,ygsm,zgsm,bxext,byext,bzext)
  if (imod.eq.1) call t96_01(imod,parmod,psi1,xgsm,ygsm,zgsm,bxext,byext,bzext)
  if (imod.eq.2) call t04_s(imod,parmod,psi1,xgsm,ygsm,zgsm,bxext,byext,bzext)
  bxgsm=bxint+bxext
  bygsm=byint+byext
  bzgsm=bzint+bzext

! B in sm
  call smgsw_08(bxsm,bysm,bzsm,bxgsm,bygsm,bzgsm,m_one)

! unit vector of B in sm
  btot=sqrt(bxsm*bxsm+bysm*bysm+bzsm*bzsm)
  fn(1)=bxsm/btot     ! fn is unit vector
  fn(2)=bysm/btot
  fn(3)=bzsm/btot

  end subroutine Btotal_SM


!*******************************************************************************
  subroutine trace_SMeq(xi,yi,zi,rlim,xf,yf,zf,iout)
!*******************************************************************************
! Routine to do field line tracing in SM coordinates at ZSM=0 

  use cfield, only: dsmax
  implicit none
  integer,parameter :: np=1000
  integer iout,L,i
  real xi,yi,zi,rlim,xf,yf,zf,dir,t0,tend,ds,x0(3),xend(3),zlast,zcurrent
  real zsign,RYZ,R2,R,xx,yy,zz,wz

      if (zi.le.0.0) dir=1.0
      if (zi.gt.0.0) dir=-1.
      t0=0.
      ds=dsmax*dir
      x0(1)=xi
      x0(2)=yi
      x0(3)=zi
      zlast=x0(3)
      iout=0
      L=0

! start tracing
  1   L=L+1
      IF (L.GT.np) then
         L=np
         iout=1
         return        ! end tracing when L>np
      ENDIF
      RYZ=x0(2)*x0(2)+x0(3)*x0(3)
      R2=X0(1)*x0(1)+RYZ
      R=SQRT(R2)
      xx=x0(1)
      yy=x0(2)
      zz=x0(3)
      IF (R.GT.rlim.OR.RYZ.GT.1600..OR.x0(1).GT.20.) then
         iout=1
         GOTO 8        ! end tracing when reaches to far point
      ENDIF
      call rk4(x0,t0,ds,xend,tend)
      zcurrent=xend(3)
      zsign=zlast*zcurrent
      if(zsign.lt.0.0) goto 99   ! end tracing when crossing the equator
      do i=1,3  
         x0(i)=xend(i)
      enddo
      t0=tend
      zlast=xend(3)
      goto 1          ! continue tracing

! Adjust the end pt and move it to the equator
 99   x0(1)=xx
      x0(2)=yy
      x0(3)=zz
      wz=-xend(3)/(x0(3)-xend(3))
      xend(1)=xend(1)-(xend(1)-x0(1))*wz
      xend(2)=xend(2)-(xend(2)-x0(2))*wz
      xend(3)=xend(3)-(xend(3)-x0(3))*wz
              
  8   xf=xend(1)
      yf=xend(2)
      zf=xend(3)

  end subroutine trace_SMeq


!*******************************************************************************
  subroutine rk4(x0,t0,ds,xend,tend)
!*******************************************************************************
! *** FOURTH-ORDER RUNGE-KUTTA ***
! Solve xend = x0 + fn*ds               ! fn is the unit vector of f

  implicit none
  integer i
  real x0(3),xend(3),x00(3),xwrk(4,3),f(3),fn(3)
  real t0,tend,ds


  call Btotal_SM(x0(1),x0(2),x0(3),f(1),f(2),f(3),fn)
  do i=1,3 
     x00(i)=x0(i)
     xwrk(1,i)=ds*fn(i)
     xend(i)=x00(i)+xwrk(1,i)/2.
     x0(i)=xend(i)
  enddo
  call Btotal_SM(x0(1),x0(2),x0(3),f(1),f(2),f(3),fn)
  do i=1,3 
     xwrk(2,i)=ds*fn(i)
     xend(i)=x00(i)+xwrk(2,i)/2.
     x0(i)=xend(i)
  enddo
  call Btotal_SM(x0(1),x0(2),x0(3),f(1),f(2),f(3),fn)
  do i=1,3 
     xwrk(3,i)=ds*fn(i)
     xend(i)=x00(i)+xwrk(3,i)
     x0(i)=xend(i)
  enddo
  call Btotal_SM(x0(1),x0(2),x0(3),f(1),f(2),f(3),fn)
  do i=1,3 
     xwrk(4,i)=ds*fn(i)
  enddo
  do i=1,3 
     xend(i)=x00(i)+(xwrk(1,i)+2.*xwrk(2,i)+2.*xwrk(3,i)+xwrk(4,i))/6.
  enddo
  tend=t0+ds

  end subroutine rk4


!*******************************************************************************
!                                fluxes
! Routine calculates and writes the fluxes at fixed energy and
! pitch angle grids
!*******************************************************************************
      subroutine fluxes(f,hour,xlati1,xlatiS1,xmltS,parmod,HRPee,HRPii,Lstar, &
                        Lstar_max,ompe,CHpower,HIpower,js,ijs,n)
      use constants
      use cgrid
      use WaveGrowth
      use cPlasmasphere,only: rppa,density
      use cfield,only: volume,xmlt,ro,xmlto,y,bo,BriN,BriS,rsb,xmltb,iba
      real y1(ik),fl(ir,ip,je,ig),xlati1(ir,ip),Lstar(ir,ip,0:ik), &
          parmod(10),HRPee(ijs,ir,ip),HRPii(ijs,ir,ip),Lstar_max(0:ik), &
          ompe(ir,ip),CHpower(ir,ip),HIpower(ir,ip), &
          xj(ik),f(ijs,ir,ip,iw,ik),xlatiS1(ir,ip),xmltS(ir,ip)
      integer js(ns)
   
! Calculate fluxes at gridy
      fl=0.
      do j=1,ip
         do i=1,ir     
            do k=1,je  
               k2=k*iw/je
               do m=1,ik
                  y1(m)=y(i,j,m)
                  sumf=f(n,i,j,k2-1,m)*dkeV(n,k2-1)+f(n,i,j,k2,m)*dkeV(n,k2)
                  x=sumf/(dkeV(n,k2-1)+dkeV(n,k2))
                  if (x.le.1.e-30) xj(m)=-30.
                  if (x.gt.1.e-30) xj(m)=log10(x)  ! log flux at (gride, y)
               enddo
               do m=1,ig
                  if (gridy(m).lt.y1(ik)) then
                     x=-30.     ! log flux = -30 near zero degree pitch angle
                  elseif (gridy(m).gt.y1(1)) then
                     x=xj(1)    ! flat dist. near 90 deg pitch angle
                  else
                     call lintp(y1,xj,ik,gridy(m),x)
                  endif
                  fl(i,j,k,m)=10.**x       ! flux at (gride, gridy)     
               enddo
            enddo      ! end of k loop
         enddo         ! end of do i=1,ir    
      enddo            ! end of do j=1,ip

! Set ompe beyond iba
  do j=1,ip
     do i=iba(j)+1,ir
        ompe(i,j)=ompe(iba(j),j)
     enddo
  enddo

! write fluxes
  write(26,'(f10.4,11f9.3," ! hour,parmod,L*max")') hour,parmod,Lstar_max(0)
  do i=1,ir             ! Write fluxes @ fixed E & y grids
     do j=1,ip
        ro1=ro(i,j)
        xmlto1=xmlto(i,j)
        write(26,'(6f8.3,1p,3e10.3,0p,i4)') &
               xlati1(i,j),xmlt(j),xlatiS1(i,j),xmltS(i,j),ro1,xmlto1, &
               BriN(i,j),BriS(i,j),bo(i,j),iba(j)
        write(26,'(1p,9e10.3,0p,2f8.3,1pe10.3)') &
               density(i,j),ompe(i,j),CHpower(i,j),HIpower(i,j), &
               denWP(n,i,j),TparaWP(n,i,j),TperpWP(n,i,j), &
               HRPee(n,i,j),HRPii(n,i,j),rppa(j),Lstar(i,j,0),volume(i,j)
        do k=1,je
           write(26,'(1p,12e11.3)') (fl(i,j,k,m),m=1,ig)
        enddo
     enddo
  enddo

  end subroutine fluxes


!****************************************************************************
  subroutine MapDistribution(f5d)
!****************************************************************************
! Routine maps distribution function at fixed Vperp grid and fits it to 
! ring and Bi-Maxwellian distributions

  use constants
  use cgrid
  use cread2
  use cfield
  use wavegrowth
  implicit none
  integer n,i,j,k,m
  real V1d(iw),f5d(ijs,ir,ip,iw,ik),Vmin,Vmax,dVlog,Vplog,y1d(ik),cosa, &
       xlogf,f2d(iw,ik),Vwp,ywp

  fWP(:,:,:,:,:)=0.
  do n=1,ijs
     if (js(n).lt.5) then       ! only of ions
        do j=1,ip
           do i=1,iba(j)

              ! Find Vper
              y1d(1:ik)=y(i,j,1:ik)
              do k=1,iw
                 Vper(n,i,j,k)=vel(n,k)*y1d(1)
              enddo
   
              ! Find Vpar 
              cosa=sqrt(1.-y1d(ik)*y1d(ik))
              Vmax=vel(n,iw)*cosa
              cosa=sqrt(1.-y1d(1)*y1d(1))   
              do k=1,iw
                 Vmin=vel(n,k)*cosa
                 dVlog=(log10(Vmax)-log10(Vmin))/(ik-1)
                 do m=1,ik
                    Vplog=log10(Vmin)+(m-1)*dVlog
                    Vpar(n,i,j,k,m)=10.**Vplog
                 enddo
              enddo

              ! Find fWP  
              do k=1,iw
                 V1d(k)=log10(vel(n,k))
                 do m=1,ik
                    f2d(k,m)=-30.
                    if (f5d(n,i,j,k,m).gt.1.e-30) f2d(k,m)=log10(f5d(n,i,j,k,m))
                 enddo
              enddo
              do k=1,iw
                 do m=1,ik
                    Vwp=sqrt(Vpar(n,i,j,k,m)**2+Vper(n,i,j,k)**2)
                    ywp=Vper(n,i,j,k)/Vwp
                    Vwp=log10(Vwp)
                    call lintp2(V1d,y1d,f2d,iw,ik,Vwp,ywp,xlogf)
                    fWP(n,i,j,k,m)=10.**xlogf
                 enddo
              enddo

           enddo
        enddo
     endif        ! end of if (js(n).lt.4)
  enddo

  end subroutine MapDistribution


!****************************************************************************
   subroutine echannel(hour,Dst,eout,esum3,bsum,xL1,xL2,n)
!****************************************************************************
!  Routine writes the energy changes by each mechanisum as a
!  function of energy in 3 bands in .ece file.
!  Input: hour,Dst,eout,esum3,n,
!         xL1,xL2 (from cgrid)
!  Output: bsum

   use cgrid
   use cread2
   use closs
   implicit none
   integer n,k,i
   real eout(ns,2,3,je+2),hour,Dst,bsum,xledi,xleei,xleci,xleli,xlebi,&
        esum3(ns,3,ip,je+2),ro4(4),xL1,xL2

  ro4(1:4)=[1.,xL1,xL2,rb]

! Write energy changes in .ece file
  write(3,*) hour,Dst,'       ! hour, Dst'
  write(3,*)' energy loss by : '
  if (js(n).lt.5) write(3,'(a)')'   driftIn     driftOut       drift      '//&
     'ChargeEx     losscone     CoulombCo    changingB      esum'
  if (js(n).eq.5) write(3,'(a)')'   driftIn     driftOut       drift      '//&
     '  wave       losscone     CoulombCo    changingB      esum'
  if (js(n).eq.6) write(3,'(a)')'   driftIn     driftOut       drift      '//&
     'StrongDiff   losscone     CoulombCo    changingB      esum'
  do i=1,3          ! 3 bands
     write(3,*) ro4(i),ro4(i+1),'    ! ro1,ro2'
     do k=1,je+2
        xledi=sum(xled(n,i,1:ip,k))
        xleei=sum(xlee(n,i,1:ip,k))
        xleci=sum(xlec(n,i,1:ip,k))
        xleli=sum(xlel(n,i,1:ip,k))
        xlebi=sum(xleb(n,i,1:ip,k))
        write(3,'(1p,8e13.4)') eout(n,1,i,k),eout(n,2,i,k),xledi,xleei, &
                               xleli,xleci,xlebi,sum(esum3(n,i,1:ip,k))
     enddo
  enddo
  write(3,*) '   '

! Find bsum
  bsum=sum(xleb(n,1:3,1:ip,je+2))

  end subroutine echannel


!***************************************************************************
!                              sume
!  Routine updates esum, psum, and calculates differences in energy and
!  particle number after a process. 
!****************************************************************************
  subroutine sume(xPre,pPre,ijs,isum)
  use cgrid
  use cfield
  use cinitial
  real e0(ir,ip,je+2),dee(ir),dpe(ir),xPre(ns,ir,ip,je+2), &
       gride1(0:je+1),ekev1,pPre(ns,ir,ip,je+2),p0(ir,ip,je+2) 

! Set up gride1(0) and gride1(je+1)
      gride1(0)=0.
      gride1(je+1)=1.e10     ! arbitrary large number

! Calculate esum, psum, etc.
  do n=1,ijs
     e0(1:ir,1:ip,1:je+2)=esum(n,1:ir,1:ip,1:je+2)
     p0(1:ir,1:ip,1:je+2)=psum(n,1:ir,1:ip,1:je+2)
     esum(n,1:ir,1:ip,1:je+2)=0.
     psum(n,1:ir,1:ip,1:je+2)=0.
     gride1(1:je)=ebound(n,1:je)
     do j=1,ip
        xlatb=xlati(iba(j),j)
        do i=1,iba(j)
           do m=1,ik
              do k=1,iw
                 ekev1=ekev(n,k)
                 weight=d4(n)*f2(n,i,j,k,m)
                 weighte=ekev1*weight
                 psum(n,i,j,je+2)=psum(n,i,j,je+2)+weight
                 esum(n,i,j,je+2)=esum(n,i,j,je+2)+weighte
                 kkloop: do kk=1,je+1
                    if (ekev1.gt.gride1(kk-1).and.ekev1.le.gride1(kk)) then
                       psum(n,i,j,kk)=psum(n,i,j,kk)+weight
                       esum(n,i,j,kk)=esum(n,i,j,kk)+weighte
                       exit kkloop
                    endif
                 enddo kkloop
              enddo
           enddo
        enddo               ! end of do i=1,iba(j)

        if (isum.gt.0) then
           do kk=1,je+2
              do i=1,iba(j)
                 xPre(n,i,j,kk)=xPre(n,i,j,kk)+esum(n,i,j,kk)-e0(i,j,kk)
                 pPre(n,i,j,kk)=pPre(n,i,j,kk)+psum(n,i,j,kk)-p0(i,j,kk)
              enddo
           enddo
        endif

     enddo                   ! end of do j=1,ip
  enddo                      ! end of do n=1,ijs

  end subroutine sume

!***************************************************************************
!                              sum3band
!  Routine updates esum3, psum3, and calculates differences in energy and
!  particle number after a process. Everything is binned into 3 bands.
!****************************************************************************
  subroutine sum3band(xL1,xL2,xle,ple)
  use cgrid
  use cread2
  use cinitial
  use cfield,only: ro,iba
  real e0(3,ip,je+2),xle(ns,3,ip,je+2),gride1(0:je+1),ekev1, &
       ple(ns,3,ip,je+2),p0(3,ip,je+2)

! Set up gride1(0) and gride1(je+1)
      gride1(0)=0.
      gride1(je+1)=1.e10     ! arbitrary large number

! Calculate esum3, psum3, etc.
    do n=1,ijs
       rbsum(n)=0.
       rcsum(n)=0.
       e0(1:3,1:ip,1:je+2)=esum3(n,1:3,1:ip,1:je+2)
       p0(1:3,1:ip,1:je+2)=psum3(n,1:3,1:ip,1:je+2)
       esum3(n,1:3,1:ip,1:je+2)=0.
       psum3(n,1:3,1:ip,1:je+2)=0.
       gride1(1:je)=ebound(n,1:je)
       do j=1,ip
          do i=1,iba(j)
             ii=2
             if (ro(i,j).le.xL1) ii=1
             if (ro(i,j).ge.xL2) ii=3
             do m=1,ik
                do k=1,iw
                   ekev1=ekev(n,k)
                   weight=d4(n)*f2(n,i,j,k,m)
                   weighte=ekev1*weight
                   psum3(n,ii,j,je+2)=psum3(n,ii,j,je+2)+weight
                   esum3(n,ii,j,je+2)=esum3(n,ii,j,je+2)+weighte
                   rbsum(n)=rbsum(n)+weighte
                   if (ro(i,j).le.6.6) rcsum(n)=rcsum(n)+weighte
                   kkloop: do kk=1,je+1
                      if (ekev1.gt.gride1(kk-1).and.ekev1.le.gride1(kk)) then
                         psum3(n,ii,j,kk)=psum3(n,ii,j,kk)+weight
                         esum3(n,ii,j,kk)=esum3(n,ii,j,kk)+weighte
                         exit kkloop
                      endif
                   enddo kkloop
                enddo
             enddo
          enddo               ! end of do i=1,iba(j)
          do i=1,3
             do kk=1,je+2
                dee=esum3(n,i,j,kk)-e0(i,j,kk)
                xle(n,i,j,kk)=xle(n,i,j,kk)+dee
                dpe=psum3(n,i,j,kk)-p0(i,j,kk)
                ple(n,i,j,kk)=ple(n,i,j,kk)+dpe
             enddo
          enddo                ! end of i=1,3
       enddo                   ! end of do j=1,ip
    enddo                      ! end of do n=1,ijs

    end subroutine sum3band


!******************************************************************************
!                             lossconeN
! Routine removes losscone particles at north with a lifetime of bounce period
!******************************************************************************
  subroutine lossconeN(f2,ijs)

  use cread1
  use cgrid,only: xk,xkb
  use cfield,only: fcone,xkcN,mcN,iba
  implicit none
  integer m,i,j,k,mc,ijs,n
  real f2(ns,ir,ip,iw,ik),xkc1,dkc,dko,coloss

  do j=1,ip
     do i=1,iba(j)       
        xkc1=xkcN(i,j)
        mc=mcN(i,j)
        if (mc.ge.1.and.mc.lt.ik) then
           dkc=(xkb(mc+1)-xkc1)/(xkb(mc+1)-xkb(mc)) ! fraction of dk in losscone
           dko=1.-dkc                               ! fraction outside losscone
        endif
        do n=1,ijs
           do k=1,iw
              do m=mc+1,ik
                 coloss=fcone(n,i,j,k,m)   ! fraction of losscone loss
                 if (m.eq.mc+1) coloss=dkc*coloss+dko  ! partial loss
                 f2(n,i,j,k,m)=f2(n,i,j,k,m)*coloss
              enddo
           enddo
        enddo
     enddo
  enddo

  end subroutine lossconeN


!******************************************************************************
!                             lossconeS
! Routine removes losscone particles at south with a lifetime of bounce period
!******************************************************************************
  subroutine lossconeS(f2,ijs)

  use cread1
  use cgrid,only: xk,xkb
  use cfield,only: fcone,xkcS,mcS,iba
  implicit none
  integer m,i,j,k,mc,ijs,n
  real f2(ns,ir,ip,iw,ik),xkc1,dkc,dko,coloss

  do j=1,ip
     do i=1,iba(j)
        xkc1=xkcS(i,j)
        mc=mcS(i,j)
        if (mc.ge.1.and.mc.lt.ik) then
           dkc=(xkb(mc+1)-xkc1)/(xkb(mc+1)-xkb(mc)) ! fraction of dk in losscone
           dko=1.-dkc                               ! fraction outside losscone
        endif
        do n=1,ijs
           do k=1,iw
              do m=mc+1,ik
                 coloss=fcone(n,i,j,k,m)   ! fraction of losscone loss
                 if (m.eq.mc+1) coloss=dkc*coloss+dko  ! partial loss
                 f2(n,i,j,k,m)=f2(n,i,j,k,m)*coloss
              enddo
           enddo
        enddo
     enddo
  enddo

  end subroutine lossconeS


!******************************************************************************
!                            Bdrift
!  Routine calculate the change of distribution function due to magnetic drift.
!  Time step is reduced if Courant number is greater 0.5.
!******************************************************************************
  subroutine Bdrift(t,dt,f2,eout,ib0,ijs)
  use cgrid
  use cfield
  use cVdrift
  use cInterFlux
  use cbound,only: fb
  integer ib0(ip)
  real f2(ns,ir,ip,iw,ik),clp(ir,ip), &
       f0(ir,ip),ekev1,gride1(0:je+1),f20(ir,ip),eout(ns,2,3,je+2)

  dtL=dt/dvarL
  dtp=dt/dphi 
  gride1(0)=0.
  gride1(je+1)=1.e10     ! arbitrary large number
  cl(0,:)=0.             ! assume no particle flows across lower L boundary

  do n=1,ijs
  gride1(1:je)=ebound(n,1:je)
  do k=1,iw
  do m=1,ik
     f20(1:ir,1:ip)=f2(n,1:ir,1:ip,k,m)         ! initial f2

     ! find nrun 
     clp(:,:)=0.     ! initial values
     do j=1,ip
        j1=j-1
        if (j1.lt.1) j1=ip
        do i=1,iba(j)
           i1=i-1
           if (i1.lt.1) i1=1
           clp(i,j)=dtL*(abs(vlB(n,i1,j,k,m))+abs(vlB(n,i,j,k,m)))+ &
                    dtp*(abs(vpB(n,i,j1,k,m))+abs(vpB(n,i,j,k,m)))
        enddo
     enddo
     cmax=maxval(clp)
     nrun=ifix(cmax)+1       ! nrun to limit the Courant number 
     dt1L=dtL/nrun
     dt1p=dtp/nrun

     ! Setup boundary fluxes and Courant numbers
     do j=1,ip
        ! boundary conditions
        ib=iba(j)
        ibo=ib0(j)
        fbL0(j)=0.                                    ! f2 at inner L boundary
        fbL1(j)=fb(n,j,k,m)*xjac(n,ir+1,j,k,m)        ! f2 at ir+1         
        if (ibo.lt.ib) f20(ibo+1:ib,j)=0.
        ! Courant numbers
        do i=1,ir   
           cl(i,j)=dt1L*vlB(n,i,j,k,m)
           cp(i,j)=dt1p*vpB(n,i,j,k,m)
        enddo
     enddo

     ! run drift nrun times
     do nn=1,nrun
        if (nrun.lt.5) call FLS_2O_2D(iba,f20)
        if (nrun.ge.5) call FLS_HO_2D(iba,f20)
        f0=f20
        do j=1,ip
           j_1=j-1
           if (j_1.lt.1) j_1=j_1+ip
           ib=iba(j)
           do i=1,ib
              f20(i,j)=f0(i,j)+cl(i-1,j)*faL(i-1,j)-cl(i,j)*faL(i,j)+ &
                                   cp(i,j_1)*fap(i,j_1)-cp(i,j)*fap(i,j) 
                 if (f20(i,j).lt.0.) then      ! use upwind scheme
                    f20(i,j)=f0(i,j)+cl(i-1,j)*fupl(i-1,j)-cl(i,j)*fupl(i,j)+ &
                                     cp(i,j_1)*fupp(i,j_1)-cp(i,j)*fupp(i,j) 
                    if (f20(i,j).lt.0.) then
                       write(*,*)' f2 < 0 in Bdrift, n,i,j,k,m ',n,i,j,k,m
                       write(*,*)' nrun,nn ',nrun,nn
                       write(*,*)'t,f2(i,j,k),f0(i,j,k) ',t,f20(i,j),f0(i,j)
                       write(*,*)'iba(j_1),iba(j) ',iba(j_1),iba(j)
                       write(*,*)'cl(i-1),cl(i) ',cl(i-1,j),cl(i,j)
                       write(*,*)'cp(j_1),cp(j) ',cp(i,j_1),cp(i,j)
                       stop
                    endif
                 endif               ! end of use upwind scheme
              enddo                  ! end of do i=1,ib

              ! Calculate gain or loss at the outer boundary
              ii=2
              if (ro(ib,j).le.xL1) ii=1
              if (ro(ib,j).ge.xL2) ii=3
              j1=j+1
              if (j1.gt.ip) j1=j1-ip
              ekev1=ekev(n,k)
              dPart=-cl(ib,j)*faL(ib,j)*d4(n)
              if (ib.gt.iba(j_1)) dPart=dPart+cp(ib,j_1)*fap(ib,j_1)*d4(n)
              if (ib.gt.iba(j1)) dPart=dPart-cp(ib,j)*fap(ib,j)*d4(n)
              dEner=ekev1*dPart
              if (dPart.gt.0.) then
                 eout(n,1,ii,je+2)=eout(n,1,ii,je+2)+dEner  ! gain
                 kloop1: do kk=1,je+1
                    if (ekev1.gt.gride1(kk-1).and.ekev1.le.gride1(kk)) then
                       eout(n,1,ii,kk)=eout(n,1,ii,kk)+dEner
                       exit kloop1
                    endif
                 enddo kloop1
              else
                 eout(n,2,ii,je+2)=eout(n,2,ii,je+2)+dEner  ! loss
                 kloop2: do kk=1,je+1
                    if (ekev1.gt.gride1(kk-1).and.ekev1.le.gride1(kk)) then
                        eout(n,2,ii,kk)=eout(n,2,ii,kk)+dEner
                        exit kloop2
                    endif
                 enddo kloop2
              endif       ! end of if (dPart.gt.0.)
        enddo          ! end of j loop
     enddo             ! end of do nn=1,nrun

     ! update f2 
     f2(n,1:ir,1:ip,k,m)=f20(1:ir,1:ip)

  enddo            ! end of m loop
  enddo            ! end of k loop
  enddo            ! end of n loop

! Update ib0
  ib0(1:ip)=iba(1:ip)

  end subroutine Bdrift


!******************************************************************************
!                            Edrift
!  Routine calculate the change of distribution function due to electric drift.
!******************************************************************************
  subroutine Edrift(t,dt,f2,eout,ib0,ijs)
  use cgrid
  use cfield
  use cVdrift
  use cInterFlux
  use cbound,only: fb
  use cread2,only: js
  real eout(ns,2,3,je+2),faLk(0:ir,ip,iw),fapk(ir,ip,iw),f20(ir,ip), &
       fww(0:iw+1),cww(0:iw), &
       f2(ns,ir,ip,iw,ik),clp(ir,ip),fawk(ir,ip,0:iw), &
       f0(ir,ip,iw),ekev1,gride1(0:je+1),f2k(ir,ip,iw), &
       fbL1k(ip,iw),fupLk(0:ir,ip,iw),fuppk(ir,ip,iw),cw(ir,ip),faw(0:iw)
  integer ib0(ip)

  dtL=dt/dvarL
  dtp=dt/dphi 
  gride1(0)=0.
  gride1(je+1)=1.e10     ! arbitrary large number
  cl(0,:)=0.             ! assume no particle flows across lower L boundary

  do n=1,ijs
  gride1(1:je)=ebound(n,1:je)
  dtW=dt/dlnp(n)
  pp3=(gridp(n,1)/gridp(n,0))**3
  do m=1,ik
     f2k(1:ir,1:ip,1:iw)=f2(n,1:ir,1:ip,1:iw,m)         ! initial f2

     ! find nrun 
     do i=1,ir
        i1=i-1
        if (i1.lt.1) i1=1
        do j=1,ip
           j1=j-1
           if (j1.lt.1) j1=ip
           clp(i,j)=dtL*(abs(vlE(i1,j))+abs(vlE(i,j)))+ &
                    dtp*(abs(vpE(i,j1))+abs(vpE(i,j)))+ &
                    dtW*2.*abs(vlnp(i,j,m))
        enddo
     enddo
     cmax=maxval(clp)
     nrun=ifix(cmax)+1       ! nrun to limit the Courant number
     dt1L=dtL/nrun
     dt1p=dtp/nrun
     dt1W=dtW/nrun    

     ! Setup boundary fluxes in L and Courant numbers
     do j=1,ip
        ! boundary conditions
        ib=iba(j)
        ibo=ib0(j)
        fbL0(j)=0.                                     ! f2(0)=0
        do k=1,iw
           fbL1k(j,k)=fb(n,j,k,m)*xjac(n,ir+1,j,k,m)   ! f2 at ir+1       
           if (ibo.lt.ib) f2k(ibo+1:ib,j,k)=0.
        enddo
        ! Courant numbers
        do i=1,ir
           cW(i,j)=dt1W*vlnp(i,j,m)
           cl(i,j)=dt1L*vlE(i,j)
           cp(i,j)=dt1p*vpE(i,j)
        enddo
     enddo

     ! run drift nrun times
     do nn=1,nrun
        do k=1,iw
           fbL1(:)=fbL1k(:,k)
           f20(:,:)=f2k(:,:,k)
           call FLS_2O_2D(iba,f20)
           faLk(:,:,k)=faL(:,:)
           fapk(:,:,k)=fap(:,:)
           fupLk(:,:,k)=fupL(:,:)
           fuppk(:,:,k)=fupp(:,:)
        enddo
        do i=1,ir
           do j=1,ip
              cww(0:iw)=cW(i,j)
              fww(1:iw)=f2k(i,j,1:iw)
              if (js(n).lt.5) fww(0)=fww(1)        ! f2(0)=f2(1) for ion
              if (js(n).ge.5) fww(0)=fww(1)/pp3    ! psd(0)=psd(1) for e- 
              fww(iw+1)=fww(iw)                    ! f2(iw+1)=f2(iw)
              call FLS_2O_1D(cww,fww,faw)
              fawk(i,j,0:iw)=faw(0:iw)
           enddo
        enddo
        f0(:,:,:)=f2k(:,:,:)
        do j=1,ip
           j_1=j-1
           if (j_1.lt.1) j_1=j_1+ip
           ib=iba(j)
           do k=1,iw
              do i=1,ib
                 f2k(i,j,k)=f0(i,j,k) &
                             +cl(i-1,j)*faLk(i-1,j,k)-cl(i,j)*faLk(i,j,k) &
                             +cp(i,j_1)*fapk(i,j_1,k)-cp(i,j)*fapk(i,j,k) &
                             +cW(i,j)*fawk(i,j,k-1)-cW(i,j)*fawk(i,j,k)
                 if (f2k(i,j,k).lt.0.) then      ! use upwind scheme
                    f2k(i,j,k)=f0(i,j,k) &
                          +cl(i-1,j)*fupLk(i-1,j,k)-cl(i,j)*fupLk(i,j,k) &
                          +cp(i,j_1)*fuppk(i,j_1,k)-cp(i,j)*fuppk(i,j,k) &
                          +cW(i,j)*fawk(i,j,k-1)-cW(i,j)*fawk(i,j,k)
                    if (f2k(i,j,k).lt.0.) then
                       write(*,*)' f2 < 0 in Edrift, n,i,j,k,m ',n,i,j,k,m
                       write(*,*)' nrun,nn ',nrun,nn
                       write(*,*)'t,f2(i,j,k),f0(i,j,k) ',t,f2k(i,j,k),f0(i,j,k)
                       write(*,*)'iba(j_1),iba(j) ',iba(j_1),iba(j)
                       write(*,*)'cl(i-1),cl(i) ',cl(i-1,j),cl(i,j)
                       write(*,*)'cp(j_1),cp(j) ',cp(i,j_1),cp(i,j)
                       write(*,*)'cW(i,j) ',cW(i,j)
                       stop
                    endif
                 endif               ! end of use upwind scheme
              enddo                  ! end of do i=1,ib

              ! Calculate gain or loss at the outer boundary
              ii=2
              if (ro(ib,j).le.xL1) ii=1
              if (ro(ib,j).ge.xL2) ii=3
              j1=j+1
              if (j1.gt.ip) j1=j1-ip
              ekev1=ekev(n,k)
              dPart=-cl(ib,j)*faLk(ib,j,k)*d4(n)
              if (ib.gt.iba(j_1)) dPart=dPart+cp(ib,j_1)*fapk(ib,j_1,k)*d4(n)
              if (ib.gt.iba(j1)) dPart=dPart-cp(ib,j)*fapk(ib,j,k)*d4(n)
              dEner=ekev1*dPart
              if (dPart.gt.0.) then
                 eout(n,1,ii,je+2)=eout(n,1,ii,je+2)+dEner  ! gain
                 kloop1: do kk=1,je+1
                    if (ekev1.gt.gride1(kk-1).and.ekev1.le.gride1(kk)) then
                       eout(n,1,ii,kk)=eout(n,1,ii,kk)+dEner
                       exit kloop1
                    endif
                 enddo kloop1
              else
                 eout(n,2,ii,je+2)=eout(n,2,ii,je+2)+dEner   ! loss
                 kloop2: do kk=1,je+1
                    if (ekev1.gt.gride1(kk-1).and.ekev1.le.gride1(kk)) then
                        eout(n,2,ii,kk)=eout(n,2,ii,kk)+dEner
                        exit kloop2
                    endif
                 enddo kloop2
              endif       ! end of if (dPart.gt.0.)
           enddo       ! end of k loop
        enddo          ! end of j loop
     enddo             ! end of do nn=1,nrun

     ! update f2
     f2(n,1:ir,1:ip,1:iw,m)=f2k(1:ir,1:ip,1:iw)

  enddo            ! end of m loop
  enddo            ! end of n loop

! Update ib0
  ib0(1:ip)=iba(1:ip)

  end subroutine Edrift


!***********************************************************************
!                         charexchange
!  routine calculates the decay of distributions due to charge exchange
!***********************************************************************
      subroutine charexchange(f2,achar,iba,js,ijs)
      use cimigrid_dim
      real f2(ns,ir,ip,iw,ik),achar(ns,ir,ip,iw,ik)
      integer iba(ip),js(ns)

      do n=1,ijs
         if (js(n).le.4) then        ! ion only
            do j=1,ip
            do i=1,iba(j)
            do m=1,ik
            do k=1,iw
               f2(n,i,j,k,m)=f2(n,i,j,k,m)*achar(n,i,j,k,m) 
            enddo
            enddo
            enddo    
            enddo    
         endif
      enddo    

      end subroutine charexchange


!****************************************************************************
  subroutine Coulomb(f2,coulie,density)
!****************************************************************************
! Routine calculates change of f2 due to Coulomb collisions with the
! plasmasphere

  use cgrid
  use cread2,only: ijs,js,dt
  use cfield,only: iba
  implicit none
  integer n,i,j,k,m
  real f2(ns,ir,ip,iw,ik),coulie(ns,0:iw),fww(0:iw+1),cww(0:iw)
  real faw(0:iw),density(ir,ip),dtw


  do n=1,ijs
     if (js(n).lt.5) then     ! ring current ion only
        dtW=dt/dlnp(n)
        do j=1,ip
           do i=1,iba(j)
              do m=1,ik
                 fww(1:iw)=f2(n,i,j,1:iw,m)      ! f2 before Coulomb collisions
                 ! find Courant numbers
                 cww(:)=dtW*density(i,j)*coulie(n,:)
                 ! setup boundary condition
                 fww(0)=fww(1)                   ! f2(0)=f2(1)
                 fww(iw+1)=fww(iw)               ! f2(iw+1)=f2(iw)
                 call FLS_2O_1D(cww,fww,faw)
                 do k=1,iw
                    f2(n,i,j,k,m)=fww(k)+cww(k-1)*faw(k-1)-cww(k)*faw(k)
                    if (f2(n,i,j,k,m).lt.0.) then
                       write(*,*)'n,i,j,k,m ',n,i,j,k,m
                       write(*,*)'fww,f2 ',fww(k),f2(n,i,j,k,m)
                       write(*,*)'cww,faw, ',cww(k),faw(k)
                       stop
                    endif
                 enddo
              enddo
           enddo
        enddo
     endif     ! end of if (js(n).lt.5)
  enddo

  end subroutine Coulomb


!******************************************************************************
  subroutine energyDep(HRPp,isum)
!******************************************************************************
! Routine calculates accumulated ring current energy deposition to the
! plasmasphere species

  use cimigrid_dim
  use cfield,only : iba
  use cgrid,only : d4,ekev
  use cread2,only : ijs,js
  use cinitial,only : f2,E2D
  implicit none
  integer isum,n,i,j,k,m
  real HRPp(ns,ir,ip),E2D0(ns,ir,ip),E2D1,dH

! Save E2D to E2D0 if isum=1
  if (isum.eq.1) E2D0(:,:,:)=E2D(:,:,:)

! Update E2D
  E2D(:,:,:)=0.
  do j=1,ip
     do i=1,iba(j)
        do n=1,ijs
           if (js(n).lt.5) then     ! ring current ion only
              do m=1,ik
                 do k=1,iw
                    E2D1=ekev(n,k)*f2(n,i,j,k,m)*d4(n)
                    E2D(n,i,j)=E2D(n,i,j)+E2D1
                 enddo
              enddo
           endif
        enddo
     enddo
  enddo

! Calculate HRPp if isum=1
  if (isum.eq.1) then
     do j=1,ip
        do i=1,iba(j)
           do n=1,ijs
              dH=E2D0(n,i,j)-E2D(n,i,j)
              HRPp(n,i,j)=HRPp(n,i,j)+dH    ! HRPp in keV
           enddo
        enddo
     enddo
  endif

  end subroutine energyDep


!****************************************************************************
!                             diffuse_VLF  
!  Routine solves electron diffusion in velocity space in (U,V)=(lnp,lnK) 
!  coordinates due to interactions with chorus and hiss waves
!****************************************************************************
  subroutine diffuse_VLF(f2)

  use constants
  use cWpower
  use cfield,only: xjac,y,tya,ro,bo,iba
  use cgrid,only: xk,dlnk,lnp,dlnp,gridp,ekeV
  use cread2,only: dt,ichor,ihiss,ijs,js
  use cPlasmasphere,only: rppa
  use waveDiffCoef
  implicit none

  integer i,j,m,k,n,L,ipc1,iph1,iLc1,iLh1,iexit,mpa
  real f2(ns,ir,ip,iw,ik),Cfactor(kLat),Upower0,Hfactor,cosa,pp3, &
       Gjac(0:iw+1,0:ik+1),f2d(iw,ik),BLbo
  real ompe1,ro1,rbo,ao,DDm,DDp,y_2,kak,dKda,dtU2,dtV2
  real Daa,Dap,Dpp,DUU(0:iw+1,0:ik+1),DVV(0:iw+1,0:ik+1),DUV(0:iw+1,0:ik+1)
  real,allocatable,dimension(:,:) :: wDaa,wDap,wDpp
  real,allocatable,dimension(:) :: wPA,DD,um,up,a1d,b1d,c1d,f0,fr,fnew

  Upower0=10000.        ! coeff based on UB chorus power of (100pT)^2
! xlam=0.5              ! implicitness in solving diffusion equation
! alam=1.-xlam          ! xlam=1: fully implicit; xlam=0.5: Crank-Nicolson

do n=1,ijs
  if (js(n).eq.5) then       ! high-energy electrons
   dtU2=dt/dlnp(n)/dlnp(n)
   dtV2=dt/dlnK/dlnK
   pp3=(gridp(n,1)/gridp(n,0))**3
   do j=1,ip
     do i=1,iba(j)
        ro1=ro(i,j)
        ompe1=ompe(i,j)                           ! fpe/fce
        Cfactor(:)=CHpowerL(i,j,:)/Cpower0
        Hfactor=HIpower(i,j)/Hpower0

        ! determine whether it is in chorus or hiss region
        iexit=1
        if (ichor.gt.0.and.ro1.gt.rppa(j)) then
           call locate1(xLc,iLc,ro1,iLc1)
           if (iLc1.eq.0.or.iLc1.eq.iLc) goto 9999    ! beyond L range
           if ((xLc(iLc1+1)-ro1).lt.(ro1-xLc(iLc1))) iLc1=iLc1+1
           BLbo=BLc(iLc1)/bo(i,j)
           if (BLbo.gt.1.) Cfactor(:)=Cfactor(:)*BLbo
           call locate1(cOmpe,ipc,ompe1,ipc1)
           if (ipc1.eq.0) goto 9999                  ! beyond fpe/fce range
           if (ipc1.lt.ipc) then
              if ((cOmpe(ipc1+1)-ompe1).lt.(ompe1-cOmpe(ipc1))) ipc1=ipc1+1
           endif
           iexit=0
           mpa=ipa
           allocate (wPA(mpa))
           allocate (wDaa(0:iw+1,mpa),wDap(0:iw+1,mpa),wDpp(0:iw+1,mpa))
           wPA(:)=cPA(:)
           do k=0,iw+1
              wDaa(k,:)=0.
              wDap(k,:)=0.
              wDpp(k,:)=0.
              do L=1,kLat
                 wDaa(k,:)=wDaa(k,:)+cDaa(ipc1,k,:,iLc1,L)*Cfactor(L)
                 wDap(k,:)=wDap(k,:)+cDap(ipc1,k,:,iLc1,L)*Cfactor(L)
                 wDpp(k,:)=wDpp(k,:)+cDpp(ipc1,k,:,iLc1,L)*Cfactor(L)
              enddo
           enddo
        endif
        if (ihiss.gt.0.and.ro1.le.rppa(j)) then
           call locate1(xLh,iLh,ro1,iLh1)
           if (iLh1.eq.0.or.iLh1.eq.iLh) goto 9999       ! beyond L range
           if ((xLh(iLh1+1)-ro1).lt.(ro1-xLh(iLh1))) iLh1=iLh1+1
           BLbo=BLh(iLh1)/bo(i,j)
           if (BLbo.gt.1.) Hfactor=Hfactor*BLbo
           call locate1(hOmpe,iph,ompe1,iph1)
           if (iph1.eq.0) goto 9999                   ! beyond fpe/fce range
           if (iph1.lt.iph) then
              if ((hOmpe(iph1+1)-ompe1).lt.(ompe1-hOmpe(iph1))) iph1=iph1+1
           endif
           iexit=0
           mpa=jpa
           allocate (wPA(mpa))
           allocate (wDaa(0:iw+1,mpa),wDap(0:iw+1,mpa),wDpp(0:iw+1,mpa))
           wPA(:)=hPA(:)
           do k=0,iw+1
              wDaa(k,:)=hDaa(iph1,k,:,iLh1)*Hfactor
              wDap(k,:)=hDap(iph1,k,:,iLh1)*Hfactor
              wDpp(k,:)=hDpp(iph1,k,:,iLh1)*Hfactor
           enddo
        endif
        if (iexit.eq.1) goto 9999

        ! Setup DUU, DVV, DUV
        rbo=ro(i,j)*re_m*sqrt(bo(i,j))
        do m=0,ik+1
           y_2=y(i,j,m)**2
           cosa=sqrt(1.-y_2)
           dKda=2.*rbo*cosa*tya(i,j,m)/y_2
           kak=dKda/xk(m)
           ao=asin(y(i,j,m))*180./pi
           if (ao.lt.wPA(1)) ao=wPA(1)
           if (ao.gt.wPA(mpa)) ao=wPA(mpa)
           do k=0,iw+1
              call lintp(wPA,wDaa(k,:),mpa,ao,Daa)
              call lintp(wPA,wDap(k,:),mpa,ao,Dap)
              call lintp(wPA,wDpp(k,:),mpa,ao,Dpp)
              DUU(k,m)=Dpp
              DUV(k,m)=kak*Dap
              DVV(k,m)=kak*kak*Daa
           enddo
        enddo
        deallocate (wPA,wDaa,wDap,wDpp)

        ! Setup 2D psd in (V,U) and Gjac
        Gjac(:,:)=xjac(n,i,j,:,:) 
        do m=1,ik
           do k=1,iw
              f2d(k,m)=f2(n,i,j,k,m)/xjac(n,i,j,k,m)     ! psd in cimi grid
           enddo
        enddo

        ! do diffusion in U(lnp)
        if (iDpp.eq.1) then
           allocate (a1d(iw),b1d(iw),c1d(iw),fr(iw),fnew(iw),f0(0:iw+1), &
                     DD(0:iw+1),um(iw),up(iw))
           do m=1,ik
              do k=0,iw+1
                 DD(k)=Gjac(k,m)*DUU(k,m)
              enddo
              do k=1,iw
                 DDm=0.5*(DD(k)+DD(k-1))
                 DDp=0.5*(DD(k)+DD(k+1))
                 um(k)=dtU2*DDm/Gjac(k,m)
                 up(k)=dtU2*DDp/Gjac(k,m)
              enddo
              do k=1,iw
                 a1d(k)=-um(k)
                 b1d(k)=1.+(um(k)+up(k))
                 c1d(k)=-up(k)
              enddo
              ! start diffusion in U 
              f0(1:iw)=f2d(1:iw,m)
              f0(0)=f0(1)             ! f0 is psd, psd(0)=psd(1) for e-
              f0(iw+1)=f0(iw)/pp3     ! f0 is psd, f2(iw+1)=f2(iw) for e-
              fr(1:iw)=f0(1:iw)
              fr(1)=fr(1)+um(1)*f0(0)
              fr(iw)=fr(iw)+up(iw)*f0(iw+1)
              call tridiagonal(a1d,b1d,c1d,fr,iw,fnew)
              f2d(:,m)=fnew(:)
           enddo
           deallocate (DD,um,up,a1d,b1d,c1d,f0,fr,fnew)
        endif   ! end of if (iDpp.eq.1)

        ! do diffusion in V(lnK)
        if (iDaa.eq.1) then
           allocate (a1d(ik),b1d(ik),c1d(ik),fr(ik),fnew(ik),f0(0:ik+1), &
                     DD(0:ik+1),um(ik),up(ik))
           do k=1,iw
              do m=0,ik+1
                 DD(m)=Gjac(k,m)*DVV(k,m)
              enddo
              do m=1,ik
                 DDm=0.5*(DD(m)+DD(m-1))
                 DDp=0.5*(DD(m)+DD(m+1))
                 um(m)=dtV2*DDm/Gjac(k,m)
                 up(m)=dtV2*DDp/Gjac(k,m)
              enddo
              do m=1,ik
                 a1d(m)=-um(m)
                 b1d(m)=1.+(um(m)+up(m))
                 c1d(m)=-up(m)
              enddo
              ! start diffusion in V
              f0(1:ik)=f2d(k,1:ik)
              f0(0)=f0(1)        ! flat distribution near 90 deg pitch angle
              f0(ik+1)=0.        ! zero distribution near zero deg pitch angle 
              fr(1:ik)=f0(1:ik)
              fr(1)=fr(1)+um(1)*f0(0)
              fr(ik)=fr(ik)+up(ik)*f0(ik+1)
              call tridiagonal(a1d,b1d,c1d,fr,ik,fnew)
              f2d(k,:)=fnew(:)
           enddo
           deallocate (DD,um,up,a1d,b1d,c1d,f0,fr,fnew)
        endif    ! end of if (iDaa.eq.1)

        ! Update f2
        do m=1,ik
           do k=1,iw
              f2(n,i,j,k,m)=f2d(k,m)*xjac(n,i,j,k,m)
           enddo
        enddo
 9999   continue

     enddo       ! end of i loop
   enddo         ! end of j loop
  endif          ! end of if (js(n).eq.5)
enddo            ! end of n loop (species)

end subroutine diffuse_VLF  


!****************************************************************************
!                             diffuse_UV_ions  
!  Routine solves diffusion in velocity space in (U,V)=(lnp,lnK) coordinates.
!****************************************************************************
  subroutine diffuse_UV_ions(f2)

  use constants
  use cWpower, only:EmicPower0,&
                    EmicHPowerInPs, EmicHePowerInPs, EmicOPowerInPs,&
                    EmicHPowerOutPs,EmicHePowerOutPs,EmicOPowerOutPs,&
                    Ompe 
  use cfield,only: xjac,y,tya,ro,bo,iba
  use cgrid,only: xk,dlnk,lnp,dlnp,gridp,ekeV
  use cread2,only: dt,ichor,ihiss,iEMICdiff,ijs,js
  use cPlasmasphere,only: rppa
  use waveDiffCoef
  use ModCurvScatt, only: iflc,Daa_flc
  implicit none

  integer i,j,m,k,n,L,ipEmic1,iexit,mpa
  real f2(ns,ir,ip,iw,ik),&
       EmicHFact,&   ! H band EMIC wave amplitude factor
       EmicHeFact,&   ! He band EMIC wave amplitude factor
       EmicOFact,&   ! O band EMIC wave amplitude factor
       cosa,pp3, &
       Gjac(0:iw+1,0:ik+1),f2d(iw,ik),BLbo(3)
  real ompe1,ro1,rbo,ao,DDm,DDp,y_2,kak,dKda,dtU2,dtV2
  real Daa,Dap,Dpp,DUU(0:iw+1,0:ik+1),DVV(0:iw+1,0:ik+1),DUV(0:iw+1,0:ik+1),&
       Daaflc(0:iw+1,0:ik+1)
  real,allocatable,dimension(:,:) :: wDaa,wDap,wDpp
  real,allocatable,dimension(:) :: wPA,DD,um,up,a1d,b1d,c1d,f0,fr,fnew

! xlam=0.5              ! implicitness in solving diffusion equation
! alam=1.-xlam          ! xlam=1: fully implicit; xlam=0.5: Crank-Nicolson
mpa=iePA
if (.not.allocated(wPA)) allocate (wPA(mpa))
if (.not.allocated(wDaa)) allocate (wDaa(0:iw+1,mpa))
if (.not.allocated(wDap)) allocate (wDap(0:iw+1,mpa))
if (.not.allocated(wDpp)) allocate (wDpp(0:iw+1,mpa))
wPA(:)=ePA(:)

do n=1,ijs
  if (js(n).lt.4) then       ! ions only
   dtU2=dt/dlnp(n)/dlnp(n)
   dtV2=dt/dlnK/dlnK
   pp3=(gridp(n,1)/gridp(n,0))**3
   do j=1,ip
     if (iEMICdiff==0) goto 999
     do i=1,iba(j)
        ro1=ro(i,j)
        ompe1=ompe(i,j)                           ! fpe/fce
        BLbo(:)=BEmic(:)/bo(i,j)

        iexit=1
        ! outside plasmasphere
        if (ro1.gt.rppa(j)) then
           ! no EMIC wave
           if (EmicHPowerOutPs(i,j)==0..and.&
               EmicHePowerOutPs(i,j)==0..and.&
               EmicOPowerOutPs(i,j)==0.) goto 9999
           call locate1(EmicOmpe,ipEmic,ompe1,ipEmic1)
           if (ipEmic1==0) goto 9999                   ! beyond fpe/fce range
           EmicHFact=EmicHPowerOutPs(i,j)/EmicPower0*BLbo(1)
           EmicHeFact=EmicHePowerOutPs(i,j)/EmicPower0*BLbo(2)
           EmicOFact=EmicOPowerOutPs(i,j)/EmicPower0*BLbo(3)
           if (ipEmic1.lt.ipEmic) then
              if ((EmicOmpe(ipEmic1+1)-Ompe1)<(Ompe1-EmicOmpe(ipEmic1))) &
                 ipEmic1=ipEmic1+1
           endif
           iexit=0
           do k=0,iw+1
              wDaa(k,:)=EmicHDaa(n,ipEmic1,k,:)*EmicHFact&
                       +EmicHeDaa(n,ipEmic1,k,:)*EmicHeFact&
                       +EmicODaa(n,ipEmic1,k,:)*EmicOFact
              wDap(k,:)=EmicHDap(n,ipEmic1,k,:)*EmicHFact&
                       +EmicHeDap(n,ipEmic1,k,:)*EmicHeFact&
                       +EmicODap(n,ipEmic1,k,:)*EmicOFact
              wDpp(k,:)=EmicHDpp(n,ipEmic1,k,:)*EmicHFact&
                       +EmicHeDpp(n,ipEmic1,k,:)*EmicHeFact&
                       +EmicODpp(n,ipEmic1,k,:)*EmicOFact
           enddo
        endif
        ! insdie plasmasphere
        if (ro1.le.rppa(j)) then
           ! no EMIC wave
           if (EmicHPowerInPs(i,j)==0..and.&
               EmicHePowerInPs(i,j)==0..and.&
               EmicOPowerInPs(i,j)==0.) goto 9999
           call locate1(EmicOmpe,ipEmic,ompe1,ipEmic1)
           if (ipEmic1==0) goto 9999                   ! beyond fpe/fce range
           EmicHFact=EmicHPowerInPs(i,j)/EmicPower0*BLbo(1)
           EmicHeFact=EmicHePowerInPs(i,j)/EmicPower0*BLbo(2)
           EmicOFact=EmicOPowerInPs(i,j)/EmicPower0*BLbo(3)
           if (ipEmic1.lt.ipEmic) then
              if ((EmicOmpe(ipEmic1+1)-Ompe1)<(Ompe1-EmicOmpe(ipEmic1))) &
                 ipEmic1=ipEmic1+1
           endif
           iexit=0
           do k=0,iw+1
              wDaa(k,:)=EmicHDaa(n,ipEmic1,k,:)*EmicHFact&
                       +EmicHeDaa(n,ipEmic1,k,:)*EmicHeFact&
                       +EmicODaa(n,ipEmic1,k,:)*EmicOFact
              wDap(k,:)=EmicHDap(n,ipEmic1,k,:)*EmicHFact&
                       +EmicHeDap(n,ipEmic1,k,:)*EmicHeFact&
                       +EmicODap(n,ipEmic1,k,:)*EmicOFact
              wDpp(k,:)=EmicHDpp(n,ipEmic1,k,:)*EmicHFact&
                       +EmicHeDpp(n,ipEmic1,k,:)*EmicHeFact&
                       +EmicODpp(n,ipEmic1,k,:)*EmicOFact
           enddo
        endif

        ! FLC scattering
        if (iflc==1) then
           Daaflc(1:iw,1:ik)=Daa_flc(n,i,j,1:iw,1:ik)
           Daaflc(   0,1:ik)=Daa_flc(n,i,j,1,1:ik)
           Daaflc(iw+1,1:ik)=Daa_flc(n,i,j,iw,1:ik)
           Daaflc(0   ,   0)=Daa_flc(n,i,j,1,1)
           Daaflc(0   ,ik+1)=Daa_flc(n,i,j,1,ik)
           Daaflc(iw+1,   0)=Daa_flc(n,i,j,iw,1)
           Daaflc(iw+1,ik+1)=Daa_flc(n,i,j,iw,ik)
           Daaflc(1:iw,   0)=Daa_flc(n,i,j,1:iw,1)
           Daaflc(1:iw,ik+1)=Daa_flc(n,i,j,1:iw,ik)
           iexit=0
        else if (iflc==0) then
           Daaflc(:,:)=0.
        endif
        if (iexit.eq.1) goto 9999

        ! Setup DUU, DVV, DUV
        rbo=ro(i,j)*re_m*sqrt(bo(i,j))
        do m=0,ik+1
           y_2=y(i,j,m)**2
           cosa=sqrt(1.-y_2)
           dKda=2.*rbo*cosa*tya(i,j,m)/y_2
           kak=dKda/xk(m)
           ao=asin(y(i,j,m))*180./pi
           if (iEMICdiff/=0.or.iflc/=0) then
              if (ao.lt.wPA(1)) ao=wPA(1)
              if (ao.gt.wPA(mpa)) ao=wPA(mpa)
           endif
           do k=0,iw+1
              if (iEMICdiff/=0.or.iflc/=0) then
                 call lintp(wPA,wDaa(k,:),mpa,ao,Daa)
                 call lintp(wPA,wDap(k,:),mpa,ao,Dap)
                 call lintp(wPA,wDpp(k,:),mpa,ao,Dpp)
              else
                 Daa=0.
                 Dap=0.
                 Dpp=0.
              endif
              DUU(k,m)=Dpp
              DUV(k,m)=kak*Dap
              DVV(k,m)=kak*kak*(Daa+Daaflc(k,m))
           enddo
        enddo
        !deallocate (wPA,wDaa,wDap,wDpp)

        ! Setup 2D psd in (V,U) and Gjac
        Gjac(:,:)=xjac(n,i,j,:,:) 
        do m=1,ik
           do k=1,iw
              f2d(k,m)=f2(n,i,j,k,m)/xjac(n,i,j,k,m)     ! psd in cimi grid
           enddo
        enddo

        ! do diffusion in U(lnp)
        if (iDpp.eq.1) then
           allocate (a1d(iw),b1d(iw),c1d(iw),fr(iw),fnew(iw),f0(0:iw+1), &
                     DD(0:iw+1),um(iw),up(iw))
           do m=1,ik
              do k=0,iw+1
                 DD(k)=Gjac(k,m)*DUU(k,m)
              enddo
              do k=1,iw
                 DDm=0.5*(DD(k)+DD(k-1))
                 DDp=0.5*(DD(k)+DD(k+1))
                 um(k)=dtU2*DDm/Gjac(k,m)
                 up(k)=dtU2*DDp/Gjac(k,m)
              enddo
              do k=1,iw
                 a1d(k)=-um(k)
                 b1d(k)=1.+(um(k)+up(k))
                 c1d(k)=-up(k)
              enddo
              ! start diffusion in U 
              f0(1:iw)=f2d(1:iw,m)
              f0(0)=f0(1)             ! f0 is psd, psd(0)=psd(1) for e-
              f0(iw+1)=f0(iw)/pp3     ! f0 is psd, f2(iw+1)=f2(iw) for e-
              fr(1:iw)=f0(1:iw)
              fr(1)=fr(1)+um(1)*f0(0)
              fr(iw)=fr(iw)+up(iw)*f0(iw+1)
              call tridiagonal(a1d,b1d,c1d,fr,iw,fnew)
              f2d(:,m)=fnew(:)
           enddo
           deallocate (DD,um,up,a1d,b1d,c1d,f0,fr,fnew)
        endif   ! end of if (iDpp.eq.1)

        ! do diffusion in V(lnK)
        if (iDaa.eq.1) then
           allocate (a1d(ik),b1d(ik),c1d(ik),fr(ik),fnew(ik),f0(0:ik+1), &
                     DD(0:ik+1),um(ik),up(ik))
           do k=1,iw
              do m=0,ik+1
                 DD(m)=Gjac(k,m)*DVV(k,m)
              enddo
              do m=1,ik
                 DDm=0.5*(DD(m)+DD(m-1))
                 DDp=0.5*(DD(m)+DD(m+1))
                 um(m)=dtV2*DDm/Gjac(k,m)
                 up(m)=dtV2*DDp/Gjac(k,m)
              enddo
              do m=1,ik
                 a1d(m)=-um(m)
                 b1d(m)=1.+(um(m)+up(m))
                 c1d(m)=-up(m)
              enddo
              ! start diffusion in V
              f0(1:ik)=f2d(k,1:ik)
              f0(0)=f0(1)        ! flat distribution near 90 deg pitch angle
              f0(ik+1)=0.        ! zero distribution near zero deg pitch angle 
              fr(1:ik)=f0(1:ik)
              fr(1)=fr(1)+um(1)*f0(0)
              fr(ik)=fr(ik)+up(ik)*f0(ik+1)
              call tridiagonal(a1d,b1d,c1d,fr,ik,fnew)
              f2d(k,:)=fnew(:)
           enddo
           deallocate (DD,um,up,a1d,b1d,c1d,f0,fr,fnew)
        endif    ! end of if (iDaa.eq.1)

        ! Update f2
        do m=1,ik
           do k=1,iw
              f2(n,i,j,k,m)=f2d(k,m)*xjac(n,i,j,k,m)
           enddo
        enddo
 9999   continue

     enddo       ! end of i loop
 999   continue
   enddo         ! end of j loop
  endif          ! end of if (js(n).eq.4)
enddo            ! end of n loop (species)

end subroutine diffuse_UV_ions 


!****************************************************************************
!                             diffuse_flc
!  Routine solves particle diffusion in pitch angle in V(lnK)
!  coordinates due to field line curvature (flc) scattering  
!****************************************************************************
  subroutine diffuse_flc(f2)

  use constants,only: re_m
  use cimigrid_dim,only: ns,ir,ip,iw,ik
  use cfield,only: xjac,y,tya,ro,bo,iba
  use cgrid,only: xk,dlnk
  use cread2,only: dt,ijs
  use ModCurvScatt, only: Daa_flc
  implicit none

  integer i,j,m,k,n
  real f2(ns,ir,ip,iw,ik),cosa,Gjac(0:ik+1),f1d(ik),f0(0:ik+1),fr(ik),fnew(ik)
  real Dsum,wDaa(0:ik+1),rbo,DDm,DDp,y_2,kak,dKda,dtV2
  real DVV(0:ik+1),DD(0:ik+1),um(ik),up(ik),a1d(ik),b1d(ik),c1d(ik)

  do n=1,ijs
     dtV2=dt/dlnK/dlnK
     do j=1,ip
        do i=1,iba(j)
           rbo=ro(i,j)*re_m*sqrt(bo(i,j))
           do k=1,iw
              wDaa(1:ik)=Daa_flc(n,i,j,k,1:ik)
              wDaa(0)=wDaa(1)
              wDaa(ik+1)=wDaa(ik)
              Dsum=sum(wDaa(1:ik))
              if (Dsum.le.0.) goto 9999          ! Daa's are 0, no FLC

              ! Setup DVV
              do m=0,ik+1
                 y_2=y(i,j,m)**2
                 cosa=sqrt(1.-y_2)
                 dKda=2.*rbo*cosa*tya(i,j,m)/y_2
                 kak=dKda/xk(m)
                 DVV(m)=kak*kak*wDaa(m)
              enddo

              ! Setup for diffusion in V(lnK)
              Gjac(:)=xjac(n,i,j,k,:)
              do m=1,ik
                 f1d(m)=f2(n,i,j,k,m)/xjac(n,i,j,k,m)     ! psd in cimi grid
              enddo
              do m=0,ik+1
                 DD(m)=Gjac(m)*DVV(m)
              enddo
              do m=1,ik
                 DDm=0.5*(DD(m)+DD(m-1))
                 DDp=0.5*(DD(m)+DD(m+1))
                 um(m)=dtV2*DDm/Gjac(m)
                 up(m)=dtV2*DDp/Gjac(m)
              enddo
              do m=1,ik
                 a1d(m)=-um(m)
                 b1d(m)=1.+(um(m)+up(m))
                 c1d(m)=-up(m)
              enddo

              ! start diffusion in V
              f0(1:ik)=f1d(1:ik)
              f0(0)=f0(1)        ! flat distribution near 90 deg pitch angle
              f0(ik+1)=0.        ! zero distribution near zero deg pitch angle
              fr(1:ik)=f0(1:ik)
              fr(1)=fr(1)+um(1)*f0(0)
              fr(ik)=fr(ik)+up(ik)*f0(ik+1)
              call tridiagonal(a1d,b1d,c1d,fr,ik,fnew)
              f1d(:)=fnew(:)

              ! Update f2
              do m=1,ik
                 f2(n,i,j,k,m)=f1d(m)*xjac(n,i,j,k,m)
              enddo
      
9999          continue
           enddo     ! end of kloop
        enddo        ! end of iloop
     enddo           ! end of jloop
  enddo              ! end of nloop

end subroutine diffuse_flc


!*****************************************************************************'
!                          EmicPower(t)  
!
! This subroutine calculates the EMIC intensity 
!  base on the CRRES data (Kersten et al. 2014)
!
! EMIC wave amplitude is categorized by Kp index
!
! L = 3.5 - 7, every 0.5
!*****************************************************************************
  subroutine EmicPower(t)
  use cimigrid_dim
  use cfield, only: ro
  use cread2, only: nKp, tKp, zkpa
  use cWpower, only: EmicHPowerInPs ,EmicHePowerInPs ,&
                     EmicHPowerOutPs,EmicHePowerOutPs  
               

  integer,parameter :: nCRRES=14
  real,intent(in) :: t
  real,dimension(nCRRES,3) :: EMIC_CRRES_H,EMIC_CRRES_He
  real,dimension(nCRRES) :: L_CRRES
  integer iKp,jKp,i,j
  
  data EMIC_CRRES_H(:,1) /0.,0.,0.,0.    ,0.    ,0.     ,5.e-5  ,6.25e-5,1.25e-5,0.     ,0.     ,0.,0.,0./
  data EMIC_CRRES_H(:,2) /0.,0.,0.,2.5e-4,4.5e-4,4.75e-4,7.75e-4,1.e-3  ,1.3e-3 ,1.18e-3,4.25e-4,0.,0.,0./
  data EMIC_CRRES_H(:,3) /0.,0.,0.,3.e-3 ,8.5e-3,7.e-3  ,7.e-3  ,2.05e-2,2.5e-2 ,2.e-2  ,1.25e-2,0.,0.,0./
  data EMIC_CRRES_He(:,1) /0.,0.,0.,2.e-4 ,3.e-5 ,2.5e-4 ,2.5e-4,5.e-4,0.001 ,5.e-4 ,5.e-4,0.,0.,0./ 
  data EMIC_CRRES_He(:,1) /0.,0.,0.,3.5e-3,5.e-3 ,0.01   ,0.01  ,0.02 ,2.5e-2,2.5e-2,0.005,0.,0.,0./ 
  data EMIC_CRRES_He(:,1) /0.,0.,0.,5.e-4 ,7.e-2 ,0.1    ,0.1   ,0.1  ,8.e-2 ,3.e-2 ,0.001,0.,0.,0./
  data L_CRRES /1.,2.,3., 3.5 ,4. ,4.5 ,5. ,5.5 ,6. ,6.5 ,7. ,8.,9.,10./


  call locate1(tKp(1:nKp),nKp,t,iKp)

  if      (zkpa(iKp)<2. ) then
     jKp=1
  else if (zkpa(iKp)>=4.) then
     jKp=3
  else!if (zkpa(iKp)>=2..and.zkpa(iKp)<4.) then
     jKp=2
  endif

  ! interpolate EMIC wave intensity
  do j=1,ip
     do i=1,ir
        call lintp(L_CRRES,EMIC_CRRES_H(:,jKp) ,nCRRES,ro(i,j),EmicHPowerInPs(i,j) )
        call lintp(L_CRRES,EMIC_CRRES_He(:,jKp),nCRRES,ro(i,j),EmicHePowerInPs(i,j))
        EmicHpowerOutPs(i,j) =EmicHPowerOutPs(i,j)
        EmicHepowerOutPs(i,j)=EmicHePowerOutPs(i,j)
     enddo
  enddo

  end subroutine EmicPower

!*******************************************************************************
   subroutine FLS_2O_2D(iba,f20)
!  Routine calculates the inter-flux for advection in L and phi. 2nd order
!  scheme is used.
!*******************************************************************************
      use cInterFlux
      integer iba(ip)
      real f20(ir,ip),fwbc(0:ir+1,ip)

      fwbc(1:ir,1:ip)=f20(1:ir,1:ip)   ! fwbc is f20 with B. C.

! Set up boundary condition
      fwbc(0,1:ip)=fbL0(1:ip)
      fwbc(ir+1,1:ip)=fbL1(1:ip)

! find f2o and fup
  fupL(0,:)=0.        ! assume no particle flows across lower L boundary
  faL(0,:)=0.         ! 
  do j=1,ip
     j_1=j-1
     j1=j+1
     j2=j+2
     if (j_1.lt.1) j_1=j_1+ip
     if (j1.gt.ip) j1=j1-ip
     if (j2.gt.ip) j2=j2-ip
     ibm=max(iba(j),iba(j1))
     do i=1,ibm
        ! find f*L
        xsign=sign(1.,cl(i,j))
        fup=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i+1,j)
        flw=0.5*(1.+cl(i,j))*fwbc(i,j)+0.5*(1.-cl(i,j))*fwbc(i+1,j)
        fupL(i,j)=fup
        faL(i,j)=fup
        x=fwbc(i+1,j)-fwbc(i,j)
        if (abs(x).gt.1.e-27.and.i.lt.ir) then
           if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i-1,j))/x
           if (xsign.eq.-1.) r=(fwbc(i+2,j)-fwbc(i+1,j))/x
           if (r.le.0.) faL(i,j)=fup
           if (r.gt.0.) then
              xlimiter=max(min(2.*r,1.),min(r,2.))       ! Superbee
              corr=flw-fup
              faL(i,j)=fup+xlimiter*corr
              if (faL(i,j).lt.0.) faL(i,j)=fup     ! faL can't be < 0
           endif
        endif

        ! find f*p
        xsign=sign(1.,cp(i,j))
        fup=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i,j1)
        flw=0.5*(1.+cp(i,j))*fwbc(i,j)+0.5*(1.-cp(i,j))*fwbc(i,j1)
        fupp(i,j)=fup
        fap(i,j)=fup
        x=fwbc(i,j1)-fwbc(i,j)
        if (abs(x).gt.1.e-27) then
           if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i,j_1))/x
           if (xsign.eq.-1.) r=(fwbc(i,j2)-fwbc(i,j1))/x
           if (r.le.0.) fap(i,j)=fup
           if (r.gt.0.) then
              xlimiter=max(min(2.*r,1.),min(r,2.))
              corr=flw-fup
              fap(i,j)=fup+xlimiter*corr
              if (fap(i,j).lt.0.) fap(i,j)=fup     ! fap can't be < 0
           endif
        endif
     enddo          ! end of iloop
  enddo             ! end of jloop

  end subroutine FLS_2O_2D


!*******************************************************************************
  subroutine FLS_2O_1D(cww,fww,faw)
!*******************************************************************************
! Routine calculates inter-flux, faw in lnp grid using 2nd order flux limited
! scheme with super-bee limiter method

  use cimigrid_dim
  implicit none
  integer k,iup
  real cww(0:iw),fww(0:iw+1),faw(0:iw),xsign,fup,flw,dfw,parR
  real xlimiter,corr

  do k=0,iw
     xsign=sign(1.,cww(k))
     fup=0.5*(1.+xsign)*fww(k)+0.5*(1.-xsign)*fww(k+1)     ! upwind
     flw=0.5*(1.+cww(k))*fww(k)+0.5*(1.-cww(k))*fww(k+1)   ! LW
     dfw=fww(k+1)-fww(k)
     iup=0
     if (k.eq.0.or.k.eq.iw.or.abs(dfw).le.1.e-27) iup=1
     faw(k)=fup
     if (iup.eq.0) then
        if (xsign.eq.1.) parR=(fww(k)-fww(k-1))/dfw
        if (xsign.eq.-1.) parR=(fww(k+2)-fww(k+1))/dfw
        if (parR.le.0.) faw(k)=fup
        if (parR.gt.0.) then
           xlimiter=max(min(2.*parR,1.),min(parR,2.))
           corr=flw-fup
           faw(k)=fup+xlimiter*corr
           if (faw(k).lt.0.) faw(k)=fup     ! faw can't be < 0
        endif
     endif
  enddo

  end subroutine FLS_2O_1D


!*******************************************************************************
  subroutine ReadWavePower
!*******************************************************************************
! Routine reads the chorus and wave power:
! icP=2 - chorus wave power provided by Qianli Ma at Boston University
! icP=3 - chorus wave power provided by Homayon Aryan at UCLA + Gaussian Fit
! icP=4 - chorus wave power provided by Homayon Aryan at UCLA + BU power
! ihP=3 - hiss wave power provided by Homayon Aryan at UCLA

  use cWpower
  use waveDiffCoef,only: kLat
  use cread2,only: ichor,icP,ihP
  implicit none
  integer i,j,k,m,n,AEL,AEU,j0,ii,i0
  real wAMP(kLat),CpowerBU0
  real,allocatable,dimension(:) :: Wpower1 
  real,allocatable,dimension(:,:) :: CpowerBU2D
  real,allocatable,dimension(:,:,:) :: CpowerBU3D
  character header*80

! Read chorus wave power from BU if icP=2 or 4
  if (icP.eq.2.or.icP.eq.4) then
     open(unit=17,file='PowerLBchorus_BU.dat',status='old')
     read(17,*) kAE,kLB,kltB
     allocate (LcBU(kLB),MLTBU(kltB+1))
     allocate (CpowerBU(kAE,kLB,kltB+1,kLat))
     allocate (CpowerBU2D(kLB,kltB),CpowerBU3D(kLB,kltB+1,kLat))
     read(17,'(a80)') header
     do k=1,kAE
        ! Read data
        do i=1,kLB
           do j=1,kltB
              read(17,*) AEL,AEU,LcBU(i),MLTBU(j),wAMP
              do n=1,kLat
                 CpowerBU3D(i,j,n)=wAMP(n)*wAMP(n)
              enddo
           enddo
           CpowerBU3D(i,kltB+1,:)=CpowerBU3D(i,1,:)
        enddo
        ! smooth data
        do n=1,kLat
           ! smooth data in MLT
           do i=1,kLB
              do j=1,kltB
                 j0=j-1
                 if (j0.lt.1) j0=j0+kltB
                 CpowerBU2D(i,j)=(CpowerBU3D(i,j0,n)+CpowerBU3D(i,j,n)+ &
                                  CpowerBU3D(i,j+1,n))/3.
              enddo
           enddo
           ! smooth data in L
           do j=1,kltB
              do i=1,kLB
                 CpowerBU0=0.
                 do ii=-5,5
                    i0=i+ii
                    if (i0.lt.1) i0=1
                    if (i0.gt.kLB) i0=kLB
                    CpowerBU0=CpowerBU0+CpowerBU2D(i0,j)
                 enddo
                 CpowerBU(k,i,j,n)=CpowerBU0/11.
              enddo
           enddo
           CpowerBU(k,:,kltB+1,n)=CpowerBU(k,:,1,n)
        enddo    ! end of do n=1,kLat
     enddo       ! end of do k=1,kAE
     MLTBU(kltB+1)=MLTBU(1)+24.
     close(17)
  endif

! Read chorus wave power from Homayon Aryan if icP.ge.3
  if (icP.ge.3) then
     open(unit=17,file='GSC_ABV.txt',status='old')
     read(17,*) kAE,kBs,kVs,kLA,kltA
     allocate (LcAr(kLA),MLTAr(kltA+1),Wpower1(kltA+1))
     allocate (CpowerAr(kAE,kBs,kVs,kLA,kltA+1,kLat))
     ! setup LcAr and MLTAr
     do i=1,kLA
        LcAr(i)=(i-1)*0.1
     enddo
     do j=1,kltA+1
        MLTAr(j)=(j-1)*1.
     enddo
     ! read wave power
     do k=1,kAE
        do m=1,kBs
           do n=1,kVs
              read(17,'(a80)') header
              do i=1,kLA
                 read(17,*) Wpower1
                 Wpower1(kltA+1)=Wpower1(1)    ! periodic boundary condition
                 CpowerAr(k,m,n,i,1:kltA+1,1)=Wpower1(1:kltA+1)
              enddo
              read(17,'(a80)') header
           enddo
        enddo
     enddo
     close(17)
     deallocate (Wpower1)
  endif

! Read hiss wave power from Homayon Aryan if ihP=3
  if (ihP.eq.3) then
     open(unit=17,file='GSH_ABV.txt',status='old')
     read(17,*) kAE,kBs,kVs,kLA,kltA
     allocate (Lhh(kLA),MLTh(kltA+1),Wpower1(kltA+1))
     allocate (HpowerAr(kAE,kBs,kVs,kLA,kltA+1))
     ! setup Lhh and MLTh
     do i=1,kLA
        Lhh(i)=(i-1)*0.1
     enddo
     do j=1,kltA+1
        MLTh(j)=(j-1)*1.
     enddo
     ! read wave power
     do k=1,kAE
        do m=1,kBs
           do n=1,kVs
              read(17,'(a80)') header
              do i=1,kLA
                 read(17,*) Wpower1
                 Wpower1(kltA+1)=Wpower1(1)    ! periodic boundary condition
                 HpowerAr(k,m,n,i,1:kltA+1)=Wpower1(1:kltA+1)
              enddo
              read(17,'(a80)') header
           enddo
        enddo
     enddo
     close(17)
  endif

  end subroutine ReadWavePower


!*******************************************************************************
  subroutine WavePower(density,Bo,ro,xmlto,iba,icP,ihP,ichor,ihiss)
!*******************************************************************************
! Routine determines Chorus and hiss wave power (pT)^2 and fpe/fce as a function
! of location
!
! Input: density,AE,Bo,ro,xmlto,iba
! Output: CHpowerL,CHpower,HIpower,ompe (through module cWpower)

  use constants
  use convect,only: AE,SWBz,SWVel
  use cWpower
  use waveDiffCoef,only: kLat
  implicit none
  integer iba(ip),iae,jae,iVs,iBz,i,j,n,icP,ihP,ichor,ihiss
  real ro(ir,ip),ro1,xmlt1,Cpower2DA(kLA,kltA+1),Hpower2D(kLA,kltA+1)
  real Cpower2DB(kLB,kltB+1)
  real density(ir,ip),Bo(ir,ip),xmlto(ir,ip),CHpower1

! Determine AE level in 3 bands: 0-100, 100-300, >300
  if (AE.lt.100.) iae=1
  if (AE.ge.100..and.AE.lt.300.) iae=2
  if (AE.ge.300.) iae=3

! Determine AE level in 3 bands: 0-100, 100-500, >500
  if (AE.lt.100.) jae=1
  if (AE.ge.100..and.AE.lt.500.) jae=2
  if (AE.ge.500.) jae=3
 
! Determine Bz level in 3 bands: Bz>-4, -8<Bz<-4, Bz<-8
  if (SWBz.gt.-4.) iBz=1
  if (SWBz.ge.-8..and.SWBz.le.-4.) iBz=2
  if (SWBz.lt.-8.) iBz=3

! Determine Vs level in 3 bands: 0-400, 400-600, >600
  if (SWVel.lt.400.) iVs=1 
  if (SWVel.ge.400..and.SWVel.le.600.) iVs=2 
  if (SWVel.gt.600.) iVs=3 

! Determine ompe (fpe/fce) and wave power (in pT^2)
  CHpowerL(:,:,:)=0.
  do n=1,kLat
     if (icP.eq.2.or.icP.eq.4) Cpower2DB(:,:)=CpowerBU(jae,:,:,n)
     if (icP.ge.3) Cpower2DA(:,:)=CpowerAr(iae,iBz,iVs,:,:,n)
     if (ihP.eq.3.and.n.eq.1) Hpower2D(:,:)=HpowerAr(iae,iBz,iVs,:,:)
     do j=1,ip
        do i=1,iba(j)
           ro1=ro(i,j)
           xmlt1=xmlto(i,j)
           if (xmlt1.lt.0.) xmlt1=xmlt1+24.
           if (n.eq.1) then
              ompe(i,j)=sqrt(density(i,j)*e_mass/epsilon0)/bo(i,j)
              if (ihiss.gt.0) then
                 if (ihP.eq.1) call HissBpower_GN(ro1,xmlt1,jae,HIpower(i,j))
                 if (ihP.gt.1.and.ro1.le.Lhh(kLA)) &
                    call lintp2(Lhh,MLTh,Hpower2D,kLA,kltA+1,ro1,xmlt1, &
                                HIpower(i,j))
              endif
           endif
           if (ichor.gt.0) then
              if (icP.eq.1) call ChorusBpower_GN(ro1,xmlt1,iae,CHpower1)
              if (icP.eq.2) &   
                 call lintp2(LcBU,MLTBU,Cpower2DB,kLB,kltB+1,ro1,xmlt1,CHpower1)
              if (icP.ge.3.and.ro1.le.LcAr(kLA)) &   
                 call lintp2(LcAr,MLTAr,Cpower2DA,kLA,kltA+1,ro1,xmlt1,CHpower1)
              if (icP.eq.3.and.ro1.gt.LcAr(kLA)) &   
                 call ChorusBpower_GN(ro1,xmlt1,iae,CHpower1)
              if (icP.eq.4.and.ro1.gt.LcAr(kLA)) &   
                 call lintp2(LcBU,MLTBU,Cpower2DB,kLB,kltB+1,ro1,xmlt1,CHpower1)
              CHpowerL(i,j,n)=CHpower1
           endif
        enddo
     enddo
  enddo
  CHpower(:,:)=CHpowerL(:,:,1)

  end subroutine WavePower


!******************************************************************************
  subroutine ChorusBpower_GN(ro1,xmlt1,iae,Bpower)
!******************************************************************************
! Routine calculates the AE dependent lower band Chorus power in (pT)^2 by
! Gaussian fits provided by Qiuhua.

  use constants
  implicit none

  integer iae
  real Bpower,ro1,xmlt1,xxbar,yybar,pi2,r12,sigx1,sigy1,qq
  real Amp(3),xbar(3),sigx(3),ybar(3),sigy(3),rxy(3)

  data Amp/18951.,54237.,111104./  
  data xbar/7.46,7.82,5.20/        ! mlt of peak power
  data ybar/6.44,6.37,5.70/        ! ro of peak power
  data sigx/5.50,4.87,4.64/
  data sigy/1.30,1.68,1.77/
  data rxy/0.,-0.1,0./

  pi2=2.*pi        
  r12=1./sqrt(1.-rxy(iae)*rxy(iae))
  sigx1=sigx(iae)
  sigy1=sigy(iae)

  xxbar=abs(xmlt1-xbar(iae))
  if (xxbar.gt.12.) xxbar=24.-xxbar
  yybar=abs(ro1-ybar(iae))
  qq=xxbar**2/sigx1**2+yybar**2/sigy1**2-2.*rxy(iae)*xxbar*yybar/sigx1/sigy1
  Bpower=Amp(iae)*r12/pi2/sigx1/sigy1/exp(0.5*r12*qq)

  end subroutine ChorusBpower_GN


!*****************************************************************************
  subroutine HissBpower_GN(ro1,xmlt1,iae,Bpower)
!*****************************************************************************'
! Routine calculates the CRRES hiss power in (pT)^2 by Gaussian fits
! provided by Qiuhua.

  use constants
  implicit none

  integer iae
  real Bpower,ro1,xmlt1,xxbar,yybar,pi2,r12,sigx1,sigy1,qq
  real Amp(3),xbar(3),sigx(3),ybar(3),sigy(3),rxy(3)

  data Amp/37964.,71459.,130742./          
  data xbar/15.20,11.82,14.46/
  data ybar/3.0,3.25,3.55/
  data sigx/6.08,4.87,5.64/
  data sigy/0.81,1.31,1.32/
  data rxy/0.,0.,0./

  pi2=2.*pi          
  r12=1./sqrt(1.-rxy(iae)*rxy(iae))
  sigx1=sigx(iae)
  sigy1=sigy(iae)

  xxbar=abs(xmlt1-xbar(iae))
  if (xxbar.gt.12.) xxbar=24.-xxbar 
  yybar=abs(ro1-ybar(iae))
  qq=xxbar**2/sigx1**2+yybar**2/sigy1**2-2.*rxy(iae)*xxbar*yybar/sigx1/sigy1
  Bpower=Amp(iae)*r12/pi2/sigx1/sigy1/exp(0.5*r12*qq)

  end subroutine HissBpower_GN


!***********************************************************************
!                          closed
! Routine for numerical integration using closed form 
!
!  S = y1dx1 + y2dx2 + .. + yidxi + .. + yn-1dxn-1 + yndxn
!                 
!  where yi is the value at the middle of dxi
!***********************************************************************
      subroutine closed(n1,n2,y,dx,s)
      real y(n2),dx(n2)

      s=0.
      do i=n1,n2
         s=s+y(i)*dx(i)
      enddo

      end subroutine closed


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
      subroutine indexx(n,arrin,indx)
!--------------------------------------------------------------------------
!  Routine arranges an array in ascending order.

      real arrin(n)
      integer indx(n)
      do 11 j=1,n
        indx(j)=j
11    continue
      l=n/2+1
      irj=n
10    continue
        if(l.gt.1)then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
        else
          indxt=indx(irj)
          q=arrin(indxt)
          indx(irj)=indx(1)
          irj=irj-1
          if(irj.eq.1)then
            indx(1)=indxt
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.irj)then
          if(j.lt.irj)then
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
          endif
          if(q.lt.arrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
          else
            j=irj+1
          endif
        go to 20
        endif
        indx(i)=indxt
      go to 10
  
      end subroutine indexx


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


!-------------------------------------------------------------------------------
!       subroutine lintp2b(x,y,v,nx,ny,x1,y1,v1)
        subroutine lintp2b(y,x,v,ny,nx,y1,x1,v1)
!-------------------------------------------------------------------------------
!  This sub program takes 2-d interpolation. x is 1-D and y is 2-D.
!
  implicit none

  integer nx,ny,i,i1,j,j1,j2,j3
! real y(nx,ny),v(nx,ny)
  real y(ny,nx),v(ny,nx)
  real x(nx),x1,x2,y1,y2,v1,y1d(ny),a,b,b1,q00,q01,q10,q11,minx,maxx,miny,maxy

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

! y1d(1:ny)=y(i,1:ny)
  y1d(1:ny)=y(1:ny,i)
  miny=minval(y1d)
  maxy=maxval(y1d)
  y2=y1
  if (y2.lt.miny) y2=miny      ! force y2 inside the y range
  if (y2.gt.maxy) y2=maxy      !
  call locate1(y1d,ny,y2,j)
  if (j.gt.(ny-1)) j=ny-1        
  if (j.lt.1) j=1                
  j1=j+1       
  b=(y2-y1d(j))/(y1d(j1)-y1d(j))

! y1d(1:ny)=y(i1,1:ny)
  y1d(1:ny)=y(1:ny,i1)
  miny=minval(y1d)
  maxy=maxval(y1d)
  if (y2.lt.miny) y2=miny      ! force y2 inside the y range
  if (y2.gt.maxy) y2=maxy      !
  call locate1(y1d,ny,y2,j2)
  if (j2.gt.(ny-1)) j2=ny-1     
  if (j2.lt.1) j2=1             
  j3=j2+1       
  b1=(y2-y1d(j2))/(y1d(j3)-y1d(j2))

! q00=(1.-a)*(1.-b)
! q01=(1.-a)*b         
! q10=a*(1.-b1)       
! q11=a*b1
! v1=q00*v(i,j)+q01*v(i,j1)+q10*v(i1,j2)+q11*v(i1,j3)
  q00=(1.-b)*(1.-a)
  q01=(1.-b1)*a         
  q10=b*(1.-a)       
  q11=b1*a
  v1=q00*v(j,i)+q01*v(j1,i1)+q10*v(j2,i)+q11*v(j3,i1)

  end subroutine lintp2b


!-----------------------------------------------------------------------------
  subroutine get_tsy_plasma(Bz,Vsw,Nsw,r,phit,t_tsy,n_tsy)
!-----------------------------------------------------------------------------
! Subroutine of Tsyganenko-Mukai plasma sheet model.
! Reference: Tsyganenko, N. A., and T. Mukai, Tail plasma sheet models derived
!               from Geotail particle data, J. Geophys. Res., 108(A3), 1136,
!               doi:10.1029/2002JA009707, 2003.

  use tsy_plasma
  implicit none
  real Bz,Bn,Bs,V,Vsw,r,rho,phit,N,Nsw,t_tsy,n_tsy,temp1,temp2,temp3

  if (Bz.ge.0.) then
     Bn=Bz/5.
     Bs=0.
  else
     Bn=0.
     Bs=-Bz/5.
  endif

  V=Vsw/500.
  rho=r/10.
  N=Nsw/10.

! temperature final fit (equation (4), page 5)
  temp1=a_tsy(1,9)*(V**a_tsy(1,15))+a_tsy(1,10)*Bn+a_tsy(1,11)*Bs
  temp1=temp1*(-1.)*(rho-1.)
  temp2=a_tsy(1,12)*(V**a_tsy(1,16))+a_tsy(1,13)*Bn+a_tsy(1,14)*Bs
  temp2=temp2*(-1.)*(rho-1.)
  temp3=a_tsy(1,5)*V+a_tsy(1,6)*Bn+a_tsy(1,7)*Bs+a_tsy(1,8)*exp(temp2)
  t_tsy=a_tsy(1,1)*V+a_tsy(1,2)*Bn+a_tsy(1,3)*Bs+a_tsy(1,4)*exp(temp1)+ &
        (temp3*sin(phit)*sin(phit))

! density final fit (equation (5), page 6)
  temp1=a_tsy(3,1)+a_tsy(3,2)*(N**a_tsy(3,10))+a_tsy(3,3)*Bn+a_tsy(3,4)*V*Bs
  temp2=a_tsy(3,5)*(N**a_tsy(3,11))+a_tsy(3,6)*Bn+a_tsy(3,7)*V*Bs
  n_tsy=temp1*(rho**a_tsy(3,8)) + temp2*(rho*a_tsy(3,9))*sin(phit)*sin(phit)

  end subroutine get_tsy_plasma  


!-----------------------------------------------------------------------------
  subroutine get_dub_plasma(BsAveN,BsAveT,BnAveT,NswAveN,VswAveT,Rdub,phit, &
                            t_dubyag,n_dubyag)
!-----------------------------------------------------------------------------
! Subroutine of Dubyagin et al. 2016 plasma sheet model.
! Reference: Dubyagin, S., Ganushkina, N. Y., Sillanpää, I., and Runov, A.
!            (2016), Solar wind-driven variations of electron plasma sheet
!            densities and temperatures beyond geostationary orbit during
!            storm times, J. Geophys. Res. Space Physics, 121, 8343– 8360,
!            doi:10.1002/2016JA022947.

  use constants
  use dub_plasma
  implicit none
  real BsAveN,BsAveT,BsAve,BnAveT,BnAve,NswAveN,NswAve,VswAveT,VswAve
  real Rdub,rho,phit,phid,t_dubyag,n_dubyag,temp1,temp2,temp3,temps

  rho=Rdub/10.
  if (phit.le.90.) phid=-phit       ! phit and phid in deg
  if (phit.gt.90..and.phit.lt.270.) phid=phit-180.
  if (phit.ge.270.) phid=360.-phit  
  phid=phid/90.

  NswAve=NswAveN/10.
  BsAve=BsAveN/2.

! density fit (equation (5), page 8353)
  temp1=a_dub(1,1)+a_dub(1,2)*rho+a_dub(1,3)*phid**2*rho
  temp2=a_dub(1,4)*phid**2+a_dub(1,5)*NswAve
  temp3=(a_dub(1,6)+a_dub(1,7)*rho)*BsAve
  n_dubyag=temp1+temp2+temp3
  if (n_dubyag.lt.0.15) n_dubyag=0.15   ! make minimum density ge 0.15

  VswAve=VswAveT/400.
  BsAve=BsAveT/2.
  BnAve=BnAveT/2.

! temperature fit (equation (8), page 8354)
  temp1=a_dub(2,1)+a_dub(2,2)*phid+a_dub(2,3)*VswAve
  temp2=(a_dub(2,4)+a_dub(2,5)*phid**2*rho)*BsAve**a_dub(2,7)
  temp3=a_dub(2,6)*rho*BnAve**a_dub(2,8)
  temps=temp1+temp2+temp3
  if (temps.lt.0.37) temps=0.37   ! make min temperature ge 0.01 keV
  t_dubyag=temps**a_dub(2,9)

  end subroutine get_dub_plasma


!--------------------------------------------------------------------------
      subroutine geocorona(igeo,t,xsm,ysm,zsm,hden)
!--------------------------------------------------------------------------
! Routine calculates geocorona density. 
!
! igeo=1: Chamberlain model of [H] fitted by an exponential function. The fit 
!         matches Rairden etal [1986] for radial distance from 1.08 to 12 Re.
! igeo=2: Hodges model
!
! input: igeo,xsm,ysm,zsm
! output: hden (H density in m^-3)
   
  use hodges_data
  use ModMate, only: mate_geocorona
  implicit none
  integer igeo
  real t,xsm,ysm,zsm,hden,radius,rexob,rr,ar,ch0,ch1,ch2

  radius=sqrt(xsm*xsm+ysm*ysm+zsm*zsm)

! Raiden's model of [H] when igeo=1
  ch0=10.692
  ch1=-4.4431
  ch2=0.715831
  if (igeo.eq.1) then  
     rexob=1.08   ! 1.08 = exobase in Rairden et al. [1986]
     rr=radius
     if (rr.lt.1.) rr=1.
     ar=log(rr/rexob)       ! rexob = exobase in Rairden et al. [1986]
     if (rr.lt.rexob) then
        ar=abs(ar)
        ch1=1.
     endif
     hden=exp(ch0+ch1*ar**ch2)  
  endif

! Hodges model of [H] when igeo=2
  if (igeo.eq.2) call hodges_geocorona(xsm,ysm,zsm,hden)

! MATE [H] density in /cm3
  if (igeo.eq.3) call mate_geocorona(t,xsm,ysm,zsm,hden)

! Convert hden in m^-3
  hden=hden*1.e6        ! H density in m^-3
            
  end subroutine geocorona


!-----------------------------------------------------------------------------
  subroutine hodges_init
!-----------------------------------------------------------------------------
! GEOCORONA MODEL: Hodges, 1994
! Hodges geocorona model is valid between 1.04 and 9.74 RE

  use hodges_data
  implicit none
  integer :: i
  character*100 :: StringHeader

  Alm=-1000.
  Blm=-1000.
  Ylm=-1000.

! Get hodges_file and defines if winter or summer for Solstice data files
  read(4,'(a)') hodges_file
  read(4,'(a6)') season
  close(4)

  open(300,FILE='Hodges_data/'//trim(hodges_file),status='old')
  read (300,*) StringHeader
  read (300,*) StringHeader
  do i = 1, nr
     read(300,*) tabular_rad(i),tabular_n(i),Alm(1,0,i),Alm(1,1,i),Blm(1,1,i), &
                 Alm(2,0,i),Alm(2,1,i),Blm(2,1,i),Alm(2,2,i),Blm(2,2,i), &
                 Alm(3,0,i),Alm(3,1,i),Blm(3,1,i),Alm(3,2,i),Blm(3,2,i), &
                 Alm(3,3,i), Blm(3,3,i)
  end do
  close(300)

  Alm(0,0,1:nr)=1.e+4
  Blm(0,0,1:nr)=0.
  Blm(1,0,1:nr)=0.
  Blm(2,0,1:nr)=0.
  Blm(3,0,1:nr)=0.

  do i=1, nr-1
     rad_deriv(i)=tabular_rad(i+1)-tabular_rad(i)
  enddo

  Alm = Alm * 1.e-4
  Blm = Blm * 1.e-4

  end subroutine hodges_init


!-----------------------------------------------------------------------------
  subroutine hodges_geocorona(x,y,z,hden)
!-----------------------------------------------------------------------------
! The subroutine calculates Hodges geocorona model from R. Hodges
! Monte Carlo simulation of the terrestrial hydrogen exosphere, JGR, 1994.
! Coordinate system - idealized when magnetic pole coincide with geographic pole
! Hodges defines 'phi_h' as longitude measured from the midnight meridian, and
! theta is colatitude measured from the summer pole.
!
! There are 8 .dat files for equinox/solstice and for different F107 index
!
! Input: x,y,z in SM
! Output hden, H density is in cm^-3

  use constants
  use hodges_data
  implicit none
  integer :: i
  real :: x,y,z,hden,rad,radxy,cos_phi,cos_theta,phi_h,theta,rad2,den2,den1, &
          rad_E

  rad = sqrt(x*x+y*y+z*z)
  radxy=  sqrt(x*x+y*y)
  rad_E=6371.
  rad2=rad*rad_E
  if (rad2.lt.tabular_rad(1)) rad2=tabular_rad(1)
  if (rad2.gt.tabular_rad(nr)) rad2=tabular_rad(nr)

!!! angle calculations
  cos_phi = -x / radxy
  cos_theta = z / rad
  if (season .eq. 'winter') cos_theta = -z / rad
  phi_h=acos(cos_phi)
  theta=acos(cos_theta)
  if (y .lt. 0) phi_h = 2.*pi - phi_h

  do i=1,nr-1
     if ( (rad2 .ge. tabular_rad(i))  .and. (rad2 .le. tabular_rad(i+1)) ) then
        call hodges_den_tabular(i,Theta,Phi_h,den1)
        call hodges_den_tabular(i+1,Theta,Phi_h,den2)
        if ((den1 .lt. 0) .or. (den2 .lt. 0) ) then
           write(*,*) 'ERROR in Hodges model: negative density'
           stop
        endif
        ! linear interpolation in r:
        hden = den1 + (rad2- tabular_rad(i))*(den2-den1)/rad_deriv(i)
     endif
  enddo
      
  end subroutine hodges_geocorona


!-----------------------------------------------------------------------------
  subroutine hodges_den_tabular(jR,Theta,Phi_h,Hdensity)
!-----------------------------------------------------------------------------
! jr is radial distance INDEX, from Hodges tables

  use constants
  use hodges_data
  implicit none
  integer :: l,m,jR
  real :: Theta, Phi_h, Hdensity, Z

  Z = 0.0
  do l = 0, 3
     do m = 0, l
        call hodges_Ylm(Theta)
        Z=Z+(Alm(l,m,jR)*cos(m*Phi_h)+Blm(l,m,jR)*sin(m*Phi_h))*Ylm(l,m)
     end do
  end do
  Hdensity = tabular_n(jR) * sqrt(4. * Pi) * Z
      
  end subroutine hodges_den_tabular


!-----------------------------------------------------------------------------
  subroutine hodges_Ylm(Theta)
!-----------------------------------------------------------------------------

  use constants
  use hodges_data
  implicit none
  real, intent(in)  :: Theta
  real ::  cos2Theta, sin2Theta, cos3Theta, sin3Theta

  cos2Theta = (cos(Theta))**2
  sin2Theta = (sin(Theta))**2
  cos3Theta = (cos(Theta))**3
  sin3Theta = (sin(Theta))**3
  Ylm(0,0) =   sqrt(1./(4.*Pi))
  Ylm(1,0) =   sqrt(3./(4.*Pi)) * cos(Theta)
  Ylm(1,1) = - sqrt(3./(8.*Pi)) * sin(Theta)
  Ylm(2,0) =   sqrt(5./(4.*Pi)) * (3./2. * cos2Theta - 0.5)
  Ylm(2,1) = - sqrt(15./(8.*Pi)) * sin(Theta) * cos(Theta)
  Ylm(2,2) = 1./4. * sqrt(15./(2.*Pi)) * sin2Theta
  Ylm(3,0) = sqrt(7./(4.*Pi)) * (5./2. * cos3Theta - 3./2. * cos(Theta))
  Ylm(3,1) = - 1./4. * sqrt(21./(4.*Pi)) * sin(Theta) * (5. * cos2Theta -1.)
  Ylm(3,2) =   1./4. * sqrt(105./(2.*Pi)) * sin2Theta * cos(Theta)
  Ylm(3,3) = - 1./4. * sqrt(35./(4.*Pi)) * sin3Theta
      
end subroutine hodges_Ylm


!-----------------------------------------------------------------------------
      function gammln(xx)
!-----------------------------------------------------------------------------
!    ln(gamma(xx))

      real cof(6)
      data cof,stp1/76.18009173,-86.50532033,24.01409822, &
          -1.231739516,0.120858003e-2,-0.536382e-5,2.50662827465/
      data half,one,fpf/0.5,1.0,5.5/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      gammln=tmp+log(stp1*ser)

      end function gammln


!-----------------------------------------------------------------------------
  function erf(xx)
!-----------------------------------------------------------------------------
! Error function

  implicit none
  integer i
  real erf,xx,xsq,tt,tn,pp,aa(5),ssf

  pp=0.3275911
  aa=[0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429]
  tt=1./(1.+pp*xx)
  xsq=xx*xx

  tn=1.
  ssf=0.
  do i=1,5
     tn=tn*tt
     ssf=ssf+aa(i)*tn
  enddo
  erf=1.-ssf/exp(xsq)

  end function erf


!-----------------------------------------------------------------------------
  function gfun(xx)
!-----------------------------------------------------------------------------
! G function in the Coulomb drag coefficient

  use constants
  implicit none
  real gfun,xx,g1
  
  g1=erf(xx)-2.*xx/sqrt(pi)*exp(-xx*xx)
  gfun=g1/2./xx/xx
      
  end function gfun


!-----------------------------------------------------------------------------
  subroutine traceF(imod,intB,xi,yi,zi,dir,dsmax,err,rlim,rmn,parmod,psi, &
                    np,xf,yf,zf,xa,ya,za,ra,ba,npf,iout)
!-----------------------------------------------------------------------------
! Routine does field line tracing from (xi,yi,zi) to rmn
!
! xi,yi,zi,xf,yf,zf,xa,ya,za are in RE and in sm coordinates
! ba is in nT

  implicit none
  external dip_08,IGRF_GSW_08,t96_01,t04_s,zeroB
  integer,parameter :: i_one=1,m_one=-1
  integer np,np1,npf,iout,m,imod,intB
  real xi,yi,zi,dir,rlim,rmn,parmod(10),psi,xf,yf,zf,rf,xa(np),ya(np),za(np), &
       xg,yg,zg,ra(np),ba(np),xa1(np),ya1(np),za1(np),bxint,byint,bzint, &
       bxext,byext,bzext,bx,by,bz,dsmax,err,rmn1
   
! Initial setup 
  np1=np-1   
  iout=0 
  call smgsw_08(xi,yi,zi,xg,yg,zg,i_one)      ! sm to gsm
  
! Start fieldline tracing
     if (imod.eq.0.and.intB.eq.0) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
                rmn,imod,parmod,zeroB,dip_08,xf,yf,zf,xa1,ya1,za1,npf,np)
     if (imod.eq.1.and.intB.eq.0) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
                rmn,imod,parmod,t96_01,dip_08,xf,yf,zf,xa1,ya1,za1,npf,np)
     if (imod.eq.2.and.intB.eq.0) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
                rmn,imod,parmod,t04_s,dip_08,xf,yf,zf,xa1,ya1,za1,npf,np)
     if (imod.eq.0.and.intB.eq.1) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
                rmn,imod,parmod,zeroB,IGRF_GSW_08,xf,yf,zf,xa1,ya1,za1,npf,np)
     if (imod.eq.1.and.intB.eq.1) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
                rmn,imod,parmod,t96_01,IGRF_GSW_08,xf,yf,zf,xa1,ya1,za1,npf,np)
     if (imod.eq.2.and.intB.eq.1) call trace_08(xg,yg,zg,dir,dsmax,err,rlim, &
                rmn,imod,parmod,t04_s,IGRF_GSW_08,xf,yf,zf,xa1,ya1,za1,npf,np)

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
        if (imod.eq.1) call t96_01(imod,parmod,psi,xa1(m),ya1(m),za1(m), &
                                   bxext,byext,bzext)
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

!  Limit the values of parmod(1:4) in t96 model (imod=1)
      if (imod.eq.1) then
         if (parmod(1).gt.10.) parmod(1)=10.             
         if (parmod(2).lt.-100.) parmod(2)=-100.
         if (parmod(2).gt.20.) parmod(2)=20.
         if (parmod(3).lt.-10.) parmod(3)=-10.
         if (parmod(3).gt.10.) parmod(3)=10.
         if (parmod(4).lt.-10.) parmod(4)=-10.
         if (parmod(4).gt.10.) parmod(4)=10.
      endif

!  parmod(5:10) for t04_s: w04(1:6) defined in Tsyganenko and Sitnov, 2005
      if (imod.eq.2) then
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
      endif        ! end of if (imod.eq.2) 

! Set limit to parmod for t04_s
      if (imod.eq.2) then
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
      endif

      end subroutine TsyParmod


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


! ************************************************************************
!                          tridiagonal
!
!  Subroutine solves inversion of tridiagonal matrix
!   through Thomas's tridiagonal matrix algorithm.
!
!  b1*fi-1 + b2*fi + b3*fi+1 = b4,  i= 1,2,3,..,n
!
!  Inputs : b1,b2,b3,b4,n
!  Outputs : f
! ************************************************************************
     subroutine tridiagonal(b1,b2,b3,b4,n,f)
     implicit none

     integer,intent(in) :: n
     real,intent(in) :: b1(n),b2(n),b3(n),b4(n)
     real,intent(out) :: f(n)
     real c0,c1(n),c2(n)
     integer i

     c1(1)=b3(1)/b2(1)
     c2(1)=b4(1)/b2(1)
     do i=2,n
        c0=b2(i)-b1(i)*c1(i-1)
        c1(i)=b3(i)/c0
        c2(i)=(b4(i)-b1(i)*c2(i-1))/c0
     enddo

     f(n)=c2(n)
     do i=n-1,1,-1
        f(i)=c2(i)-c1(i)*f(i+1)
     enddo

     end subroutine tridiagonal


!-------------------------------------------------------------------------------
  subroutine setgrid(nlp,npp,varLp,mphip)
!-------------------------------------------------------------------------------
! Routine set the grid for the plasmasphere calculation
! Input: varL,mphi,nlp,npp
! Output: varLp,mphip,Brip,xlatp,xlatpS

  use constants
  use cgrid, only: ir,ip,varL,mphi,xlati,ksai
  use cPlasmasphere_new, only: Brip,xlatp,xlatpS,ksaip
  use cfield, only: BriN,xlatiS
  implicit none

  integer nlp,npp,i,j
  real varL1(ir),varLp(nlp),mphip(npp),dvarLp,dphip
  real mphi1(ip+1),BriN1(ir,ip+1),xlati1(ir,ip+1),xlatiS1(ir,ip+1)
  real ksai1(ir,ip+1)

! Setup the varLp grid
  dvarLp=(varL(ir)-varL(1))/(nlp-1)
  do i=1,nlp
     varLp(i)=varL(1)+(i-1)*dvarLp
  enddo

! Setup the mphip grid
  dphip=2.*pi/npp
  do j=1,npp
     mphip(j)=(j-1)*dphip          ! magnetic longitude, 0 to 2pi   
  enddo

! Add point ip+1
  mphi1(1:ip)=mphi(1:ip)
  mphi1(ip+1)=mphi1(1)+2.*pi
  BriN1(1:ir,1:ip)=BriN(1:ir,1:ip)
  BriN1(1:ir,ip+1)=BriN1(1:ir,1)
  xlati1(1:ir,1:ip)=xlati(1:ir,1:ip)
  xlati1(1:ir,ip+1)=xlati1(1:ir,1)
  xlatiS1(1:ir,1:ip)=xlatiS(1:ir,1:ip)
  xlatiS1(1:ir,ip+1)=xlatiS1(1:ir,1)
  ksai1(1:ir,1:ip)=ksai(1:ir,1:ip)
  ksai1(1:ir,ip+1)=ksai1(1:ir,1)

! Find ionospheric parameter in fine grid: Brip, xlatp, xlatpS
  varL1(1:ir)=varL(1:ir)
  do j=1,npp
    do i=1,nlp   
       call lintp2(varL1,mphi1,BriN1,ir,ip+1,varLp(i),mphip(j),Brip(i,j))
       call lintp2(varL1,mphi1,xlati1,ir,ip+1,varLp(i),mphip(j),xlatp(i,j))
       call lintp2(varL1,mphi1,xlatiS1,ir,ip+1,varLp(i),mphip(j),xlatpS(i,j))
       call lintp2(varL1,mphi1,ksai1,ir,ip+1,varLp(i),mphip(j),ksaip(i,j))
    enddo
  enddo

  end subroutine setgrid


!-------------------------------------------------------------------------------
  subroutine mapgrid
!-------------------------------------------------------------------------------
! Routine finds local time phip, phipS, flux tube volume volp,
! equatorial rp, ibp at the plasmasphere grid.
! Input: ir,ip,volume,ro
! Output: rp,volp,ibp

  use constants, only: pi
  use cfield, only: ir,ip,volume,ro,phi
  use cgrid, only: varL,mphi
  use cPlasmasphere_new
  implicit none

  integer i,j,ibp1,imap
  real varL1(ir),pi2,mphi1(ip+1),ro1(ir,ip+1),volume1(ir,ip+1),deltp

  pi2=2.*pi

! Find phip and phipS
  deltp=phi(1)-mphi(1)
  do j=1,npp
     phip(j)=mphip(j)+deltp
     if (phip(j).gt.pi2) phip(j)=phip(j)-pi2
  enddo
  phipS(1:npp)=phip(1:npp)   ! assume south foot pts have same MLT as north

! Setup CIMI parameters for mapping
  varL1(1:ir)=varL(1:ir)
  do j=1,ip
     mphi1(j)=mphi(j)
  enddo
  mphi1(ip+1)=mphi1(1)+pi2
  do j=1,ip
     ro1(1:ir,j)=ro(1:ir,j)
     volume1(1:ir,j)=volume(1:ir,j)
  enddo
  ro1(1:ir,ip+1)=ro1(1:ir,1)
  volume1(1:ir,ip+1)=volume1(1:ir,1)

! Find volp, rp and ibp if imap=1
  do j=1,npp
     do i=1,nlp   
        call lintp2(varL1,mphi1,volume1,ir,ip+1,varLp(i),mphip(j),volp(i,j))
        call lintp2(varL1,mphi1,ro1,ir,ip+1,varLp(i),mphip(j),rp(i,j))
        if (i.gt.1.and.rp(i,j).gt.rp(i-1,j)) ibp1=i
     enddo
     ibp(j)=ibp1
  enddo

  end subroutine mapgrid   


!-------------------------------------------------------------------------------
  subroutine mapPotent
!-------------------------------------------------------------------------------
! Routine maps convection velocity into the plasmasphere grid.
! Input: ir,ip,potent 
! Output: vlEp,vpEp

  use constants, only: pi
  use cgrid, only: ir,ip,varL,mphi
  use cVdrift, only: vlE,vpE
  use cPlasmasphere_new
  implicit none

  integer i,j
  real varL1(ir),mphi1(ip+1),vlE1(ir,ip+1),vpE1(ir,ip+1)

! Setup CIMI parameters for mapping
  varL1(1:ir)=varL(1:ir)
  do j=1,ip
     mphi1(j)=mphi(j)
  enddo
  mphi1(ip+1)=mphi1(1)+pi*2.
  do j=1,ip
     vlE1(1:ir,j)=vlE(1:ir,j)
     vpE1(1:ir,j)=vpE(1:ir,j)
  enddo
  vlE1(1:ir,ip+1)=vlE1(1:ir,1)
  vpE1(1:ir,ip+1)=vpE1(1:ir,1)

! Find vlEp, vpEp
  do j=1,npp
     do i=1,nlp
        call lintp2(varL1,mphi1,vlE1,ir,ip+1,varLp(i),mphip(j),vlEp(i,j))
        call lintp2(varL1,mphi1,vpE1,ir,ip+1,varLp(i),mphip(j),vpEp(i,j))
     enddo
  enddo

  end subroutine mapPotent

!-------------------------------------------------------------------------------
  subroutine get_density(nlp,npp,Nion,volp,varLp,mphip,ir,ip,varL,mphi)
!-------------------------------------------------------------------------------
! Routine extracts density from high resolution plasmasphere density and 
! locate the plasmapause, rppa
! Input: nlp,npp,Nion,volp,ir,ip
! Output: density (in m^-3),rppa

  use constants,only: pi
  use cPlasmasphere,only: density,densityP,rppa
  use cfield,only: iba,ro
  implicit none

  integer nlp,npp,ir,ip,i,j,ipp(ip),ipp1
  real varLp(nlp),mphip(npp),Nion(nlp,npp),volp(nlp,npp),denp(nlp,npp+1), &
       varL(ir),mphi(ip),phip(npp),mphip1(npp+1)

! Get density in plasmasphere grid
  do i=1,nlp
     do j=1,npp
        denp(i,j)=Nion(i,j)/volp(i,j)
     enddo
     denp(i,npp+1)=denp(i,1)
  enddo
  mphip1(1:npp)=mphip(1:npp)
  mphip1(npp+1)=mphip1(1)+2.*pi

! Get plasmasphere density (m^-3) in CIMI grid
  do j=1,ip 
  do i=1,ir     
     call lintp2(varLp,mphip1,denp,nlp,npp+1,varL(i),mphi(j),density(i,j))
  enddo
  enddo

! Find plasmapause location as a function of mlt by densityP
  do j=1,ip
     ipp(j)=iba(j)-1    ! initial value
     findPP1: do i=3,iba(j)-1
        if (densityP.gt.density(i,j)) then
           ipp(j)=i-1
           exit findPP1 
        endif
     enddo findPP1 
     ipp1=ipp(j)+1
     rppa(j)=0.5*(ro(ipp(j),j)+ro(ipp1,j))
  enddo

  end subroutine get_density


!-------------------------------------------------------------------------------
  subroutine EMICGrowthRate(hour,js,ijs)
!-------------------------------------------------------------------------------
! Routine calculates EMIC growth rates with real ion distribution.

  use constants
  use cfield
  use WaveGrowth
  use cPlasmasphere,only: Tcold,Rcold,density
  use cgrid
  implicit none
  integer i,j,k,m,mg,n,ijs,js(ns),nj
  real Wpg(ir,ip),empc,Wwpp(ir,ip),emp,delta(ir,ip),Vgr(nw,ir,ip),chmass, &
       hour,aMass,cpw,wwc,DD,Dpar,Rpar(nw,ir,ip),kpar(nw,ir,ip), &
       V21,Vrat,dVp(ijs,ir,ip,iw),Ireal(nw,ir,ip),pWpg, &
       denWH(ir,ip),denCH(ir,ip),par0,par1,par2,par3,par4,par5,par6,epar2, &
       gammaReal(nw,ir,ip),denjw,cosa,Rhc(ir,ip),Nhot, &
       TparaWH(ir,ip),Iedp(nw,ir,ip),gammaEDp(nw,ir,ip),Vth,denCP

! Find the ratio of hot and cold ion density
  Rhc(:,:)=0.
  do i=1,ir
     do j=1,ip    
        Nhot=0.
        do n=1,ijs
           if (js(n).lt.5) Nhot=Nhot+denWP(n,i,j)    ! sum ions only
        enddo
        if (density(i,j).gt.0.) Rhc(i,j)=Nhot/density(i,j)
     enddo
  enddo

! Find the warm proton density, denWH, and cold density, denCH
  do n=1,ijs
     if (js(n).eq.1) then
        denCH(1:ir,1:ip)=density(1:ir,1:ip)*Rcold(1)
        denWH(1:ir,1:ip)=denWP(n,1:ir,1:ip)
        TparaWH(1:ir,1:ip)=TparaWP(n,1:ir,1:ip)
     endif
  enddo

! Calculate proton gyrofrequency, Wpg, warm proton plasma frequency, Wwpp, and
! delta, ratio of cold and warm plasma density
  empc=echarge/xmp
  emp=echarge/sqrt(xmp*epsilon0)
  do j=1,ip
     do i=1,iba(j)
        if (denWP(1,i,j).le.0.) then
           write(*,*) 'denWP(1,i,j) = 0, i,j => ',i,j
           stop
        endif
        Wpg(i,j)=empc*bo(i,j)
        Wwpp(i,j)=emp*sqrt(denWH(i,j))
        delta(i,j)=denCH(i,j)/denWH(i,j)
     enddo
     do i=iba(j)+1,ir
        Wpg(i,j)=Wpg(iba(j),j)
     enddo
  enddo

! Calculate group velocity, Vgr, R parameter, Rpar, and k parameter, kpar
  do j=1,ip
     do i=1,iba(j)
        cpw=2.*EM_speed*Wpg(i,j)/Wwpp(i,j)
        wwc=Wwpp(i,j)/EM_speed
        do m=1,nw
           DD=(1.+delta(i,j))/(1.-Xpar(m))
           Rpar(m,i,j)=DD*(2.-Xpar(m))/(1.-Xpar(m))
           do nj=2,3    ! sum over heavy ion, O+(nj=2) and He+(nj=3)
              denjw=0.
              do n=1,ijs
                 if (js(n).eq.nj) denjw=denWP(n,i,j)
              enddo
              par0=1.-Xpar(m)*xmass1(nj)
              if (par0.ne.0.) then
                 par1=(denjw+density(i,j)*Rcold(nj))/denWH(i,j)
                 par2=xmass1(nj)/par0
                 par3=par2*(2.-Xpar(m)*xmass1(nj))/par0
                 DD=DD+par1*par2
                 Rpar(m,i,j)=Rpar(m,i,j)+par1*par3
              endif
           enddo
           Dpar=0.
           if (DD.gt.0.) Dpar=sqrt(DD)
           kpar(m,i,j)=wwc*Xpar(m)*Dpar
           Vgr(m,i,j)=cpw*Dpar/Rpar(m,i,j)        ! Vgr in m/s
        enddo
     enddo
     do i=iba(j)+1,ir
        Vgr(1:nw,i,j)=Vgr(1:nw,iba(j),j)
     enddo
  enddo

! Calculate dVp for warm ions
  do n=1,ijs
     do j=1,ip
        do i=1,iba(j)
           v21=Vper(n,i,j,2)/Vper(n,i,j,1)
           Vrat=(V21-1.)/sqrt(V21)
           do k=1,iw
              dVp(n,i,j,k)=Vrat*Vper(n,i,j,k)
           enddo
        enddo
     enddo
  enddo

! Calculate I function for EMIC damping
  Iedp(:,:,:)=0.
  do n=1,3              ! 3 ion species (H+, O+, He+) in the plasmasphere
     aMass=xmass1(n)
     chmass=1.6e-19/aMass/xmp
     Vth=sqrt(2.*Tcold*chmass)
     do j=1,ip
        do i=1,iba(j)
           denCP=density(i,j)*Rcold(n)
           par0=denCP/aMass/aMass/Vth/pi**1.5
           do m=1,nw
              if (kpar(m,i,j).gt.0.) then
                 par1=Wpg(i,j)*(Xpar(m)*aMass-1.)/aMass/kpar(m,i,j)/Vth
                 par2=par1*par1
                 if (par2.lt.80.) then
                    epar2=exp(par2)
                    Iedp(m,i,j)=Iedp(m,i,j)-par0*Xpar(m)*aMass/epar2
                 endif
              endif
           enddo
        enddo
     enddo
  enddo

! Calculate numeric I function with real distribution, Ireal
  call Ifunction(kpar,Wpg,dVp,Ireal)

! Calculate growth rate with real distribution
  gammaReal(:,:,:)=0.
  gammaEDp(:,:,:)=0.
  par3=sqrt(pi)/2./EM_speed/EM_speed
  do j=1,ip
     do i=1,iba(j)
        pWpg=pi*Wpg(i,j)
        par4=par3*Wpg(i,j)*denWH(i,j)*Wwpp(i,j)*Wwpp(i,j)/density(i,j)
        par5=pWpg*pWpg/denWH(i,j)
        do m=1,nw
           if (kpar(m,i,j).gt.0.) then
              par6=par5/kpar(m,i,j)/Xpar(m)/Rpar(m,i,j)
              gammaReal(m,i,j)=par6*Ireal(m,i,j)
              gammaEDp(m,i,j)=par6*Iedp(m,i,j)
           endif
        enddo
     enddo
  enddo

! write results
  write(9,*) hour,'    ! hour'
  write(9,*) ' ro ---'
  write(9,'(11f8.3)') ro
  write(9,*) ' xmlto ---'
  write(9,'(11f8.3)') xmlto
  write(9,*) ' ratio of hot and cold density ---'
  write(9,'(11f8.3)') Rhc   
  write(9,*) ' Hot H+ Tpara (keV) ---'
  write(9,'(1p,8e11.3)') TparaWH
  write(9,*) ' H+ gyrofrequency ---'
  write(9,'(1p,8e11.3)') Wpg
  write(9,*) ' EMIC wave group velcoity (m/s) ---'
  write(9,'(1p,8e11.3)') Vgr
  write(9,*) ' EMIC wave number: k (m^-1) ---'
  write(9,'(1p,8e11.3)') kpar
  write(9,*) ' EMIC wave growth rate with real distribution ---'
  write(9,'(1p,8e11.3)') gammaReal
  write(9,*) ' EMIC damping rate ---'
  write(9,'(1p,8e11.3)') gammaEDp  

  end subroutine EMICGrowthRate


!-------------------------------------------------------------------------------
  subroutine Ifunction(kpar,Wpg,dVp,Ifunc)
!-------------------------------------------------------------------------------
! Routine calculates the I function in EMIC growth with distribution fWP 
! Input: kpar,Wpg,dVp,fWP 
! Output: Ifunc

  use constants
  use cgrid
  use cread2
  use cfield
  use WaveGrowth
  implicit none
  integer n,i,j,k,m,mg,k0,k1
  real Ifunc(nw,ir,ip),kpar(nw,ir,ip),aMass,kMWp,Wpg(ir,ip),Ur,Ur1,Usign, &
       cosa,Vpar1(ik+1),f1d(ik+1), &
       fUr(iw),dfdVper,dfdVpar,par3,par4,dVp(ijs,ir,ip,iw)

  Ifunc(:,:,:)=0.
  do j=1,ip
     do i=1,iba(j)
        Xloop: do m=1,nw
           do n=1,ijs
              if (kpar(m,i,j).gt.0..and.js(n).lt.5) then    ! Sum over ions
                 aMass=xmass(n)/xmp
                 kMWp=kpar(m,i,j)*aMass/Wpg(i,j)
                 Ur1=(Xpar(m)*aMass-1.)/kMWp            ! resonant velocity
                 Ur=abs(Ur1)
                 Usign=1.
                 if (Ur.gt.0.) Usign=Ur1/Ur
                 ! Find fUr, distribution at Vpar=Ur
                 do k=1,iw
                    do mg=1,ik
                       Vpar1(mg+1)=Vpar(n,i,j,k,mg)
                       f1d(mg+1)=fWP(n,i,j,k,mg)
                    enddo
                    Vpar1(1)=-Vpar1(2)  ! add one point of negative Vparallel
                    f1d(1)=f1d(2)       !
                    fUr(k)=0.
                    if (Ur.le.Vpar1(ik+1)) call lintp(Vpar1,f1d,ik+1,Ur,fUr(k))
                 enddo
                 ! integrate over Vper to get Ifunc
                 do k=1,iw
                   ! find dfdVper
                   k0=max(k-1,1)
                   k1=min(k+1,iw)
                   dfdVper=(fUr(k1)-fUr(k0))/(Vper(n,i,j,k1)-Vper(n,i,j,k0))
                   ! find dfdVpar
                   Vpar1(2:ik+1)=Vpar(n,i,j,k,1:ik)
                   f1d(2:ik+1)=fWP(n,i,j,k,1:ik)
                   Vpar1(1)=-Vpar1(2)  ! add one point of negative Vparallel
                   f1d(1)=f1d(2)       !
                   call locate1(Vpar1,ik+1,Ur,mg)
                   if (mg.gt.ik) mg=ik
                   dfdVpar=Usign*(f1d(mg+1)-f1d(mg))/(Vpar1(mg+1)-Vpar1(mg))
                   par3=(Vper(n,i,j,k)/aMass)**2
                   par4=(dfdVper+kMWp*Vper(n,i,j,k)*dfdVpar)*par3
                   Ifunc(m,i,j)=Ifunc(m,i,j)+par4*dVp(n,i,j,k)
                 enddo
              endif       ! end if (kpar(m,i,j).gt.0..and.js(n).lt.5)
           enddo
        enddo Xloop
     enddo
  enddo

  end subroutine Ifunction


!-------------------------------------------------------------------------------
  subroutine derivative_3pt(x0,x1,x2,f0,f1,f2,xj,dfdx)
!-------------------------------------------------------------------------------
! Routine calculates derivative df/dx at xj using 3-point formula
! (R. Burden and J. Faires, Numerical Analysis, Prindle, Weber & Schmidt, 1985)

  implicit none
  real x0,x1,x2,f0,f1,f2,xj,dfdx,der0,der1,der2

  der0=f0*(2.*xj-x1-x2)/((x0-x1)*(x0-x2))
  der1=f1*(2.*xj-x0-x2)/((x1-x0)*(x1-x2))
  der2=f2*(2.*xj-x0-x1)/((x2-x0)*(x2-x1))
  dfdx=der0+der1+der2

  end subroutine derivative_3pt


!-------------------------------------------------------------------------------
  subroutine hybridsetup
!-------------------------------------------------------------------------------
! Routine setup for using output from Yu Lin's hybrid model

  use hybrid
  implicit none
  integer i
  real,pointer,dimension(:,:) :: temp2D
  character*40 ccc,hyfile

  open(unit=20,file='hybrid.dat')
  read(20,'(a40)') hyfile

  open(unit=21,file=trim(hyfile))
  do i=1,11
     read(21,*) ccc   ! to skip i=1-11
  enddo
 
! Read plasma data at equator. Unit: N #/cm^3, V km/s, T-eV
  read(21,*) nxp,nyp        ! at equator
  write(*,*) 'nxp,nyp ',nxp,nyp
  allocate(xpp(nxp),ypp(nyp))
  allocate(temp2D(nxp,nyp),eqdata(nxp,nyp,6))   ! eqdata : (N,Vx,Vy,Vz,T||,Tper)
  read(21,'(5e12.4)') xpp
  read(21,'(5e12.4)') ypp
  do i=1,6
     read(21,'(5e12.4)') temp2D
     eqdata(:,:,i)=temp2D
  enddo

! Read 3-D Magnetic field in Cartisian coordinate. Unit: B-nT
  read(21,*) nxb,nyb,nzb    ! 3-D data
  write(*,*) 'nxb,nyb,nzb ',nxb,nyb,nzb
  allocate(xbb(nxb),ybb(nyb),zbb(nzb))
  allocate(bxx(nxb,nyb,nzb),byy(nxb,nyb,nzb),bzz(nxb,nyb,nzb))
  read(21,'(5e12.4)') xbb
  read(21,'(5e12.4)') ybb
  read(21,'(5e12.4)') zbb
  read(21,'(5e12.4)') bxx
  read(21,'(5e12.4)') byy
  read(21,'(5e12.4)') bzz

! Read electric potential at ionosphere. Unit: potenthy(volt)
  read(21,*) mih,nih
  write(*,*) 'mih,nih ',mih,nih
  allocate(thetah(mih),phih(nih),potenthy(mih,nih))
  read(21,'(5e12.4)') thetah
  read(21,'(5e12.4)') phih
  read(21,'(5e12.4)') potenthy

! write(*,*) 'mih,thetah(1,2,mih) ',mih,thetah(1),thetah(2),thetah(mih)
! write(*,*) 'nih,phih(1,2,nih) ',nih,phih(1),phih(2),phih(nih)

  close(20) 
  close(21) 

  end subroutine hybridsetup


!*******************************************************************************
!                                StDiTime
!  Routine calculate the strong diffusion lifetime for sub-keV electrons.
!*******************************************************************************
  subroutine StDiTime(ijs,js,dt,volume,xlati,rc,re_m,xme,iba)

  use cStDiTime
  use cgrid,only : vel
  implicit none

  integer ijs,js(ns),iba(ip),i,j,k,m,n
  real volume(ir,ip),xlati(ir,ip),dt,rc,re_m,xme,eb,xmer3,sinlat2,Bi,vBe,SDtime1

  eb=0.25                         ! fraction of back scatter e-
  xmer3=xme/(rc*re_m)**3

  SDtime(:,:,:,:)=0.
  do n=1,ijs
     if (js(n).eq.6) then        ! sub-keV electrons
        do j=1,ip
           do i=1,iba(j)
              sinlat2=sin(xlati(i,j))*sin(xlati(i,j))
              Bi=xmer3*sqrt(3.*sinlat2+1.)      ! magnetic field at ionosphere
              vBe=2.*volume(i,j)*Bi/(1.-eb)
              do k=1,iw
                 do m=1,ik
                    SDtime1=vBe/vel(n,k) ! strong diff T, (gamma*mo/p=1/v)
                    SDtime(i,j,k,m)=exp(-dt/SDtime1)
                 enddo
              enddo
           enddo
        enddo
     endif
  enddo

  end subroutine StDiTime


!***********************************************************************
!                            StrongDiff
!  Routine calculate the change of electron psd (f2) by strong diffusion
!***********************************************************************
  subroutine StrongDiff(f2,iba,js,ijs)

  use cStDiTime
  implicit none

  integer ijs,iba(ip),js(ns),i,j,k,m,n
  real f2(ns,ir,ip,iw,ik)

      do n=1,ijs
         if (js(n).eq.6) then        ! low-energy electrons only
            do j=1,ip
            do i=1,iba(j)
            do m=1,ik
            do k=1,iw
               f2(n,i,j,k,m)=f2(n,i,j,k,m)*SDtime(i,j,k,m)
            enddo
            enddo
            enddo
            enddo
         endif
      enddo

  end subroutine StrongDiff


!***********************************************************************
  subroutine conductance_model_v2(fac0,facm1,facp1,imlt,sigmap,sigmah)
!***********************************************************************
! Conductance model provided by Bob Robinson on 4 December 2016.
! Input the fac at any magnetic latitude (fac0). pos. -> upward current
! the fac one degree equatorward (facm1)
! the fac one degree poleward (facp1)
! and the magnetic local time as an integer from 0 to 24 (imlt)
! Output:
!  sigmap=Pedersen conductance
!  hoverp=Hall to Pdersen ratio
!  sigmah=Hall conductance
!  enflux=precipitating particle energy flux in ergs/cm2-s

  implicit none

  integer imlt
  real spm0(0:24),spm1(0:24),spp0(0:24),spp1(0:24),hpm(0:24),hpp(0:24),pedm1, &
       pedp1,fac0,facm1,facp1,sigmap,hoverp,sigmah,pedm,pedp,avge,enflux

!  set the linear fit coefficients for the conductance model:
 spm0=[1.90641, 2.25115, 2.57007, 2.85271, 3.09056, 3.27706, 3.40765, 3.47969, &
       3.49254, 3.44751, 3.34787, 3.19886, 3.00767, 2.78347, 2.53740, 2.28253, &
       2.03393, 1.80862, 1.62557, 1.50574, 1.47203, 1.54932, 1.76444, 2.14620, &
       1.90641]
 spm1=[-2.48804, -2.43332, -2.57739, -2.85514, -3.20905, -3.58924, -3.95343, &
       -4.26700, -4.50291, -4.64178, -4.67182, -4.58889, -4.39646, -4.10561, &
       -3.73508, -3.31119, -2.86790, -2.44681, -2.09712, -1.87565, -1.84686, &
       -2.08283, -2.66325, -3.67544, -2.48804]
 spp0=[0.321099, 0.447197, 0.517054, 0.542299, 0.533535, 0.500332, 0.451231, &
       0.393743, 0.334348, 0.278496, 0.230609, 0.194077, 0.171259, 0.163485, &
       0.171057, 0.193242, 0.228282, 0.273386, 0.324733, 0.377474, 0.425726, &
       0.462579, 0.480095, 0.469298, 0.321099]
 spp1=[13.0462, 17.5440, 20.3049, 21.6573, 21.9016, 21.3101, 20.1273, 18.5696, &
       16.8254, 15.0550, 13.3910, 11.9378, 10.7718, 9.94146, 9.46723, 9.34159, &
       9.52900, 9.96593, 10.5609, 11.1943, 11.7188, 11.9588, 11.7108, 10.7433,&
       13.0462]
 hpm=[1.00137, 1.01065, 1.01222, 1.00684, 0.995311, 0.978486, 0.957258, & 
      0.932567,0.905402, 0.876798, 0.847838, 0.819651, 0.793412, 0.770343, &
      0.751715,0.738842, 0.733089, 0.735864, 0.748624, 0.772872, 0.810158, &
      0.862080,0.930280,1.01645, 0.0905402]
 hpp=[0.944843, 1.04363, 1.11876, 1.17208, 1.20554, 1.22109, 1.22077, 1.20662, &
     1.18077, 1.14538, 1.10265, 1.05483, 1.00424, 0.953226, 0.904178, 0.859546,&
     0.821822, 0.793546, 0.777307, 0.775739, 0.791526, 0.827396, 0.886130, &
     0.970551, 0.944843]

    pedm1=0.
    pedp1=0.
!  if fac is changing sign, use the average values on either side
    if (fac0*facm1.lt.0.) then 
       if (facm1.lt.0.) pedm1=spm0(imlt)+spm1(imlt)*facm1
       if (facm1.ge.0.) pedm1=spp0(imlt)+spp1(imlt)*facm1
       if (facp1.lt.0.) pedp1=spm0(imlt)+spm1(imlt)*facp1
       if (facp1.ge.0.) pedp1=spp0(imlt)+spp1(imlt)*facp1
       sigmap=0.5*(pedm1+pedp1)
       hoverp=0.5*(hpm(imlt)+hpp(imlt)) ! Hall cond. by average of ratios for 
       sigmah=sigmap*hoverp             ! +ve and -ve fac
    else 
       if (fac0.lt.-0.1) then 
          sigmap=spm0(imlt)+spm1(imlt)*fac0
          sigmah=sigmap*hpm(imlt)
          hoverp=hpm(imlt)
       endif
       if (fac0.gt.0.1) then 
          sigmap=spp0(imlt)+spp1(imlt)*fac0
          sigmah=sigmap*hpp(imlt)
          hoverp=hpp(imlt)
       endif
       if (fac0.ge.-0.1.and.fac0.le.0.1) then 
          pedm=spm0(imlt)+spm1(imlt)*fac0
          pedp=spp0(imlt)+spp1(imlt)*fac0
          sigmap=0.5*(pedm+pedp)
          hoverp=0.5*(hpm(imlt)+hpp(imlt))
          sigmah=sigmap*hoverp
       endif
    endif
    avge=((1/.45)*hoverp)**(1./0.85)
    if (avge.gt.0) enflux=(sigmap/(40.*avge/(16.+avge*avge)))**2
    if (avge.le.0) enflux=0.

  end subroutine conductance_model_v2


!*******************************************************************************
!                            FLS_HO_2D
!
!  Routine calculates the inter-flux, f(i+cl,j) and f(i,j+cp),
!  using up to 7th order Lagrangian polynomial interpolation
!  and ULTIMATE universal limiter as in Leonard [1991].
!  where cl=cl(i,j), cp=cp(i,j)
!
!  Note: Use odd-order scheme.
!        odd-order scheme is more accurate then even order scheme.
!
!  fhol,fhop,flw: higher order inter flux
!  fupl,fupp,fup: 1st order inter flux
!  a00,b00 : modified Pascal's progression
!
!  By S.-B. Kang, Code 673, at NASA/GSFC, August 2017
!
!  Modification history
!  Jan 10, 2020 by Mei-Ching Fok
!   - Change fwbc(-2:0,j)=f(1,j) to fwbc(0,j)=fbL0(j)
!   - In populate fwbc, change do i=iba(j)+1,ir+4 to fwbc(ir+1,j)=fbL1(j)
!   - At i=1,2,ir-1 and ir, lower-order schemes are used.
!*******************************************************************************

      subroutine FLS_HO_2D(iba,f)
      use cimigrid_dim, only: ir,ip
      use cVdrift, only: ihol,ihop
      use cInterFlux, only: fbL0,fbL1,cl,cp,fhol=>faL,fhop=>fap,fupL,fupp
      implicit none

      integer,intent(in) :: iba(ip)
      real,intent(in) ::  f(ir,ip)
      real fwbc(0:ir+1,ip),df,cl1,cp1,cl2,cp2,c1,&
           xsignl(ir,ip),xsignp(ir,ip),xsign,fup,flw2,flw3,flw4,flw5,flw6,flw7,&
           x,r,xlimiter,corr,d1,d2,fc,fu,fd,ref,del,adel,acurv,flim
      real a21,a22,a23,a24,&
           a41,a42,a43,a44,a45,a46,&
           a61,a62,a63,a64,a65,a66,a67,a68,&        ! a_(order)_(ordinal)
           b31,b32,b33,b34,&
           b51,b52,b53,b54,b55,b56,&
           b71,b72,b73,b74,b75,b76,b77,b78          ! b_(order)_(ordinal)
      integer i0,i1,ib,ib2,ibm
      integer i01,i02,i03,i10,i20,i30,i40,&! i_(positive)_(negative)
              j01,j02,j03,j10,j20,j30,j40  ! j_(positive)_(negative)
      integer i,j,k         ! dummy index

      i0=1

! Set aOO,bOO
      a21=1./12.
      a22=-1./12.
      a23=a22
      a24=a21
      a41=1./240.
      a42=-3./240.
      a43=2./240.
      a44=a43
      a45=a42
      a46=a41
      a61=1./10080.
      a62=-5./10080.
      a63=9./10080.
      a64=-5./10080.
      a65=a64
      a66=a63
      a67=a62
      a68=a61

      b31=1./24.
      b32=-3./24.
      b33=3./24.
      b34=-1./24.
      b51=1./720.
      b52=-5./720.
      b53=10./720.
      b54=-10./720.
      b55=5./720.
      b56=-1./720.
      b71=1./10080.
      b72=-7./10080.
      b73=21./10080.
      b74=-35./10080.
      b75=35./10080.
      b76=-21./10080.
      b77=7./10080.
      b78=-1./10080.

      fwbc(1:ir,1:ip)=f(1:ir,1:ip)   ! flwbc is f with boundary condition

! Set up boundary condition
      do j=1,ip
         fwbc(0,j)=fbL0(j)     
         fwbc(ir+1,j)=fbL1(j)
      enddo

! Calculate xsign
      do i=i0,ir
         do j=1,ip
            xsignl(i,j)=sign(1.,cl(i,j))
            xsignp(i,j)=sign(1.,cp(i,j))
         enddo
      enddo

! find 1st order fup and higher order inter-flux
      fupl(0,:)=0.        ! assume no particle flows across lower L boundary
      fhol(0,:)=0.        ! 
      do j=1,ip
         j01=j-1
         j02=j-2
         j03=j-3
         j10=j+1
         j20=j+2
         j30=j+3
         j40=j+4
         if (j01.lt.1) j01=j01+ip
         if (j02.lt.1) j02=j02+ip
         if (j03.lt.1) j03=j03+ip
         if (j10.gt.ip) j10=j10-ip
         if (j20.gt.ip) j20=j20-ip
         if (j30.gt.ip) j30=j30-ip
         if (j40.gt.ip) j40=j40-ip
         ib=max(iba(j),iba(j10))
         do i=1,ib
            cl1=cl(i,j)
            cp1=cp(i,j)
            cl2=cl1**2
            cp2=cp1**2
            i01=i-1
            i02=i-2
            i03=i-3
            i10=i+1
            i20=i+2
            i30=i+3
            i40=i+4

       ! find f*l
       ! 1st order
            xsign=xsignl(i,j)
            if (xsign.eq.1.) then
               fup=fwbc(i,j)
            else
               fup=fwbc(i10,j)
            endif
            fupl(i,j)=fup
  !    ! 2nd order (Lax-Wendroff)
  !         flw2=0.5*((fwbc(i10,j)+fwbc(i,j))+cl1*(fwbc(i,j)-fwbc(i10,j)))   
  !         if (ihol.eq.2.or.i20.gt.ir+1) then
  !            ! Superbee limter
  !            x=fwbc(i10,j)-fwbc(i,j)
  !            if (abs(x).le.1.e-27) fhol(i,j)=fup
  !            if (abs(x).gt.1.e-27.and.i.lt.ir) then
  !               if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i01,j))/x
  !               if (xsign.eq.-1.) r=(fwbc(i20,j)-fwbc(i10,j))/x
  !               if (r.le.0.) fhol(i,j)=fup
  !               if (r.gt.0.) then
  !                  xlimiter=max(min(2.*r,1.),min(r,2.))
  !                  corr=flw2-fup
  !                  fhol(i,j)=fup+xlimiter*corr
  !                  if (fhol(i,j).lt.0.) fhol(i,j)=fup     ! fsbp can't be < 0
  !               endif
  !            endif
  !            goto 7777        ! end of finding fhol
  !         endif
  !    ! 3rd order (QUICKEST)
  !         if (ihol.ge.3) then
  !            d1=a21*fwbc(i20,j)+a22*fwbc(i10,j)+a23*fwbc(i,j)+a24*fwbc(i01,j) 
  !            d2=b31*fwbc(i20,j)+b32*fwbc(i10,j)+b33*fwbc(i,j)+b34*fwbc(i01,j)
  !            c1=cl2-1.
  !            if (ihol.eq.3.or.i02.lt.0.or.i30.gt.ir+1) then
  !               d2=d2*2.
  !               flw3=flw2+c1*(d1-xsign*d2)
  !               fhol(i,j)=flw3
  !               goto 7777        ! end of finding fhol
  !            endif
  !         endif
  !    ! 4th order
  !         if (ihol.ge.4) then
  !            flw4=flw2+c1*(d1-cl1*d2)   ! 4th order inter-flux
  !            if (ihol.eq.4) then
  !               fhol(i,j)=flw4
  !               goto 7777        ! end of finding fhol
  !            endif
  !         endif
  !    ! 5th order
  !         if (ihol.ge.5) then
  !            d1=a41*fwbc(i30,j)+a42*fwbc(i20,j)+a43*fwbc(i10,j)+a44*fwbc(i,j)&
  !              +a45*fwbc(i01,j)+a46*fwbc(i02,j) 
  !            d2=b51*fwbc(i30,j)+b52*fwbc(i20,j)+b53*fwbc(i10,j)+b54*fwbc(i,j)&
  !              +b55*fwbc(i01,j)+b56*fwbc(i02,j)
  !            c1=c1*(cl2-4.) 
  !            if (ihol.eq.5.or.i03.lt.0.or.i40.gt.ir+1) then
  !               d2=d2*3.
  !               flw5=flw4+c1*(d1-xsign*d2)   ! 5th order interflux
  !               fhol(i,j)=flw5
  !               goto 7777        ! end of finding fhol
  !            endif
  !         endif
  !    ! 6th order
  !         if (ihol.ge.6) then
  !            flw6=flw4+c1*(d1-cl1*d2)   ! 6th order interflux
  !            if (ihol.eq.6) then
  !               fhol(i,j)=flw6
  !               goto 7777        ! end of finding fhol
  !            endif
  !         endif
  !    ! 7th order
  !         if (ihol.eq.7) then
  !            d1=a61*fwbc(i40,j)+a62*fwbc(i30,j)+a63*fwbc(i20,j)+a64*fwbc(i10,j)&
  !              +a65*fwbc(i,j)+a66*fwbc(i01,j)+a67*fwbc(i02,j)+a68*fwbc(i03,j)
  !            d2=b71*fwbc(i40,j)+b72*fwbc(i30,j)+b73*fwbc(i20,j)+b74*fwbc(i10,j)&
  !              +b75*fwbc(i,j)+b76*fwbc(i01,j)+b77*fwbc(i02,j)+b78*fwbc(i03,j)
  !            c1=c1*(cl2-9.) 
  !            flw7=flw6+c1*(d1-xsign*d2)   ! 7th order interflux
  !            fhol(i,j)=flw7
  !         endif
       ! 2nd order
            flw2=0.5*((fwbc(i10,j)+fwbc(i,j))+cl1*(fwbc(i,j)-fwbc(i10,j)))
            if (i20.gt.ir+1) then
               fhol(i,j)=flw2
               goto 7777        ! end of finding fhol
            endif
       ! 3rd order
            d1=a21*fwbc(i20,j)+a22*fwbc(i10,j)+a23*fwbc(i,j)+a24*fwbc(i01,j)
            d2=b31*fwbc(i20,j)+b32*fwbc(i10,j)+b33*fwbc(i,j)+b34*fwbc(i01,j)
            if (ihol.eq.3.or.i02.lt.0.or.i30.gt.ir+1) then
               flw3=flw2+(cl2-1.)*(d1-xsign*2.*d2)   ! 3rd order flux
               fhol(i,j)=flw3
               goto 7777        ! end of finding fhol
            endif
       ! 4th order
            flw4=flw2+(cl2-1.)*(d1-cl1*d2)   ! 4th order flux
       ! 5th order
            d1=a41*fwbc(i30,j)+a42*fwbc(i20,j)+a43*fwbc(i10,j)+a44*fwbc(i,j)&
              +a45*fwbc(i01,j)+a46*fwbc(i02,j)
            d2=b51*fwbc(i30,j)+b52*fwbc(i20,j)+b53*fwbc(i10,j)+b54*fwbc(i,j)&
              +b55*fwbc(i01,j)+b56*fwbc(i02,j)
            if (ihol.eq.5.or.i03.lt.0.or.i40.gt.ir+1) then
               flw5=flw4+(cl2-1.)*(cl2-4.)*(d1-xsign*3.*d2)   ! 5th order inter flux
               fhol(i,j)=flw5
               goto 7777        ! end of finding fhol
            endif
       ! 6th order
            flw6=flw4+(cl2-1.)*(cl2-4.)*(d1-cl1*d2)   ! 6th order inter flux
       ! 7th order
            d1=a61*fwbc(i40,j)+a62*fwbc(i30,j)+a63*fwbc(i20,j)+a64*fwbc(i10,j)&
              +a65*fwbc(i,j)+a66*fwbc(i01,j)+a67*fwbc(i02,j)+a68*fwbc(i03,j)
            d2=b71*fwbc(i40,j)+b72*fwbc(i30,j)+b73*fwbc(i20,j)+b74*fwbc(i10,j)&
              +b75*fwbc(i,j)+b76*fwbc(i01,j)+b77*fwbc(i02,j)+b78*fwbc(i03,j)
            flw7=flw6+(cl2-1.)*(cl2-4.)*(cl2-9.)*(d1-xsign*d2)
            fhol(i,j)=flw7
7777        continue        ! end of finding fhol      

       ! find f*p
       ! 1st order
            xsign=xsignp(i,j)
            if (xsign.eq.1.) then
               fup=fwbc(i,j)
            else
               fup=fwbc(i,j10)
            endif
            fupp(i,j)=fup
  !    ! 2nd order (Lax-Wendroff)
  !         flw2=0.5*((fwbc(i,j10)+fwbc(i,j))+cp1*(fwbc(i,j)-fwbc(i,j10)))   
  !         if (ihop.eq.2) then
  !         ! Superbee limter
  !            x=fwbc(i,j10)-fwbc(i,j)
  !            if (abs(x).le.1.e-27) fhop(i,j)=fup
  !            if (abs(x).gt.1.e-27) then
  !               if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i,j01))/x
  !               if (xsign.eq.-1.) r=(fwbc(i,j20)-fwbc(i,j10))/x
  !               if (r.le.0.) fhop(i,j)=fup
  !               if (r.gt.0.) then
  !                  xlimiter=max(min(2.*r,1.),min(r,2.))
  !                  corr=flw2-fup
  !                  fhop(i,j)=fup+xlimiter*corr
  !                  if (fhop(i,j).lt.0.) fhop(i,j)=fup     ! fsbp can't be < 0
  !               endif
  !            endif
  !         endif
  !    ! 3rd order (QUICKEST)
  !         if (ihop.ge.3) then
  !            d1=a21*fwbc(i,j20)+a22*fwbc(i,j10)+a23*fwbc(i,j)+a24*fwbc(i,j01) 
  !            d2=b31*fwbc(i,j20)+b32*fwbc(i,j10)+b33*fwbc(i,j)+b34*fwbc(i,j01)
  !            c1=cp2-1. 
  !            if (ihop.eq.3) then
  !               d2=2.*d2
  !               flw3=flw2+c1*(d1-xsign*d2)   ! 3th order inter-flux
  !               fhop(i,j)=flw3
  !            endif
  !         endif
  !    ! 4th order
  !         if (ihop.ge.4) then
  !            flw4=flw2+c1*(d1-cp1*d2)   ! 4th order flux
  !            if (ihop.eq.4) fhop(i,j)=flw4
  !         endif
  !    ! 5th order
  !         if (ihop.ge.5) then
  !            d1=a41*fwbc(i,j30)+a42*fwbc(i,j20)+a43*fwbc(i,j10)+a44*fwbc(i,j)&
  !              +a45*fwbc(i,j01)+a46*fwbc(i,j02) 
  !            d2=b51*fwbc(i,j30)+b52*fwbc(i,j20)+b53*fwbc(i,j10)+b54*fwbc(i,j)&
  !              +b55*fwbc(i,j01)+b56*fwbc(i,j02) 
  !            c1=c1*(cp2-4.) 
  !            if (ihop.eq.5) then
  !               d2=3.*d2
  !               flw5=flw4+c1*(d1-xsign*d2)   ! 6th order interflux
  !               fhop(i,j)=flw5
  !            endif
  !         endif
  !    ! 6th order
  !         if (ihop.ge.6) then
  !            flw6=flw4+c1*(d1-cp1*d2)   ! 6th order interflux
  !            if (ihop.eq.6) fhop(i,j)=flw6
  !         endif
  !    ! 7th order
  !         if (ihop.eq.7) then
  !            d1=a61*fwbc(i,j40)+a62*fwbc(i,j30)+a63*fwbc(i,j20)+a64*fwbc(i,j10)&
  !              +a65*fwbc(i,j)+a66*fwbc(i,j01)+a67*fwbc(i,j02)+a68*fwbc(i,j03)
  !            d2=b71*fwbc(i,j40)+b72*fwbc(i,j30)+b73*fwbc(i,j20)+b74*fwbc(i,j10)&
  !              +b75*fwbc(i,j)+b76*fwbc(i,j01)+b77*fwbc(i,j02)+b78*fwbc(i,j03)
  !            c1=c1*(cp2-9.) 
  !            flw7=flw6+c1*(d1-xsign*d2)   ! 7th order interflux
  !            fhop(i,j)=flw7
  !         endif
       ! 2nd order
            flw2=0.5*((fwbc(i,j10)+fwbc(i,j))+cp1*(fwbc(i,j)-fwbc(i,j10)))
       ! 3rd order
            d1=a21*fwbc(i,j20)+a22*fwbc(i,j10)+a23*fwbc(i,j)+a24*fwbc(i,j01)
            d2=b31*fwbc(i,j20)+b32*fwbc(i,j10)+b33*fwbc(i,j)+b34*fwbc(i,j01)
            if (ihop.eq.3) then
               flw3=flw2+(cp2-1.)*(d1-xsign*2.*d2)   ! 3rd order flux
               fhop(i,j)=flw3
               goto 8888
            endif
       ! 4th order
            flw4=flw2+(cp2-1.)*(d1-cp1*d2)   ! 4th order flux
       ! 5th order
            d1=a41*fwbc(i,j30)+a42*fwbc(i,j20)+a43*fwbc(i,j10)+a44*fwbc(i,j)&
              +a45*fwbc(i,j01)+a46*fwbc(i,j02)
            d2=b51*fwbc(i,j30)+b52*fwbc(i,j20)+b53*fwbc(i,j10)+b54*fwbc(i,j)&
              +b55*fwbc(i,j01)+b56*fwbc(i,j02)
            if (ihop.eq.5) then
               flw5=flw4+(cp2-1.)*(cp2-4.)*(d1-xsign*3.*d2)   ! 5th order 
               fhop(i,j)=flw5
               goto 8888
            endif
       ! 6th order
            flw6=flw4+(cp2-1.)*(cp2-4.)*(d1-cp1*d2)   ! 6th order inter flux
       ! 7th order
            d1=a61*fwbc(i,j40)+a62*fwbc(i,j30)+a63*fwbc(i,j20)+a64*fwbc(i,j10)&
              +a65*fwbc(i,j)+a66*fwbc(i,j01)+a67*fwbc(i,j02)+a68*fwbc(i,j03)
            d2=b71*fwbc(i,j40)+b72*fwbc(i,j30)+b73*fwbc(i,j20)+b74*fwbc(i,j10)&
              +b75*fwbc(i,j)+b76*fwbc(i,j01)+b77*fwbc(i,j02)+b78*fwbc(i,j03)
            flw7=flw6+(cp2-1.)*(cp2-4.)*(cp2-9.)*(d1-xsign*d2)   ! 7th order
            fhop(i,j)=flw7
8888        continue        ! end of finding fhop      
         enddo              ! end of do i=1,ir
      enddo                 ! end of do j=1,ip

! ULTIMATE universial limiter
      do j=1,ip
         j01=j-1
         j10=j+1
         j20=j+2
         if (j01.lt.1) j01=j01+ip
         if (j10.gt.ip) j10=j10-ip
         if (j20.gt.ip) j20=j20-ip
         ib=max(iba(j),iba(j10))
         ib=min(ib,ir-1)
         do i=1,ib
         ! f*l
            if (ihol.ge.3) then
               cl1=abs(cl(i,j))
               cp1=abs(cp(i,j))
               i01=i-1
               i10=i+1
               i20=i+2
               if (cl(i,j).gt.0.) then
                  fu=fwbc(i01,j)
                  fc=fwbc(i,j)
                  fd=fwbc(i10,j)
               else
                  fu=fwbc(i20,j)
                  fc=fwbc(i10,j)
                  fd=fwbc(i,j)
               endif
               del=fd-fu
               adel=abs(del)
               acurv=abs(fd+fu-fc-fc)
               if (acurv.ge.adel) then
                  fhol(i,j)=fc
               else
                  ref=fu+(fc-fu)/cl1
                  if (del.gt.0.) then
                     flim=max(fhol(i,j),fc)
                     fhol(i,j)=min(flim,min(ref,fd))
                  else
                     flim=max(fhol(i,j),max(ref,fd))
                     fhol(i,j)=min(flim,fc)
                  endif
               endif
            endif
         ! f*p
            if (ihop.ge.3) then
               if (cp(i,j).gt.0.) then
                  fu=fwbc(i,j01)
                  fc=fwbc(i,j)
                  fd=fwbc(i,j10)
               else
                  fu=fwbc(i,j20)
                  fc=fwbc(i,j10)
                  fd=fwbc(i,j)
               endif
               del=fd-fu
               adel=abs(del)
               acurv=abs(fd+fu-fc-fc)
               if (acurv.ge.adel) then
                  fhop(i,j)=fc
               else
                  ref=fu+(fc-fu)/cp1
                  if (del.gt.0.) then
                     flim=max(fhop(i,j),fc)
                     fhop(i,j)=min(flim,min(ref,fd))
                  else
                     flim=max(fhop(i,j),max(ref,fd))
                     fhop(i,j)=min(flim,fc)
                  endif
               endif
            endif
         enddo              ! end of do i=1,ir
      enddo                 ! end of do j=1,ip

      end subroutine FLS_HO_2D


!******************************************************************************
  subroutine make_igrfvol(nlg,mlon,volg,xlatg)
!******************************************************************************

  use constants
  use cimigrid_dim
  use cread1, only: rc,outname
  use cread2, only: hlosscone,iyear,iday
  implicit none
  external IGRF_GSW_08,zeroB
  integer,parameter :: itmax=40  ! max no. of iteration in finding zero dip
  integer,parameter :: np=1000
  real,parameter :: errb=1.e-5      ! tolerance in magnetic field
  real,parameter :: errd=1.e-8      ! tolerance in dip latitude
  real,parameter :: dir=1.,dsmax=0.2,err=0.0001
  integer j,n,npf,npf1,m,nmax(ip),nm,nlg
  real xlatN,xlatS,xlat0,gcolat,glon,Btheta,Bphi, &
       Br,mlon(ip),xlatdip(ip),volg(nlg),parmod(10), &
       xmag,ymag,zmag,xgeo,ygeo,zgeo,xlat(nlg), &
       xlatg(nlg,ip),xL1,xL2,xL,dL,coslat,xsm,ysm,zsm,xgsm,ygsm,zgsm,rf, &
       xf,yf,zf,xa(np),ya(np),za(np),ba(np),dss(np),yint(np),sso,bx,by,bz, &
       b_mid,volumeg(nlg,ip),xlat1,rlim

! Initial setup
  parmod(:)=1.               ! parmod's are not used in calculation
  open(unit=24,file=trim(outname)//'.igrfvol')

! Find magnetic dip at mlon, xlatdip(ip)
  do j=1,ip
     xlatN=40.12     ! Northern bound
     xlatS=-39.      ! Southern bound
     xlat0=0.5*(xlatN+xlatS)
     find_dip: do n=1,itmax
        xmag=rc*cosd(xlat0)*cosd(mlon(j))
        ymag=rc*cosd(xlat0)*sind(mlon(j))
        zmag=rc*sind(xlat0)
        call GEOMAG_08(xgeo,ygeo,zgeo,xmag,ymag,zmag,m_one)
        gcolat=acos(zgeo/rc)
        glon=atan2(ygeo,xgeo)
        call igrf_geo_08(rc,gcolat,glon,Br,Btheta,Bphi)
        if (abs(Br).lt.errb) exit find_dip
        if (Br.lt.0.) xlatN=xlat0
        if (Br.gt.0.) xlatS=xlat0
        xlat0=0.5*(xlatN+xlatS)
        if (abs(xlat0-xlatN).lt.errd) exit find_dip
     enddo find_dip
     xlatdip(j)=xlat0
  enddo

! Calculate xlat(1:nlg)
  xL1=1.15
  xL2=15.
  dL=(xL2-xL1)/(nlg-1)
  do n=1,nlg
     xL=xL1+(n-1)*dL
     coslat=sqrt(1./xL)
     xlat(n)=acos(coslat)*180./pi
  enddo

! Calculate flux tube volume as a function of xlat along mlon longitude grid
  rlim=1.5*xL2
  do j=1,ip
     nmax(j)=nlg
     do n=1,nlg
        xlat0=xlatdip(j)+xlat(n)
        xmag=rc*cosd(xlat0)*cosd(mlon(j))
        ymag=rc*cosd(xlat0)*sind(mlon(j))
        zmag=rc*sind(xlat0)
        call magsm_08(xmag,ymag,zmag,xsm,ysm,zsm,i_one)           ! mag to sm
        call SMGSW_08(xsm,ysm,zsm,xgsm,ygsm,zgsm,i_one)           ! sm to gsm
        call trace_08(xgsm,ygsm,zgsm,dir,dsmax,err,rlim,rc,i_one,parmod, &
                      zeroB,IGRF_GSW_08,xf,yf,zf,xa,ya,za,npf,np)
        rf=sqrt(xf*xf+yf*yf+zf*zf)
        if (rf.gt.xL1) then     ! tracing doesn't stop at south ionosphere
           nmax(j)=n-1
           volumeg(n:nlg,j)=volumeg(n-1,j)
           goto 777      ! do next j
        endif
        do m=1,npf
           call IGRF_GSW_08(xa(m),ya(m),za(m),bx,by,bz)
           ba(m)=sqrt(bx*bx+by*by+bz*bz)*1.e-9    ! B field in Tesla
        enddo
        npf1=npf-1
        do m=1,npf1
           dss(m)=sqrt((xa(m)-xa(m+1))**2+(ya(m)-ya(m+1))**2+(za(m)-za(m+1))**2)
           b_mid=0.5*(ba(m)+ba(m+1))
           yint(m)=1./b_mid
        enddo
        call closed(i_one,npf1,yint,dss,sso)      ! use closed form
        volumeg(n,j)=sso*re_m  ! volume per unit flux
     enddo
777  continue
  enddo

! Find constant volume volg
  volg(1)=maxval(volumeg(1,:))
  volg(nlg)=minval(volumeg(nlg,:))
  xL1=volg(1)**0.25
  xL2=volg(nlg)**0.25
  dL=(xL2-xL1)/(nlg-1)
  do n=2,nlg-1
     xL=xL1+(n-1)*dL
     volg(n)=xL**4
  enddo

! Find xlatg at volume volg
  do j=1,ip
     nm=nmax(j)
     do n=1,nlg
        call lintp(volumeg(1:nm,j),xlat(1:nm),nm,volg(n),xlat1)
        xlatg(n,j)=xlatdip(j)+xlat1
     enddo
  enddo

! Write results in outname.igrfvol file
  write(24,*) iyear,iday,'   ! iyear,iday'
  write(24,*) nlg,ip,'   ! nlg,ip'
  write(24,*) 'mlon(1:ip), magnetic longitude in degree'
  write(24,*) mlon
  write(24,*) 'xlatdip(1:ip), magnetic dip in degree'
  write(24,*) xlatdip
  write(24,*) 'volg(1:nlg), Flux tube volume/magnetic flux in m/Tesla'
  write(24,*) volg
  write(24,*) 'xlatg(1:nlg,1:ip), MLAT in degree with volume=volg'
  write(24,*) xlatg
  close(24)

  end subroutine make_igrfvol


!******************************************************************************
  subroutine corPotentgeo(cor,nlg,npg,thetag,xlongg,potentg)
!******************************************************************************
! Subroutine calculates corotation potential on fixed geographic grid,
! (thetag,xlongg)  

  use constants, only: re_m,pi
  use cread1, only: rc
  implicit none
  integer nlg,npg,i,j
  real cor,thetag(nlg),xlongg(npg),potentg(nlg,npg)
  real rc2,thmax,dtheta,dphig,thetai,Bri1,Btheta,Bphi,ksaig

  rc2=rc*rc*re_m*re_m

! Setup the fixed geograhic grid
  thmax=1.3    ! max colat in radian
  dtheta=thmax/(nlg-1)
  do i=1,nlg
     thetag(i)=(i-1)*dtheta
  enddo
  dphig=2.*pi/(npg-1)
  do j=1,npg
     xlongg(j)=-pi+(j-1)*dphig     ! geographic longitude from -pi -> pi
  enddo

! Find potentg
  do j=1,npg
     potentg(1,j)=0.
     do i=2,nlg
        thetai=0.5*(thetag(i-1)+thetag(i))
        call IGRF_GEO_08(rc,thetai,xlongg(j),Bri1,Btheta,Bphi)        
        ksaig=1.e-9*Bri1*rc2*sin(thetai)
        potentg(i,j)=potentg(i-1,j)+cor*ksaig*dtheta
     enddo
  enddo

  end subroutine corPotentgeo
