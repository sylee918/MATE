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
  integer,parameter :: nPA=ig
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
  
  module cgrid
      use cimigrid_dim
      ! Scalar variables
      real :: dlnk, dvarL, dphi
      ! 1D arrays
      real, dimension(ir+1) :: varL
      real, dimension(irh) :: xlath
      real, dimension(ik) :: xkb
      real, dimension(0:ik+1) :: xk
      real, dimension(ip) :: mphi, mlon
      real, dimension(ns) :: xmass, dlnp, d4
      real, dimension(ig) :: gridy, cosSq, sinSq
      real, dimension(5) :: xmass1
      real, dimension(je) :: gride_e, gride_i
      ! 2D arrays
      real, dimension(ir,ip) :: xlati, potentc
      real, dimension(ir+1,ip) :: dlati, ksai
      real, dimension(0:ir+1,ip) :: xlatd
      real, dimension(ns,je) :: gride
      real, dimension(ns,0:je) :: ebound
      real, dimension(ns,0:iw+1) :: gridp, lnp, ekev
      real, dimension(ns,iw) :: dkeV, pcEo, vel
      ! 5D arrays
      real, dimension(ir,ip,je,ig) :: fl
      real, dimension(ir,ip) :: density
  end module

  module cfield
    use cimigrid_dim
    ! Scalar variables
    integer, dimension(ip) :: iba
    integer, dimension(ir,ip) :: mcN, mcS
    integer :: iday2, ihour, jnoon
    real :: Dst, DstRC, xL1, xL2, psi1, zkp, F107, Apt, dsmax, err
    real, dimension(ip) :: phi, xmlt, rsb, xmltb, bob
    real, dimension(10) :: parmod0, parmod
    real, dimension(ir,ip) :: bo, ro, xmlto, xmltS, xo, yo, volume, xkcN, xkcS, xlatiS, phiS, BiN, BiS, BriN, BriS, mlonS, sini_f
    real, dimension(ir,ip,ik) :: Hdens, dmu, dlnBmdL, dlnBmdp, dlnBmdL1, dlnBmdp1, lnbm
    real, dimension(ir,ip,0:ik+1) :: y, bm, tya
    real, dimension(ns,ir,ip,iw,ik) :: fcone, Tbounce
    real, dimension(ns,ir+1,ip,0:iw+1,0:ik+1) :: xjac
  end module

  module useless
    use cimigrid_dim
    real, dimension(ir,ip) :: density, ompe, CHpower, HIpower
    real, dimension(ns,ir,ip) :: denWP, TparaWP, TperpWP, HRPee, HRPii
    real, dimension(ip) :: rppa
    real :: Lstar_max(0:ik), Lstar(ir,ip,0:ik)
  end module

  Module MATEgrid
    use constants
    integer, parameter :: nRadial=17, nLon=72, nLat=37, nt=24  ! 0.5 Re, 5-deg, 1-hour resolution
    real, parameter :: drMATE=0.5, dangleMATE=5.0, dtMATE=1.0
    ! r: 2.0 to 10.0 Re with 0.5 Re intervals (17 values)
    real, parameter, dimension(nRadial) :: rMATE = [(real(i-1)*drMATE+2.0, i=1,nRadial)]
    ! lon: 0 to 360 degrees with 5-degree intervals (73 values)
    real, parameter, dimension(nLon) :: lonMATE = [(real(i-1)*dangleMATE, i=1,nLon)]
    ! lat: -90 to 90 degrees with 5-degree intervals (37 values)
    real, parameter, dimension(nLat) :: latMATE = [(real(i-1)*dangleMATE - 90.0, i=1,nLat)]
    real, parameter, dimension(nt) :: tMATE = [(real(i-1)*dtMATE, i=1,nt)]
    real, dimension(nRadial,nLon,nLat,nt) :: betaMATE, nPS_MATE
    real, dimension(nRadial,nLon,nLat) :: weight1
  End Module