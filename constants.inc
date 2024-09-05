!  Physical Constants
      real*8, parameter :: pi = 3.141592653589793
      real*8, parameter :: c = 2.99792458d8                 ! Speed of light [m/s]
      real*8, parameter :: e = 1.60217646d-19               ! Elementary charge [C]
      real*8, parameter :: kb = 1.38065030d-23              ! Boltzman constant [J/K]
      real*8, parameter :: kbeV = kb/e                      ! Boltzman constant [eV/K]
      real*8, parameter :: h = 6.626d-34                    ! Planck constant
      real*8, parameter :: mH = 1.6735575e-27               ! Hydrogen atom mass [kg]
      real*8, parameter :: mp = 1.67262158d-27              ! Proton mass [kg]
      real*8, parameter :: me = 9.10938188d-31              ! Electron mass [kg]
      real*8, parameter :: Re = 6.371009d6                  ! Earth radius [m]
      real*8, parameter :: mEarth = 5.9722e24               ! Earth mass [kg]
      real*8, parameter :: constG = 6.6743e-11              ! Gravitational constant
      real*8, parameter :: GM = constG * mEarth             ! For convenience
      real*8, parameter :: arad = 0.1774d-2                 ! For radiation pressure [m/s^2]
      real*8, parameter :: Wrot = 1.9910d-7                 ! Earth's angular speed [rad/s]

!  Simulational Constants (e.g. Number of grids)
!      integer, parameter :: nTheta = 91                    ! # of theta grids
!      integer, parameter :: N_vel_directions = 10318       ! # of velocity directions determined by nTheta
!      integer, parameter :: nEnergy = 121                   ! # of energy grid
      integer, parameter :: nTheta = 19                     ! # of theta grids
      integer, parameter :: N_vel_directions = 412         ! # of velocity directions determined by nTheta
      integer, parameter :: nRadial = 27                    ! # of grids for radial distance (17 for 10 Re & 37 for 20 Re)
      integer, parameter :: nEnergy = 61                   ! # of energy grid
      integer, parameter :: nLat = 7                        ! # of latitude grid for Northern hemisphere including lat=0
      integer, parameter :: nLat_NS = nLat*2-1              ! # of latitude grid for both hemisphere
      integer, parameter :: nLong = 24                      ! # of longitude grid (24 for '' )
      integer, parameter :: nLon = nLong                    ! To avoid confusion...^^;

      integer, parameter :: start_ydoy = 2014174
      integer, parameter :: end_ydoy   = 2014200
      integer, parameter :: start_year = start_ydoy/1000
      integer, parameter :: end_year   = end_ydoy/1000
      integer, parameter :: start_doy  = mod(start_ydoy,1000)
      integer, parameter :: end_doy    = mod(end_ydoy,1000)
      integer, parameter :: ndays = end_ydoy - start_ydoy +1            ! Revise it if end_year =/ start_year
      integer, parameter :: ntmax = 60                                  ! # of days of tmax for tracing (n_tmax, not nt_max)
      integer, parameter :: nt_bwd_bc = 10                              ! # of days earlier than start_doy for BC

      integer, parameter :: previous_year = end_year - 1                ! For daily-varying indices (e.g. Lya, F10.7),
      integer, parameter :: start_ydoy_index = previous_year*1000 + 1   ! Load two-year data
      integer, parameter :: end_ydoy_index = end_year*1000 + 366        ! By assuming the simulation perious is < 1 year.

      integer, parameter :: tres = 60*60                         ! Time resolution for nH (30-min resolution = 1800s)
      integer, parameter :: ntperday = 86400/tres                ! Number of time grid per day
      integer, parameter :: nt = ntperday*ndays                  ! Total number of time grid for nH

!  Grid for Exobase Boundary Condition
      integer, parameter :: bc_res = 5                            ! 5-degree resolution for both LON & LAT
      integer, parameter :: tb_res = 5*60                         ! Time resolution for BC (5-min resolution = 300s)
      integer, parameter :: nbx = 360/bc_res, nby = 180/bc_res    ! 360-longitude and 180-latitude degree
      integer, parameter :: nbtperday = 86400/tb_res              ! Total number of time grid for BC

!  Setting
      character*10, parameter :: exobaseBC_type = "TIMEGCM"       ! "MSIS", "TIMEGCM", "WACCMX"

! Directories...
      character*70, parameter :: Lya_dir = "indices_for_MSIS_1963-2023_2.txt"
      character*70, parameter :: BC_dir ="./exobaseBC/"
      character*70, parameter :: outdir = "/nobackup/slee122/MATE/0826/GRC/"
