MODULE CONSTANTS
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
      real*8, parameter :: Re2 = Re**2                      ! Square of Earth radius [m^2]
      real*8, parameter :: mEarth = 5.9722e24               ! Earth mass [kg]
      real*8, parameter :: constG = 6.6743e-11              ! Gravitational constant
      real*8, parameter :: GM = constG * mEarth             ! For convenience
      real*8, parameter :: arad = 0.1774d-2                 ! For radiation pressure [m/s^2]
      real*8, parameter :: Wrot = 1.9910d-7                 ! Earth's angular speed [rad/s]

END MODULE CONSTANTS



MODULE SETTING

   USE CONSTANTS
   IMPLICIT NONE
!!********************<< User Setting >>********************!!
!!**** Recommended to change the parameters for your run ***!!

!! TAG for an Unique Runname !! ex. output_filename = "MATE_nH_GRC_{tag}_2008174.data" !! Example
   character*20, parameter :: Runname_in_10char              = "CXtest1"

!! TIME SETTING !!
   integer, parameter ::      Start_Time_in_YYYYDOY          = 1000010
   integer, parameter ::      End_Time_in_YYYYDOY            = 1000010
   integer, parameter ::      Output_Time_Interval_in_Minute = 60*24

!! 3-D SPATIAL RESOLUTIONS !!
   integer, parameter ::      GEO_Resolution_in_Degree       = 15
   real*8, parameter  ::      RadialRange_min                = 2.0                         ! unit Re
   real*8, parameter  ::      RadialRange_max                = 10.0
   real*8, parameter  ::      dR                             = 0.5

!! Number of Particle Direction !!
   integer, parameter ::      dV_theta                         = 6                         ! Angular Resolution of Velocity direction in Degree: DEFAULT = 6  (Recommended Values = [2, 3, 4, 6, 10])
   integer, parameter ::      nTheta                           = 180/dV_theta + 1          ! # of theta grids

!! (FIX ME!!) Energy Grid !!
   integer, parameter ::      nEnergy                          = 61                        ! # of energy grid

!! DIRECTORIES SETTING !!
   character*70, parameter :: Lya_dir                          = "OMNI_extended.txt"         ! 1964 - 2024 (Oct)
   character*70, parameter :: BC_dir                           = "/nobackup/slee122/MATE/MSIS/BC/"
   character*70, parameter :: outdir                           = "/nobackup/slee122/MATE/0728/"
!   character*70, parameter :: outdir                           = "/home/sylee/exospherecode/MATE/output/0728/"

!! Including Physics !!
   integer, parameter ::      i_EarthGravity                   = 1                          ! 0 for turn off / 1 for turn on
   integer, parameter ::      i_SolarRadiationPressure         = 1                          ! If you set these values as N, then the force is N-times stronger as a coefficient.
   integer, parameter ::      i_CoriolisForce_GSE              = 1
   integer, parameter ::      i_Photoionization                = 0
   integer, parameter ::      i_ChargeExchange                 = 1                          !

!! Exobase Boundary Condition (BC) Setting !!
   character*10, parameter :: ExobaseBC_Model_Name             = "CONST"                  ! "MSIS", "TIMEGCM", "WACCMX", "CONST"
   integer, parameter ::      BC_GEO_Resolution_in_Degree      = 5
   integer, parameter ::      BC_Time_Resolution_in_Minute     = 5
   integer, parameter ::      BC_Start_Time_in_YYYYDOY         = 1000010                    ! BC covers from this time: A few days before {Start_Time_in_YYYYDOY} for tracing.

!! Particle Tracing Range !!
   real*8, parameter :: inner_boundary                         = Re + 500.d3                ! Exobase location
   real*8, parameter :: outer_boundary                         = 100*Re
   real*8, parameter :: Max_Travel_Time_in_Days                = 30                         ! DEFAULT = 60

   integer, parameter :: nstep = 50000 ! 10000                                               ! Maximum number of steps for particle tracing (typical maximum is 3e4)


!!**************<<Rearrange - Not recommended to revise below>>**************!!
!!***** The parameters below are determined by the values defined above *****!!

   integer, parameter :: start_ydoy = Start_Time_in_YYYYDOY
   integer, parameter :: end_ydoy   = End_Time_in_YYYYDOY
   character*20, parameter :: tag0 = Runname_in_10char

   integer, parameter :: nRadial = (RadialRange_max-RadialRange_min)/dR + 1  

   integer, parameter :: geores = GEO_Resolution_in_Degree
   integer, parameter :: nLat = 90/geores +1
   integer, parameter :: nLat_NS = nLat*2-1              ! # of latitude grid for both hemisphere
   integer, parameter :: nLong = 360/geores
   integer, parameter :: nLon = nLong                    ! To avoid confusion...^^;

   integer, parameter :: start_year = start_ydoy/1000
   integer, parameter :: end_year   = end_ydoy/1000
   integer, parameter :: start_doy  = mod(start_ydoy,1000)
   integer, parameter :: end_doy    = mod(end_ydoy,1000)
   integer, parameter :: ndays = end_ydoy - start_ydoy +1            ! Revise it if end_year =/ start_year
   
   integer, parameter :: bc_res = BC_GEO_Resolution_in_Degree      ! 5-degree resolution for both LON & LAT
   integer, parameter :: tb_res = BC_Time_Resolution_in_Minute*60  ! Time resolution for BC (5-min resolution = 300s)
   integer, parameter :: nbx = 360/bc_res, nby = 180/bc_res        ! 360-longitude and 180-latitude degree
   integer, parameter :: nbtperday = 86400/tb_res                  ! Total number of time grid for BC
   integer, parameter :: nt_bwd_bc = Start_Time_in_YYYYDOY - BC_Start_Time_in_YYYYDOY    ! # of days earlier than start_doy for BC

   integer, parameter :: time_resolution = 60*Output_Time_Interval_in_Minute     ! Unit second
   integer, parameter :: ntperday = 86400/time_resolution              ! Number of time grid per day
   integer, parameter :: nt = ntperday*ndays                           ! Total number of time grid for nH
   integer, parameter :: ntmax = Max_Travel_Time_in_Days
   real*8,  parameter :: tmax = ntmax * 86400.d0          ! 60 days
   integer, parameter :: n_physics = 5                               ! # of physics (Gravity, ..., ChargeExchange)

   integer, parameter :: previous_year = end_year - 1                ! For daily-varying indices (e.g. Lya, F10.7),
   integer, parameter :: start_ydoy_index = previous_year*1000 + 1   ! Load two-year data
   integer, parameter :: end_ydoy_index = end_year*1000 + 366        ! By assuming the simulation perious is < 1 year.

   ! Not parameters
   integer, save :: nvel   ! Number of velocity direction defined in SET_VELOCITY_DIRECTION


END MODULE SETTING



