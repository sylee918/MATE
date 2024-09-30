      include "constants.inc"

!!********************<< User Setting >>********************!!
!!**** Recommended to change the parameters for your run ***!!

   !! TAG for an Unique Runname !!
      character*10, parameter :: TAG_within_10_characters       = "test"                      ! ex: output_filename = "MATE_nH_GRC_{tag}_2008174.data"

   !! TIME SETTING !!
      integer, parameter ::      Start_Time_in_YYYYDOY          = 2008164
      integer, parameter ::      End_Time_in_YYYYDOY            = 2008164
      integer, parameter ::      Output_Time_Interval_in_Minute = 60

   !! 3-D SPATIAL RESOLUTIONS !!
      integer, parameter ::      GEO_Resolution_in_Degree       = 5
      real*8, parameter  ::      RadialRange_min                = 2.0                         ! unit Re
      real*8, parameter  ::      RadialRange_max                = 15.0
      real*8, parameter  ::      dR                             = 0.5

   !! Number of Particle Direction !!
      integer, parameter ::      dV_theta                         = 6                         ! Angular Resolution of Velocity direction in Degree: DEFAULT = 6  (Recommended Values = [2, 3, 4, 6, 10])
      integer, parameter ::      nTheta                           = 180/dV_theta + 1          ! # of theta grids
      integer, parameter ::      N_vel_directions                 = 1148                      ! # of velocity directions. [10318, 4586, 2586, 1148, 412] for [2, 3, 4, 6, 10]

   !! (FIX ME!!) Energy Grid !!
      integer, parameter :: nEnergy                               = 61                        ! # of energy grid

   !! DIRECTORIES SETTING !!
      character*70, parameter :: Lya_dir                          = "indices_for_MSIS_1963-2023_2.txt"
      character*70, parameter :: BC_dir                           = "./exobaseBC/"
      character*70, parameter :: outdir                           = "/nobackup/slee122/MATE/0826/GRC/"

   !! Including Physics !!
      integer, parameter ::      i_EarthGravity                   = 1                          ! 0 for turn off / 1 for turn on
      integer, parameter ::      i_SolarRadiationPressure         = 1                          ! If you set these values as N, then the force is N-times stronger as a coefficient.
      integer, parameter ::      i_CoriolisForce_GSE              = 1
      integer, parameter ::      i_Photoionization                = 1
      integer, parameter ::      i_ChargeExchange                 = 0                          ! Not developed yet

   !! Exobase Boundary Condition (BC) Setting !!
      character*10, parameter :: ExobaseBC_Model_Name             = "TIMEGCM"                  ! "MSIS", "TIMEGCM", "WACCMX"
      integer, parameter ::      BC_GEO_Resolution_in_Degree      = 5
      integer, parameter ::      BC_Time_Resolution_in_Minute     = 5
      integer, parameter ::      BC_Start_Time_in_YYYYDOY         = 2008164                    ! BC covers from this time: A few days before {Start_Time_in_YYYYDOY} for tracing.

   !! Particle Tracing Range !!
      real*8, parameter :: inner_boundary                         = Re + 500.d0                ! Exobase location
      real*8, parameter :: outer_boundary                         = 100*Re
      real*8, parameter :: Max_Travel_Time_in_Days                = 60                         ! DEFAULT = 60



!!**************<<Rearrange - Not recommended to revise below>>**************!!
!!***** The parameters below are determined by the values defined above *****!!

      integer, parameter :: start_ydoy = Start_Time_in_YYYYDOY
      integer, parameter :: end_ydoy   = End_Time_in_YYYYDOY
      character*10, parameter :: tag0 = TAG_within_10_characters

      integer, parameter :: nRadial = (RadialRange_max-RadialRange_min)/dR   

      integer, parameter :: geores = GEO_Resolution_in_Degree
      integer, parameter :: nLat = 180/geores +1
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
      integer, parameter :: n_physics = 5                               ! # of physics (Gravity, ..., ChargeExchange)

      integer, parameter :: previous_year = end_year - 1                ! For daily-varying indices (e.g. Lya, F10.7),
      integer, parameter :: start_ydoy_index = previous_year*1000 + 1   ! Load two-year data
      integer, parameter :: end_ydoy_index = end_year*1000 + 366        ! By assuming the simulation perious is < 1 year.
