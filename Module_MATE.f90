   Module INTEGRATED_INITIALIZATION

      USE SET_VELOCITY_DIRECTION
      USE GRID_PARAMETERS
      USE EXOBASE_BC
      USE SOLAR_LYMAN_ALPHA
      USE PHYSICS_TAG

   contains

      Subroutine Initialize_Setting
         USE SETTING, only: i_Photoionization

         call gen_points_for_NV
         call Init_Parameter
         call Get_exobaseBC
         call read_Lya_Bph  ;  if (i_Photoionization .eq. 0) then; bph = 0.d0; endif
         call Physics_tag
         
      end Subroutine
   
   END MODULE INTEGRATED_INITIALIZATION
   


   Module MPI_MATE
      integer nprocs, ierr
      integer :: rank
         real*8 :: rad, lon, lat   ! RANK dependent variables
         integer :: nR_loc, ilon, ilat, il
         integer :: nR_loc_MPI(nprocs)
   contains

      Subroutine Initialize_MPI
         call MPI_INIT(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

         call Calculate_Local_nRadial
      end Subroutine

      Subroutine Calculate_Local_nRadial

         use SETTING, only: nRadial
         IMPLICIT NONE

         nR_loc = nRadial / nprocs
         if (rank .lt. mod(nRadial, nprocs)) then
            nR_loc = nR_loc + 1
         endif
         nR_loc_MPI(rank+1) = nR_loc

      end Subroutine

   END MODULE MPI_MATE

 
   
   Module GRID_PARAMETERS

      USE SETTING
      IMPLICIT NONE

      real*8 radial_distance_range(nRadial), energy_range(nEnergy)
      real*8 longitude_range(nLong), latitude_range(nLat), latitudeNS_range(nLat_NS)
      real*8 radial_boundary(2)

   contains

      Subroutine Init_Parameter

         integer iR, iE, ilon, ilat
         integer index_Emin

         do iR=1,nRadial   ;  radial_distance_range(iR) = RadialRange_min + (iR-1)*dR         ;  enddo
         radial_distance_range = radial_distance_range * Re

         if (nEnergy .eq. 121) then  ! for 0.0025 - 10 eV
!            index_Emin = -33              ! Emin at 0.001 eV
               index_Emin = -20              ! Emin at 0.0025 eV
               do iE=index_Emin,index_Emin+nEnergy-1
                     energy_range(iE+1-index_Emin) = 0.01 * 1000.d0 ** (iE/100.d0)
               enddo
         else
               if (nEnergy .eq. 61) then  ! for 0.001 - 1 eV
                     do iE=1,nEnergy 
                           energy_range(iE) = 10.d0 ** ((iE-1)/20.d0 - 3.d0)
                     enddo
               else
                     print*, 'ERROR: Set a proper nEnergy'
               endif
         endif

         do ilon=1,nLong   ;  longitude_range(ilon)  = (ilon-1)*pi/180.d0*geores              ;  enddo
         do ilat=1,nLat    ;  latitude_range(ilat)   = (ilat-1)*pi/180.d0*geores              ;  enddo
         do ilat=1,nLat_NS ;  latitudeNS_range(ilat) = ((ilat-1.d0)*geores-90.d0) *pi/180.d0  ;  enddo

         radial_boundary(1) = inner_boundary
         radial_boundary(2) = outer_boundary

         return
      End

   END MODULE GRID_PARAMETERS



   Module SET_VELOCITY_DIRECTION
      USE SETTING
      IMPLICIT NONE

   contains

      Subroutine gen_points_for_NV

         IMPLICIT NONE

         integer :: nsize, n2
         integer :: i, j
         real*8 piset(nTheta), piset2(0:nTheta)

         do i=1,nTheta
            piset(i)=i*1.d0/(nTheta-1)*pi
            piset2(i)=i*2.d0/(nTheta-1)*pi
         enddo
         piset2(0)=0.d0;

         nsize=2
         do j=1,nTheta-2
            n2=nint((nTheta*2-2)*dsin(piset(j)))
            do i=0,n2-1
               nsize=nsize+1
            enddo
         enddo

         nvel = nsize

         return
      End


      Subroutine gen_points(vel_dir)
         ! Generates points on a sphere at n values of theta.
         ! Points are roughly evenly spaced on the sphere \
         ! n : number of theta values for points, odd n are preferred for most even spacing \
         !  output : boolean; if 0, returns points as normal;
         !           if 1, returns array specifying how many points are at each theta value"
         IMPLICIT NONE

         integer :: nsize, n2
         integer :: i, j, i0
         real*8 piset(nTheta), piset2(0:nTheta), piset2_i
         real*8 :: co1, co2, co3, co4
         real*8, dimension(nvel,3) :: vel_dir         
         real*8, dimension(:,:), allocatable :: coords
         real*8 :: t1,t2,t3

         do i=1,nTheta
            piset(i)=i*1.d0/(nTheta-1)*pi
            piset2(i)=i*2.d0/(nTheta-1)*pi
         enddo
         piset2(0)=0.d0;

         nsize=2
         do j=1,nTheta-2
            n2=nint((nTheta*2-2)*dsin(piset(j)))
            do i=0,n2-1
               nsize=nsize+1
            enddo
         enddo

         if (nsize .ne. nvel) then
            print*, "gen_points: <nvel> is not equal to <nsize>"
            print*, "nsize = ", nsize
            stop
         endif

         allocate(coords(nsize,3)); coords=0.d0

         i0=1+1 ! index for fortran
         do j=1,nTheta-2
            if (piset(j).eq.pi/2) then
               co1=1.d0; co2=0.d0
            else
               co1=dsin(piset(j)) ; co2=dcos(piset(j))
            endif

            n2=nint((nTheta*2-2)*co1)
            do i=0,n2-1
               piset2_i = 2.d0*i/n2*pi

               if (piset2_i.eq.pi/2 .or. piset2_i.eq.pi .or. piset2_i.eq.1.5*pi) then
                  co3=nint(dsin(piset2_i)); co4=nint(dcos(piset2_i))
               else
                  co3=dsin(piset2_i); co4=dcos(piset2_i)
               endif
               coords(i0,1) = co4*co1
               coords(i0,2) = co3*co1
               coords(i0,3) = co2

               i0=i0+1
         enddo
         enddo
         coords(1,3)=1.d0
         coords(nsize,3)=-1.d0

         t1=0.d0; t2=0.d0; t3=0.d0
         do i=1,nsize
            t1=t1+coords(i,1)
            t2=t2+coords(i,2)
            t3=t3+coords(i,3)
         enddo

         vel_dir = coords

         deallocate(coords)
  
         return
      End


      Subroutine gen_points_for_each_row(row)
         ! "gen_points" with output=1 in python code.
         IMPLICIT NONE

         integer :: i, j, k, n2
         integer row(0:nTheta)
         real*8 piset(nTheta), piset2(0:nTheta)

         row = 0
         do i=1,nTheta
            piset(i)=i*1.d0/(nTheta-1)*pi
            piset2(i)=i*2.d0/(nTheta-1)*pi
         enddo
         piset2(0)=0.d0

         row(0) = 1
         row(1) = 1
         k=2
         do j=1,nTheta-2
            n2=nint((nTheta*2-2)*dsin(piset(j)))
            do i=0,n2-1
               row(k) = row(k) + 1
            enddo
            k=k+1
         enddo
         row(k)=1

         return
      End



   END MODULE SET_VELOCITY_DIRECTION



   MODULE VOLUME_ELEMENT

      USE SETTING
      IMPLICIT NONE

      contains






   Module PHYSICS_TAG

      USE SETTING
      IMPLICIT NONE

      character*10 tag_phys
      public :: Physics_tag

   contains
      
      Subroutine Physics_tag
         ! Example: tag = "GRCPX" or "GRC"
         integer i
         character(len=1), dimension(n_physics) :: phy_name=''

         if (i_EarthGravity .eq. 1)           phy_name(1)='G'
         if (i_SolarRadiationPressure .eq. 1) phy_name(2)='R'
         if (i_CoriolisForce_GSE .eq. 1)      phy_name(3)='C'
         if (i_Photoionization .eq. 1)        phy_name(4)='P'
         if (i_ChargeExchange .eq. 1)         phy_name(5)='X'

         tag_phys = ''
         do i=1,n_physics
            tag_phys = trim(tag_phys) // trim(phy_name(i))
         enddo

         return
      End

   END MODULE PHYSICS_TAG


