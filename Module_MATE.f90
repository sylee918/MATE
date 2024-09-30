      Module Module_for_NVelocityDirection
         integer :: N_vel_directions
      contains

         Subroutine gen_points_for_NV()

            include "Setting.inc"
            integer :: nsize, n2
            integer :: i, j, i0
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

            N_vel_directions = nsize

            return
         End

      End Module


      Module Module_Physics_tag
         character*10 tag_phys
      contains
         
         Subroutine Physics_tag()
            ! Generate "tag_phys".
            ! Example: tag = "GRCPX" or "GRC"
            include "Setting.inc"
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

      End Module