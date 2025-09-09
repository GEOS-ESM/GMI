!=======================================================================
!
! $Id: setkin_ibcb.h $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_ibcb.h
!
! DESCRIPTION
!   This include file contains information about treatment of surface
!   boundary conditions.
!
!  Input mechanism:        GeosCCM_Combo_Minimal2_Mechanism.txt
!  Reaction dictionary:    GMI_reactions_JPL19.db
!  Setkin files generated: Wed Mar  5 20:39:38 2025
!
!=======================================================================
!
!.... Set default boundary condition types
!
!.... Type 1 means fixed concentration
!.... surface boundary condition, Type 2 means
!.... fixed flux surface boundary condition
!
      ibcb(:)           = 2

      ibcb(NACT+1:IGAS) = 1
!
!.... Reset boundary condition type for special cases
!                C2H6
!      ibcb(10) = 1
!                C3H8
!      ibcb(11) = 1
!                CCl4
!      ibcb(12) = 1
!                CF2ClBr
!      ibcb(13) = 1
!                CF3Br
!      ibcb(14) = 1
!                CFC113
!      ibcb(15) = 1
!                CFC11
!      ibcb(16) = 1
!                CFC12
!      ibcb(17) = 1
!                CH3Br
!      ibcb(19) = 1
!                CH3CCl3
!      ibcb(20) = 1
!                CH3Cl
!      ibcb(21) = 1
!                CH4
!      ibcb(22) = 1
!                CO
!      ibcb(28) = 1
!                H2
!      ibcb(31) = 1
!                H2O
!      ibcb(33) = 1
!                HCFC141b
!      ibcb(35) = 1
!                HCFC142b
!      ibcb(36) = 1
!                HCFC22
!      ibcb(37) = 1
!                HCl
!      ibcb(38) = 1
!                ISOP
!      ibcb(48) = 1
!                N2O
!      ibcb(55) = 1
!                                  --^--
