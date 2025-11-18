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
!  Input mechanism:        StratTrop_JPL19.txt
!  Reaction dictionary:    GMI_reactions_JPL19.db
!  Setkin files generated: Thu Nov 13 18:18:38 2025
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
!                C2Cl4
!      ibcb(11) = 1
!                C2H4Cl2
!      ibcb(12) = 1
!                C2H6
!      ibcb(13) = 1
!                C3H8
!      ibcb(14) = 1
!                CCl4
!      ibcb(15) = 1
!                CF2Br2
!      ibcb(16) = 1
!                CF2ClBr
!      ibcb(17) = 1
!                CF3Br
!      ibcb(18) = 1
!                CFC113
!      ibcb(19) = 1
!                CFC114
!      ibcb(20) = 1
!                CFC115
!      ibcb(21) = 1
!                CFC11
!      ibcb(22) = 1
!                CFC12
!      ibcb(23) = 1
!                CHCl3
!      ibcb(24) = 1
!                CH2Cl2
!      ibcb(26) = 1
!                CH3Br
!      ibcb(28) = 1
!                CH3CCl3
!      ibcb(29) = 1
!                CH3Cl
!      ibcb(30) = 1
!                CH4
!      ibcb(31) = 1
!                CO
!      ibcb(38) = 1
!                H2
!      ibcb(46) = 1
!                H2402
!      ibcb(47) = 1
!                H2O
!      ibcb(49) = 1
!                HCFC141b
!      ibcb(52) = 1
!                HCFC142b
!      ibcb(53) = 1
!                HCFC22
!      ibcb(54) = 1
!                HCl
!      ibcb(55) = 1
!                ISOP
!      ibcb(73) = 1
!                N2O
!      ibcb(88) = 1
!                                  --^--
