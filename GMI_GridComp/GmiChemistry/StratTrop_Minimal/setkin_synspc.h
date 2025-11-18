!=======================================================================
!
! $Id: setkin_synspc.h $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_synspc.h
!
! DESCRIPTION
!   This include file contains information about synthetic stratospheric
!   sources for ozone and nitrogen oxides.
!
!  Input mechanism:        StratTrop_Minimal.txt
!  Reaction dictionary:    GMI_reactions_JPL19.db
!  Setkin files generated: Thu Nov 13 18:20:37 2025
!
!=======================================================================


      logical, parameter :: USE_SYNOZ = .true.
      logical, parameter :: USE_NODOZ = .true.

      integer, parameter :: MAXNODOZ_ELEM = 10

      integer, parameter :: NOX_ELEM_MAP(MAXNODOZ_ELEM) = &
     &  (/ INO, INO2, INO3, IN2O5, IN2O5,  0,  0,  0,  0,  0 /)

      integer, parameter :: NOY_ELEM_MAP(MAXNODOZ_ELEM) = &
     &  (/  IHNO2,  IHNO3,  IHNO4,  0,  0,  0,  0,  0,  0,  0 /)

!                                  --^--

