
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   gmi_sad_constants.h
!
! DESCRIPTION
!   This include file contains the constants that map a particular type of
!   aerosol SAD (i.e., surface area density) to a numeric value for array
!   indexing.
!
!=============================================================================

! Indices for NSAD  (value in setkin_par.h):
      integer, parameter ::  &
        ILBSSAD  = 1,        & ! index for liquid binary sulfate       SADs
        ISTSSAD  = 2,        & ! index for supercooled ternary sulfate SADs
        INATSAD  = 3,        & ! index for "NAT" SADs
        IICESAD  = 4,        & ! index for ice   SADs
        IPYROSAD = 5           ! index for Pyro  SADs

! Indices for NSADdust  (value in setkin_par.h):
!  ????

! Indices for NSADaer   (value in setkin_par.h):
      integer, parameter ::  &
        ISO4SAD   = 1,       & ! index for sulfate           SADs
        IBCSAD    = 2,       & ! index for black carbon      SADs
        IOCSAD    = 3,       & ! index for organic carbon    SADs
        ISSaccSAD = 4,       & ! index for sea salt (accum)  SADs
        ISScrsSAD = 5,       & ! index for sea salt (coarse) SADs
        ISO4vSAD  = 6          ! index for volcanic sulfate  SADs

