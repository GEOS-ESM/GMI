!=======================================================================
!
! $Id: setkin_surf_bc.h $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_surf_bc.h
!
! DESCRIPTION
!   This include file contains information about surface boundary
!   conditions.
!
!  Input mechanism:        GeosCCM_Combo_2020_HFC_S_VSLCL_JPL19.txt
!  Reaction dictionary:    GMI_reactions_JPL19.db
!  Setkin files generated: Tue Oct 22 16:51:18 2024
!
!=======================================================================

!     -----------------------
!     Parameter declarations.
!     -----------------------

!     -------------------------------------------------------
!     NUM_SBC : number of surface bounday conditions to reset
!     K_SBC   : max k to which species will be adjusted
!     -------------------------------------------------------

      integer, parameter :: NUM_SBC = 1

      integer, parameter :: K_SBC   = 2

!     ------------------
!     Integer variables.
!     ------------------

      integer, save :: sbc_map(NUM_SBC) = &
     &  (/ 0 /)

!                                  --^--

