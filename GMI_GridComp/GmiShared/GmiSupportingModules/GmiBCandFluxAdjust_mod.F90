#include "MAPL_Generic.h"

module GmiBCandFluxAdjust_mod

    use ESMF
    use MAPL
    use GmiSpeciesRegistry_mod,        only : getSpeciesIndex

    implicit none

    INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

    private
    public  :: GetBrAdjustments

!=============================================================================
!
! CODE DEVELOPERS
!   Steve Steenrod  NASA/GSFC Code 614
!   Mike Manyin     NASA/GSFC Code 614
!
!=============================================================================

     contains

!-----------------------------------------------------------------------------
!
! ROUTINE
!   GetBrAdjustments
!
! DESCRIPTION
!   This routine ...
!
! ARGUMENTS
!
!-----------------------------------------------------------------------------

      subroutine GetBrAdjustments           (  &
           add_ch3br_opt,   add_ch2br2_opt,    &
          mult_ch3br_opt,  mult_ch2br2_opt,    &
          rc                                )

!     ----------------------
!     Argument declarations.
!     ----------------------

      real(kind=DBL), optional, intent(out)  ::   add_ch3br_opt,   add_ch2br2_opt  ! Suitable for Forced Boundary Conditions
      real(kind=DBL), optional, intent(out)  ::  mult_ch3br_opt,  mult_ch2br2_opt  ! Suitable for Flux   Boundary Conditions
      integer,        optional, intent(out)  ::  rc                                ! Error return code


!     ----------------------
!     Variable declarations.
!     ----------------------

      character(len=*), parameter :: IAm = 'GetBrAdjustments'

      integer        :: ich3br
      integer        :: ichbr3
      integer        :: ich2br2
      integer        :: ichbr2cl
      integer        :: ichbrcl2
      integer        :: ich2brcl

      real(kind=DBL) ::     add_ch3br,  add_ch2br2
      real(kind=DBL) ::    mult_ch3br, mult_ch2br2


!     ----------
!     Initialize
!     ----------
      add_ch3br    = 0.0d0
      add_ch2br2   = 0.0d0
      mult_ch3br   = 1.0d0
      mult_ch2br2  = 1.0d0


!... Because we have mechanisms with incomplete VSLBr species
!... we need to compensate in other VSLBr species, for now CH3Br and CH2Br2
      ich3br   = getSpeciesIndex('CH3Br',  .TRUE.)   ! .TRUE. arg=>don't error exit if species not there
      ichbr3   = getSpeciesIndex('CHBr3',  .TRUE.)
      ich2br2  = getSpeciesIndex('CH2Br2', .TRUE.)
!... other VSLBr species that could be included in a future mechanism
      ichbr2cl = getSpeciesIndex('CHBr2Cl',.TRUE.)
      ichbrcl2 = getSpeciesIndex('CHBrCl2',.TRUE.)
      ich2brcl = getSpeciesIndex('CH2BrCl',.TRUE.)




!... Because we have mechanisms with incomplete VSLBr chemistry
!...  we need to compensate in the other existing VSLBr species if they are fixed BCs
!... species importance in order:
!... (assume that species will be added in this order)
!...  CH3Br   - unless have CHBr3 and CH2Br2 in mech, add extra 5.0 ppt
!...  CHBr3   - subtract 3.6 from CH3Br's extra 5 ppt
!...  CH2Br2  - make CH3Br extra 0 ppt and add extra 1 ppt to CH2Br2
!...  CHBr2Cl - subtract 0.6 from CH2Br2 extra 1 ppt
!...  CHBrCl2 - subtract 0.3 from CH2Br2 extra 1 ppt
!...  CH2BrCl - subtract 0.1 from CH2Br2 extra 1 ppt
      if(ich3br .gt. 0) then
         add_ch3br = 5.0d-12
        mult_ch3br = -999.0 * 1.0d0     ! need value for this
        if(ichbr3 .gt. 0) then
           add_ch3br = add_ch3br-3.6d-12
          mult_ch3br = -999.0 * 1.0d0   ! need value for this
          if(ich2br2 .gt. 0) then
             add_ch3br = 0.0d-12
            mult_ch3br = 1.0d0
             add_ch2br2 = 1.0d-12   ! NOTE: Only used if CH2Br2 is a Forced Boundary Condition
            mult_ch2br2 = 1.8d0     ! NOTE: Only used if CH2Br2 is an Emission  [Qing email 4/3/24]
!... other VSLBr species that could be included (CHBr2Cl, CHBrCl2, CH2BrCl)
              if(ichbr2cl .gt. 0) then
                PRINT*,TRIM(Iam)//' Running a new section of code, please review it'
                PRINT*, " "
                VERIFY_(101)
                 add_ch2br2 = add_ch2br2-0.6d-12
                mult_ch2br2 = -999.0 * 1.0d0   ! need value for this
                if(ichbrcl2 .gt. 0) then
                  PRINT*,TRIM(Iam)//' Running a new section of code, please review it'
                  PRINT*, " "
                  VERIFY_(102)
                   add_ch2br2 = add_ch2br2-0.3d-12
                  mult_ch2br2 = -999.0 * 1.0d0   ! need value for this
                  if(ich2brcl .gt. 0) then
                    PRINT*,TRIM(Iam)//' Running a new section of code, please review it'
                    PRINT*, " "
                    VERIFY_(103)
                     add_ch2br2 = add_ch2br2-0.1d-12
                    mult_ch2br2 = -999.0 * 1.0d0   ! need value for this
                  endif
                endif
              endif
          endif
        endif
      endif

!     ---------------------------
!     Fill the Optional Arguments
!     ---------------------------
      if ( present(add_ch3br_opt  ) )  add_ch3br_opt   = add_ch3br
      if ( present(add_ch2br2_opt ) )  add_ch2br2_opt  = add_ch2br2
      if ( present(mult_ch3br_opt ) )  mult_ch3br_opt  = mult_ch3br
      if ( present(mult_ch2br2_opt) )  mult_ch2br2_opt = mult_ch2br2

      if ( present(rc) ) rc = ESMF_SUCCESS
      return

      end subroutine GetBrAdjustments

end module GmiBCandFluxAdjust_mod
