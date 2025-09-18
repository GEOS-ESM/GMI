!-------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadSolarCycle_mod
!
! !INTERFACE:
!
    module ReadSolarCycle_mod
!
   USE MAPL
!
    implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: readSolarCycleData
!
! !DESCRIPTION:
! Routine to set or read solar cycle scaling factor data.
!
! !REVISION HISTORY:
! Luke Oman 4 Aug 2016 Added Solar Cycle Read file
!
!EOP
!-------------------------------------------------------------------------

  CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readSolarCycleData
!
! !INTERFACE:
!
      subroutine readSolarCycleData ( s_cycle_dates,s_cycle,lym_cycle, sc_infile_name)
!
! !USES:
!
      use GmiASCIIoperations_mod,   only : AsciiOpenRead
!
      implicit none
!
#include "setkin_par.h"
#include "GmiParameters.h"
#include "parm_MIE_fastJX65.h"
      integer, parameter :: lym_bin = 5
!      integer, parameter :: hdr_lines = 2
      integer  :: hdr_lines
!
! !INPUT PARAMETERS:
      real, intent(inOut), dimension(2628)    :: s_cycle_dates        ! 2628 months : 1882 - 2100
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: sc_infile_name
      real, intent(inOut), dimension(W_,2628) :: s_cycle              ! 2628 months : 1882 - 2100
      real, intent(inOut), dimension(lym_bin,2628) :: lym_cycle             ! 2628 months : 1882 - 2100

! !DESCRIPTION:
! This routine sets/reads the solar cycle scaling values.
!
! !LOCAL VARIABLES:
      integer :: ic, im
      integer :: lun
      integer :: fjx_bin = W_
!
      character*100 :: hdr
!
!-----------------------------------------------------------------------------
!BOC

         call AsciiOpenRead (lun, sc_infile_name)

!         do ic = 1, 2628
!                  Read (lun, 900) s_cycle_dates(ic), (s_cycle(im,ic), im = 1, fjx_bin)
!         end do
!
! 900     format (19f10.5)
!
         IF(MAPL_AM_I_ROOT()) print *,'Reading solar cycle coefficients from : ',TRIM(sc_infile_name)
         hdr_lines = 0
         do ic = 1, 100
           Read (lun, *) hdr
!... move leading ' ' to be trailing
           hdr = ADJUSTL(hdr)
           if(SCAN(hdr,'#').ne.1) exit
           hdr_lines = hdr_lines+1
           IF(MAPL_AM_I_ROOT()) print *,hdr
         enddo
!... reset input file
         Close (lun)
         call AsciiOpenRead (lun, sc_infile_name)
!... skip header lines
         do ic = 1, hdr_lines
           Read (lun, 901) hdr
         enddo
!
!... read in coefficients
         do ic = 1, 5258/2-1
           Read (lun, 902) s_cycle_dates(ic), (s_cycle(im,ic), im = 1, fjx_bin)
           Read (lun, 903) (lym_cycle(im,ic), im = 1, lym_bin)
         end do

 901     format (a80)
 902     format (19f10.5)
 903     format (5f10.5)

         Close (lun)

      return

      end subroutine readSolarCycleData
!EOC
!---------------------------------------------------------------------------

    end module ReadSolarCycle_mod
