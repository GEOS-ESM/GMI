!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP      
! !MODULE: GmiUpdateChemistry_mod
!
      module GmiUpdateChemistry_mod
!     
! !USES:
      use GmiSolverInterface_mod,      ONLY : Update_Smv2chem
      use GmiArrayBundlePointer_mod,   ONLY : t_GmiArrayBundle
      USE GmiThermalRateConstants_mod, ONLY : Accum_Qqjk
      USE GmiSavedVariables_mod,       ONLY : t_ChemistrySaved
!    
      implicit none            
!     
      private                  
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: updateChemistry
!
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
! !AUTHOR:
!   Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
! !REVISION HISTORY:
!   - September 8, 2010 * Jules Kouatchou
!     GMI Chemistry basically became the solver only.
!EOP  
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateChemistry
!
! !INTERFACE:
!
      subroutine updateChemistry ( savedVars, rootProc, do_ftiming,            &
                       metdata_name_org, metdata_name_model, do_qqjk_inchem,   &
                       do_qqjk_reset, pr_qqjk, surfEmissForChem, press3c,      &
                       press3e, pr_smv2, pr_nc_period, mass, concentration,    &
                       qjgmi, qkgmi, kel, humidity, pctm2, qqjgmi, qqkgmi, yda,&
                       qqkda, qqjda, do_smv_reord, do_synoz,                   &
                       do_semiss_inchem, do_wetchem, nymd, nhms, gmi_sec, tdt, &
                       pr_diag, loc_proc, synoz_threshold, chem_cycle,         &
                       chem_mask_klo, chem_mask_khi, imgas_num, initrogen_num, &
                       ioxygen_num, isynoz_num, num_species, num_qks, num_qjs, &
                       num_qjo, num_sad, num_molefrac, num_chem, num_active,   &
                       ilong, ilat, ivert, itloop, ilo, ihi, julo, jhi, i1, i2,&
                       ju1, j2, k1, k2)
!
      implicit none

#     include "setkin_par.h"
#     include "gmi_AerDust_const.h"
!
! !INPUT PARAMETERS:
      character (len=*) ,intent(in) :: metdata_name_org   ! first  part of metdata_name, e.g., "NCAR"
      character (len=*) ,intent(in) :: metdata_name_model ! second part of metdata_name, e.g., "MATCH"
      logical, intent(in) :: pr_diag
      logical, intent(in) :: do_smv_reord, do_synoz, do_semiss_inchem
      logical, intent(in) :: do_wetchem
      logical, intent(in) :: rootProc, do_ftiming
      logical, intent(in) :: do_qqjk_inchem
      logical, intent(in) :: pr_qqjk
      logical, intent(in) :: pr_smv2
      integer, intent(in) :: loc_proc
      integer, intent(in) :: nymd, nhms
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: ilong, ilat, ivert, itloop
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_species, num_qks, num_qjs, num_qjo, num_sad
      integer, intent(in) :: num_chem, num_active, num_molefrac
      real*8 , intent(in) :: tdt
      integer, intent(in) :: imgas_num
      integer, intent(in) :: initrogen_num, ioxygen_num, isynoz_num
      integer, intent(in) :: chem_mask_klo, chem_mask_khi
      real*8 , intent(in) :: gmi_sec
      real*8,  intent(in) :: synoz_threshold
      real*8,  intent(in) :: chem_cycle
      real*8 , intent(in) :: pr_nc_period
      real*8 , intent(in) :: humidity(i1:i2,   ju1:j2,   k1:k2) 
      real*8 , intent(in) :: pctm2   (ilo:ihi, julo:jhi)  
      real*8 , intent(in) :: mass    (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in) :: kel    (ilo:ihi, julo:jhi, k1:k2) 
                             ! atmospheric pressure at the center of 
                             ! each grid box (mb)
      real*8 , intent(in) :: press3c(ilo:ihi, julo:jhi, k1:k2)
                             ! atmospheric pressure at the edge of 
                             ! each grid box (mb)
      real*8 , intent(in) :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      real*8,  intent(in)  :: surfEmissForChem(i1:i2, ju1:j2, num_species)
!
! !OUTPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_ChemistrySaved), intent(inOut) :: savedVars
      logical, intent(inOut) :: do_qqjk_reset
!      real*8 , intent(inOut) :: qqjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)
!      real*8 , intent(inOut) :: qqkgmi(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inOut) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inOut) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inOut) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
      type (t_GmiArrayBundle), intent(in)    :: qjgmi(num_qjo)
      type (t_GmiArrayBundle), intent(in)    :: qkgmi(num_qks)
      type (t_GmiArrayBundle), intent(inout) :: qqjgmi(num_qjo)
      type (t_GmiArrayBundle), intent(inout) :: qqkgmi(num_qks)
!
! !DESCRIPTION: 
! This routine updates the chemistry. It is mainly the solver only.
!
! !LOCAL VARIABLES:
      integer :: ix
      integer :: num_loops
      real*8  :: chemintv
      type (t_GmiArrayBundle) :: tm1Conc(num_species)
!
! !REVISION HISTORY:
!  Initial code.

!EOP
!-------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'updateChemistry called by ', loc_proc

      chemintv = tdt * chem_cycle

      if (chem_cycle < 1.0d0) then
         num_loops = Nint (1.0d0 / chem_cycle)
      else
         num_loops = 1
      end if
!
!... save initial concentration for calc of qqj/qqk diags
      if (pr_qqjk .and. .not. do_qqjk_inchem) then
        do ix = 1, num_species
          allocate ( tm1Conc(ix)%pArray3D (i1:i2, ju1:j2, k1:k2) )
          tm1Conc(ix)%pArray3D(:,:,:) = concentration(ix)%pArray3D(:,:,:)
        enddo
      end if

      do ix = 1, num_loops
         call Update_Smv2chem (savedVars, chemintv, surfEmissForChem,            &
                   humidity, qjgmi, qkgmi, press3e, pctm2, kel, concentration,   &
                   pr_diag, pr_qqjk, pr_smv2, do_smv_reord, do_synoz,            &
                   do_qqjk_inchem, do_semiss_inchem, imgas_num, initrogen_num,   &
                   ioxygen_num, isynoz_num, yda, qqkda, qqjda, pr_nc_period,     &
                   tdt, chem_mask_klo, chem_mask_khi, loc_proc, synoz_threshold, &
                   ilong, ilat, ivert, itloop, i1, i2, ju1, j2, k1, k2, ilo,     &
                   ihi, julo, jhi, num_molefrac,                                 &
                   num_qjo, num_qks, num_qjs, num_active, num_species, rootProc)
      end do

      if (pr_qqjk .and. .not. do_qqjk_inchem) then
        call Accum_Qqjk (do_qqjk_reset, imgas_num, tm1Conc,  &
                   qjgmi, qkgmi,  qqjgmi, qqkgmi, num_molefrac,     &
                   num_species, num_qks, num_qjs, num_qjo,          &
                   pr_diag, loc_proc, ilong, i1, i2, ju1, j2, k1, k2)
        do ix = 1, num_species
          deallocate ( tm1Conc(ix)%pArray3D )
        enddo
      end if

      return

      end subroutine updateChemistry
!EOC
!-------------------------------------------------------------------------
      end module GmiUpdateChemistry_mod

