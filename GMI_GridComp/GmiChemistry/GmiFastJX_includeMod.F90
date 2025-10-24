
      module GmiFastJX_includeMod

      private
      PUBLIC  :: t_fastJXbundle

# include "parm_MIE_fastJX65.h"

TYPE t_fastJXbundle
!... solar cycle capability

       REAL*8 :: fjx_solar_cycle_param(W_)
       REAL*8 :: lym_solar_cycle_param(5)

!----------------------------------------------------------------------------
       integer :: num_CCM_WL, num_CCM_aers, num_CCM_mom
       real*8, pointer :: CCM_WL(:)
       real*8, pointer :: CCM_SSALB(:,:,:,:,:)
       real*8, pointer :: CCM_OPTX (:,:,:,:,:)
       real*8, pointer :: CCM_SSLEG(:,:,:,:,:,:)

END TYPE t_fastJXbundle
      end module GmiFastJX_includeMod
