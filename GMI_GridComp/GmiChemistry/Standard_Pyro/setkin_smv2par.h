!=======================================================================
!
! $Id: $
!
! FILE
!   setkin_smv2par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    Standard_Pyro.txt
!  Reaction dictionary:     GMI_reactions_JPL19.db
!  Setkin files generated:  Wed Jul  8 20:30:04 2026
!
!========1=========2=========3=========4=========5=========6=========7==

      integer &
     &  SK_IGAS &
     & ,SK_IPHOT &
     & ,SK_ITHERM &
     & ,SK_NACT

      parameter (SK_IGAS   = 125)
      parameter (SK_IPHOT  =  81)
      parameter (SK_ITHERM = 307)
      parameter (SK_NACT   = 121)

!                                  --^--

