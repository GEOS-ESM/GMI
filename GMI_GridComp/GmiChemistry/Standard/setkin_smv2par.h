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
!  Chemistry input file:    Standard_h2o.txt
!  Reaction dictionary:     GMI_reactions_JPL19.db
!  Setkin files generated:  Wed Aug 20 21:51:37 2025
!
!========1=========2=========3=========4=========5=========6=========7==

      integer &
     &  SK_IGAS &
     & ,SK_IPHOT &
     & ,SK_ITHERM &
     & ,SK_NACT

      parameter (SK_IGAS   = 123)
      parameter (SK_IPHOT  =  81)
      parameter (SK_ITHERM = 297)
      parameter (SK_NACT   = 119)

!                                  --^--

