!=======================================================================
!
! $Id: $
!
! FILE
!   setkin_fastj.h - lookup matrix for GMI vs FastJX cross-sections
!             (setkin_fastj.h)
!   Oct 2018 - SDS
!
! DESCRIPTION
!   Include file that provides ascii strings identifying
!   GMI vs FastJX photo reaction names
!
!  Input mechanism file:    StratTrop_Minimal.txt
!  Reaction dictionary:     GMI_reactions_JPL19.db
!  Setkin files generated:  Thu Nov 13 18:20:37 2025
!
!=======================================================================
!
!.... 
      integer idxph
!
!... Reference label for FastJX
      character*8 :: fastj_lookup(NUM_J)
      data (fastj_lookup(idxph), idxph=1,NUM_J) / &
     & "O2      ", "O3      ", "O3(1D)  ", "NO      ", "N2O     ",  & 
     & "NO2     ", "H2O2    ", "CH3OOH  ", "H2COa   ", "H2COb   ",  & 
     & "HNO3    ", "HNO2    ", "HNO4    ", "HNO4    ", "NO3     ",  & 
     & "NO3     ", "N2O5    ", "N2O5    ", "Cl2     ", "OClO    ",  & 
     & "Cl2O2   ", "HOCl    ", "ClNO3a  ", "ClNO3b  ", "BrCl    ",  & 
     & "BrO     ", "HOBr    ", "BrNO3   ", "BrNO3   ", "CH3Cl   ",  & 
     & "CCl4    ", "MeCCl3  ", "CFCl3   ", "CF2Cl2  ", "F113    ",  & 
     & "F141b   ", "F142b   ", "CH3Br   ", "H1301   ", "H1211   ",  & 
     & "ActAld  ", "ActAld  ", "PAN     ", "PrAld   ", "Acet-a  ",  & 
     & "CH3OOH  ", "CH3OOH  ", "CH3OOH  ", "CH3OOH  " /
!... Quantum yield for photolysis
      real*8 :: fastj_yield(NUM_J)
      data (fastj_yield(idxph), idxph=1,NUM_J) / &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,  &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,  &
     &   1.000D+00,   1.000D+00,   2.000D-01,   8.000D-01,   8.860D-01,  &
     &   1.140D-01,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,  &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,  &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,  &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,  &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,  &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,  &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00 /
!
!... end setkin_fastj.h
