!=======================================================================
!
! $Id: $
!
! ROUTINE
!   kcalc - GMIMOD 3D model (setkin_kcalc.F90)
!   1 AUG 02 - PSC
!
! DESCRIPTION
!   Calculates and returns rate constants for thermal
!   reactions for the temperatures and pressures supplied
!
! ARGUMENTS
!  INPUT
!   pressure    : mb - profile of pressures
!   temperature : K - profile of temperatures
!   adcol       : molecules cm-3 - profile of total number density
!   specarr     : molecules cm-3 - profiles of species concentrations
!  OUTPUT
!   rcarr       : cm3 molecule-1 s-1 - rate constant values in units as
!                 appropriate
!
!  Input mechanism:        StratTrop_Minimal.txt
!  Reaction dictionary:    GMI_reactions_JPL19.db
!  Setkin files generated: Thu Nov 13 18:20:37 2025
!
!=======================================================================
      subroutine kcalc( npres0,sadcol,sadcol2,pressure,ptrop,cPBLcol, &
     &  temperature,fcld,lwc,adcol,specarr,rcarr,radA,FRH)

      implicit none

#     include "gmi_phys_constants.h"
#     include "setkin_par.h"
#     include "setkin_mw.h"
#     include "setkin_depos.h"
#     include "gmi_sad_constants.h"

!.... Argument declarations

      INTEGER, INTENT (IN)  :: npres0
      INTEGER, INTENT (IN)  :: cPBLcol(npres0)

      REAL*8,  INTENT (IN)  :: ptrop

      REAL*8,  INTENT (IN)  :: adcol       (       npres0)
      REAL*8,  INTENT (IN)  :: fcld        (       npres0)
      REAL*8,  INTENT (IN)  :: lwc         (       npres0)
      REAL*8,  INTENT (IN)  :: pressure    (       npres0)
      REAL*8,  INTENT (IN)  :: temperature (       npres0)
      REAL*8,  INTENT (IN)  :: sadcol      (NSAD  ,npres0)
!      REAL*8,  INTENT (IN)  :: sadreff     (NSAD  ,npres0)
      REAL*8,  INTENT (IN)  :: sadcol2     (NSADaer+NSADdust,npres0)
      REAL*8,  INTENT (IN)  :: radA        (NSADaer+NSADdust,npres0)
      REAL*8,  INTENT (IN)  :: specarr     (NMF   ,npres0)
      REAL*8,  INTENT (IN)  :: FRH         (       npres0)

      REAL*8,  INTENT (OUT) :: rcarr       (NUM_K ,npres0)

!.... Local variable declarations

      INTEGER               :: naltmax

      real*8 &
     &  nitrogen (npres0) &
     & ,oxygen   (npres0) &
     & ,water    (npres0)

! The sad_* variables are not used in the Tropospheric Mechanism.
      real*8 &
     &  sad_ice  (npres0) &
     & ,sad_lbs  (npres0) &
     & ,sad_nat  (npres0) &
     & ,sad_pyro (npres0) &
     & ,sad_sts  (npres0)

      real*8 mw(NSP), gammas_code
!
      real*8, DIMENSION (npres0) :: wt_h2so4, g_clono2, g_clono2_hcl, g_clono2_h2o, g_hocl_hcl
!... effective radii of stratospheric aerosols
      real*8 reff_lbs   ! reff_sts, reff_nat, reff_ice, reff_pyro
!
      mw(:) = mw_data(:)

      naltmax     = npres0
!===================================================
! Variables needed for gas/heterogeneous chemistry.
! adcol = air density (molec/cm3)
! FRH = relative humidity fraction (0-1)
! radA = effective radius of aerosol (cm)
! sadcol = surface area of aerosol/volume of air (cm2/cm3)
!===================================================

!....          Define molecular nitrogen and oxygen number densities
!
      nitrogen(:) = adcol(:) * MXRN2
      oxygen(:)   = adcol(:) * MXRO2
      water(:)    = specarr(33 ,:)
!... * 0.0d0
      reff_lbs  = 0.221d-4
!     reff_sts  = 0.221d-4
!     reff_nat  = 0.221d-4
!     reff_ice  = 0.221d-4
!     reff_pyro = 0.221d-4
!
      sad_lbs(:)  = sadcol(ILBSSAD, :)
      sad_sts(:)  = sadcol(ISTSSAD, :)
      sad_nat(:)  = sadcol(INATSAD, :)
      sad_ice(:)  = sadcol(IICESAD, :)
      sad_pyro(:) = sadcol(IPYROSAD,:)
!... Use GOCART pyro (BC+OC) for this reaction
!.old      sad_pyro(:) = sadcol2(NSADdust+2,:)+sadcol2(NSADdust+3,:)
!
!
!....          Start thermal rate constants
!
!....           O + O2 = O3
!
      rcarr(1,:) = skterlp(  6.000D-34 ,2.40D+00 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           O + O3 = 2 O2
!
      rcarr(2,:) = skarr(  8.000D-12 ,2060.0D+00 ,temperature)
!
!....           N2 + O1D = N2 + O
!
      rcarr(3,:) = skarr(  2.150D-11 ,-110.0D+00 ,temperature)
!
!....           O1D + O2 = O + O2
!
      rcarr(4,:) = skarr(  3.300D-11 ,-55.0D+00 ,temperature)
!
!....           O1D + O3 = 2 O2
!
      rcarr(5,:) = skarr(  1.200D-10 ,0.0D+00 ,temperature)
!
!....           O1D + O3 = 2 O + O2
!
      rcarr(6,:) = skarr(  1.200D-10 ,0.0D+00 ,temperature)
!
!....           H2O + O1D = 2 OH
!
      rcarr(7,:) = skarr(  1.630D-10 ,-60.0D+00 ,temperature)
!
!....           H2 + O1D = H + OH
!
      rcarr(8,:) = skarr(  1.200D-10 ,0.0D+00 ,temperature)
!
!....           N2O + O1D = N2 + O2
!
      rcarr(9,:) = skarr(  4.640D-11 ,-20.0D+00 ,temperature)
!
!....           N2O + O1D = 2 NO
!
      rcarr(10,:) = skarr(  7.260D-11 ,-20.0D+00 ,temperature)
!
!....           CH4 + O1D = MO2 + OH
!
      rcarr(11,:) = skarr(  1.310D-10 ,0.0D+00 ,temperature)
!
!....           CH4 + O1D = CH2O + H + HO2
!
      rcarr(12,:) = skarr(  3.500D-11 ,0.0D+00 ,temperature)
!
!....           CH4 + O1D = CH2O + H2
!
      rcarr(13,:) = skarr(  8.750D-12 ,0.0D+00 ,temperature)
!
!....           CFC12 + O1D = 2 Cl
!
      rcarr(14,:) = skarr(  1.400D-10 ,-25.0D+00 ,temperature)
!
!....           CFC113 + O1D = 3 Cl
!
      rcarr(15,:) = skarr(  2.320D-10 ,0.0D+00 ,temperature)
!
!....           HCFC22 + O1D = Cl
!
      rcarr(16,:) = skarr(  1.020D-10 ,0.0D+00 ,temperature)
!
!....           HCFC141b + O1D = 2 Cl
!
      rcarr(17,:) = skarr(  2.600D-10 ,0.0D+00 ,temperature)
!
!....           HCFC142b + O1D = Cl
!
      rcarr(18,:) = skarr(  2.000D-10 ,0.0D+00 ,temperature)
!
!....           H + O2 = HO2
!
      rcarr(19,:) = sktroe(  5.300D-32 ,1.80D0 & 
     &                     , 9.500D-11 ,-0.40D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           H + O3 = O2 + OH
!
      rcarr(20,:) = skarr(  1.400D-10 ,470.0D+00 ,temperature)
!
!....           O + OH = H + O2
!
      rcarr(21,:) = skarr(  1.800D-11 ,-180.0D+00 ,temperature)
!
!....           HO2 + O = O2 + OH
!
      rcarr(22,:) = skarr(  3.000D-11 ,-200.0D+00 ,temperature)
!
!....           H + HO2 = 2 OH
!
      rcarr(23,:) = skarr(  7.200D-11 ,0.0D+00 ,temperature)
!
!....           NO + O3 = NO2 + O2
!
      rcarr(24,:) = skarr(  3.000D-12 ,1500.0D+00 ,temperature)
!
!....           O3 + OH = HO2 + O2
!
      rcarr(25,:) = skarr(  1.700D-12 ,940.0D+00 ,temperature)
!
!....           HO2 + O3 = 2 O2 + OH
!
      rcarr(26,:) = skarr(  1.000D-14 ,490.0D+00 ,temperature)
!
!....           NO2 + O3 = NO3 + O2
!
      rcarr(27,:) = skarr(  1.200D-13 ,2450.0D+00 ,temperature)
!
!....           OH + OH = H2O + O
!
      rcarr(28,:) = skarr(  1.800D-12 ,0.0D+00 ,temperature)
!
!....           OH + OH = H2O2
!
      rcarr(29,:) = sktroe(  6.900D-31 ,1.00D0 & 
     &                     , 2.600D-11 ,0.00D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HO2 + OH = H2O + O2
!
      rcarr(30,:) = skarr(  4.800D-11 ,-250.0D+00 ,temperature)
!
!....           H2O2 + OH = H2O + HO2
!
      rcarr(31,:) = skarr(  1.800D-12 ,0.0D+00 ,temperature)
!
!....           HO2 + NO = NO2 + OH
!
      rcarr(32,:) = skarr(  3.440D-12 ,-260.0D+00 ,temperature)
!
!....           HO2 + HO2 = H2O2 + O2
!
      rcarr(33,:) = skho2dis (temperature ,adcol)
!
!....           H2O + HO2 + HO2 = H2O + H2O2 + O2
!
      rcarr(34,:) = skho2h2o (temperature ,adcol)
!
!....           H2 + OH = H + H2O
!
      rcarr(35,:) = skarr(  2.800D-12 ,1800.0D+00 ,temperature)
!
!....           CO + OH = H
!
      rcarr(36,:) = skohco (temperature ,adcol)
!
!....           CH4 + OH = H2O + MO2
!
      rcarr(37,:) = skohch4 (temperature)
!
!....           ClO + MO2 = CH2O + Cl + HO2 + O2
!
      rcarr(38,:) = skarr(  1.800D-11 ,600.0D+00 ,temperature)
!
!....           MO2 + NO = CH2O + HO2 + NO2
!
      rcarr(39,:) = skarr(  2.800D-12 ,-300.0D+00 ,temperature)
!
!....           HO2 + MO2 = MP + O2
!
      rcarr(40,:) = skarr(  4.100D-13 ,-750.0D+00 ,temperature)
!
!....           MO2 + MO2 = CH2O + O2
!
      rcarr(41,:) = skmo2dis_1 (temperature)
!
!....           MO2 + MO2 = 2 CH2O + 2 HO2
!
      rcarr(42,:) = skmo2dis_2 (temperature)
!
!....           MP + OH = H2O + MO2
!
      rcarr(43,:) = skarr(  2.660D-12 ,-200.0D+00 ,temperature)
!
!....           MP + OH = CH2O + H2O + OH
!
      rcarr(44,:) = skarr(  1.140D-12 ,-200.0D+00 ,temperature)
!
!....           CH2O + OH = CO + H2O + HO2
!
      rcarr(45,:) = skarr(  5.500D-12 ,-125.0D+00 ,temperature)
!
!....           N + O2 = NO + O
!
      rcarr(46,:) = skarr(  3.300D-12 ,3150.0D+00 ,temperature)
!
!....           N + NO = N2 + O
!
      rcarr(47,:) = skarr(  2.100D-11 ,-100.0D+00 ,temperature)
!
!....           NO2 + O = NO + O2
!
      rcarr(48,:) = skono2_d (temperature ,adcol)
!
!....           NO3 + O = NO2 + O2
!
      rcarr(49,:) = skarr(  1.300D-11 ,0.0D+00 ,temperature)
!
!....           NO2 + OH = HNO3
!
      rcarr(50,:) = sktroe(  1.800D-30 ,3.00D0 & 
     &                     , 2.800D-11 ,0.00D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HNO3 + OH = H2O + NO3
!
      rcarr(51,:) = skohhno3_j19 (temperature ,adcol)
!
!....           NO + OH = HNO2
!
      rcarr(52,:) = sktroe(  7.100D-31 ,2.60D0 & 
     &                     , 3.600D-11 ,0.10D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HNO2 + OH = H2O + NO2
!
      rcarr(53,:) = skarr(  3.000D-12 ,-250.0D+00 ,temperature)
!
!....           HO2 + NO2 = HNO4
!
      rcarr(54,:) = sktroe(  1.900D-31 ,3.40D0 & 
     &                     , 4.000D-12 ,0.30D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HNO4 = HO2 + NO2
!
      rcarr(55,:) = sktroe(  9.050D-05 ,3.4D0 &
     &                     , 1.900D+15 ,0.30D0 ,10900.0D+00 &
     &                     ,temperature ,adcol)
!
!....           HNO4 + OH = H2O + NO2 + O2
!
      rcarr(56,:) = skarr(  4.500D-13 ,-610.0D+00 ,temperature)
!
!....           HO2 + NO3 = NO2 + O2 + OH
!
      rcarr(57,:) = skarr(  3.500D-12 ,0.0D+00 ,temperature)
!
!....           NO + NO3 = 2 NO2
!
      rcarr(58,:) = skarr(  1.700D-11 ,-125.0D+00 ,temperature)
!
!....           NO3 + OH = HO2 + NO2
!
      rcarr(59,:) = skarr(  2.000D-11 ,0.0D+00 ,temperature)
!
!....           NO2 + NO3 = N2O5
!
      rcarr(60,:) = sktroe(  2.400D-30 ,3.00D0 & 
     &                     , 1.600D-12 ,-0.10D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           N2O5 = NO2 + NO3
!
      rcarr(61,:) = sktroe(  4.140D-04 ,3.0D0 &
     &                     , 2.760D+14 ,-0.10D0 ,10840.0D+00 &
     &                     ,temperature ,adcol)
!
!....           HCOOH + OH = H2O + HO2
!
      rcarr(62,:) = skarr(  4.000D-13 ,0.0D+00 ,temperature)
!
!....           NO2 + NO3 = NO + NO2 + O2
!
      rcarr(63,:) = skarr(  4.350D-14 ,1335.0D+00 ,temperature)
!
!....           CH2O + NO3 = CO + HNO3 + HO2
!
      rcarr(64,:) = skarr(  5.800D-16 ,0.0D+00 ,temperature)
!
!....           Cl + O3 = ClO + O2
!
      rcarr(65,:) = skarr(  2.300D-11 ,200.0D+00 ,temperature)
!
!....           Cl + H2 = H + HCl
!
      rcarr(66,:) = skarr(  3.050D-11 ,2270.0D+00 ,temperature)
!
!....           Cl + H2O2 = HCl + HO2
!
      rcarr(67,:) = skarr(  1.100D-11 ,980.0D+00 ,temperature)
!
!....           Cl + HO2 = HCl + O2
!
      rcarr(68,:) = skarr(  1.400D-11 ,-270.0D+00 ,temperature)
!
!....           Cl + HO2 = ClO + OH
!
      rcarr(69,:) = skarr(  3.600D-11 ,375.0D+00 ,temperature)
!
!....           ClO + O = Cl + O2
!
      rcarr(70,:) = skarr(  2.800D-11 ,-85.0D+00 ,temperature)
!
!....           ClO + OH = Cl + HO2
!
      rcarr(71,:) = skarr(  7.400D-12 ,-270.0D+00 ,temperature)
!
!....           ClO + OH = HCl + O2
!
      rcarr(72,:) = skarr(  6.000D-13 ,-230.0D+00 ,temperature)
!
!....           ClO + HO2 = HOCl + O2
!
      rcarr(73,:) = skarr(  2.600D-12 ,-290.0D+00 ,temperature)
!
!....           ClO + NO = Cl + NO2
!
      rcarr(74,:) = skarr(  6.400D-12 ,-290.0D+00 ,temperature)
!
!....           ClO + NO2 = ClONO2
!
      rcarr(75,:) = sktroe(  1.800D-31 ,3.40D0 & 
     &                     , 1.500D-11 ,1.90D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           ClO + ClO = 2 Cl + O2
!
      rcarr(76,:) = skarr(  3.000D-11 ,2450.0D+00 ,temperature)
!
!....           ClO + ClO = Cl2 + O2
!
      rcarr(77,:) = skarr(  1.000D-12 ,1590.0D+00 ,temperature)
!
!....           ClO + ClO = Cl + OClO
!
      rcarr(78,:) = skarr(  3.500D-13 ,1370.0D+00 ,temperature)
!
!....           ClO + ClO = Cl2O2
!
      rcarr(79,:) = sktroe(  1.900D-32 ,3.60D0 & 
     &                     , 3.700D-12 ,1.60D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           Cl2O2 = 2 ClO
!
      rcarr(80,:) = sktroe(  8.800D-06 ,3.6D0 &
     &                     , 1.710D+15 ,1.60D0 ,8537.0D+00 &
     &                     ,temperature ,adcol)
!
!....           HCl + OH = Cl + H2O
!
      rcarr(81,:) = skarr(  1.800D-12 ,250.0D+00 ,temperature)
!
!....           HOCl + OH = ClO + H2O
!
      rcarr(82,:) = skarr(  3.000D-12 ,500.0D+00 ,temperature)
!
!....           ClONO2 + O = ClO + NO3
!
      rcarr(83,:) = skarr(  3.600D-12 ,840.0D+00 ,temperature)
!
!....           ClONO2 + OH = HOCl + NO3
!
      rcarr(84,:) = skarr(  1.200D-12 ,330.0D+00 ,temperature)
!
!....           Cl + ClONO2 = Cl2 + NO3
!
      rcarr(85,:) = skarr(  6.500D-12 ,-135.0D+00 ,temperature)
!
!....           Br + O3 = BrO + O2
!
      rcarr(86,:) = skarr(  1.600D-11 ,780.0D+00 ,temperature)
!
!....           Br + HO2 = HBr + O2
!
      rcarr(87,:) = skarr(  4.800D-12 ,310.0D+00 ,temperature)
!
!....           Br + CH2O = CO + HBr + HO2
!
      rcarr(88,:) = skarr(  1.700D-11 ,800.0D+00 ,temperature)
!
!....           BrO + O = Br + O2
!
      rcarr(89,:) = skarr(  1.900D-11 ,-230.0D+00 ,temperature)
!
!....           BrO + HO2 = HOBr + O2
!
      rcarr(90,:) = skarr(  4.500D-12 ,-460.0D+00 ,temperature)
!
!....           BrO + NO = Br + NO2
!
      rcarr(91,:) = skarr(  8.800D-12 ,-260.0D+00 ,temperature)
!
!....           BrO + NO2 = BrONO2
!
      rcarr(92,:) = sktroe(  5.500D-31 ,3.10D0 & 
     &                     , 6.600D-12 ,2.90D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           BrO + ClO = Br + OClO
!
      rcarr(93,:) = skarr(  9.500D-13 ,-550.0D+00 ,temperature)
!
!....           BrO + ClO = Br + Cl + O2
!
      rcarr(94,:) = skarr(  2.300D-12 ,-260.0D+00 ,temperature)
!
!....           BrO + ClO = BrCl + O2
!
      rcarr(95,:) = skarr(  4.100D-13 ,-290.0D+00 ,temperature)
!
!....           BrO + BrO = 2 Br + O2
!
      rcarr(96,:) = skbrodis (temperature)
!
!....           HBr + OH = Br + H2O
!
      rcarr(97,:) = skarr(  5.500D-12 ,-200.0D+00 ,temperature)
!
!....           CH2O + O = CO + HO2 + OH
!
      rcarr(98,:) = skarr(  3.400D-11 ,1600.0D+00 ,temperature)
!
!....           CH4 + Cl = HCl + MO2
!
      rcarr(99,:) = skarr(  7.100D-12 ,1270.0D+00 ,temperature)
!
!....           CH2O + Cl = CO + HCl + HO2
!
      rcarr(100,:) = skarr(  8.100D-11 ,30.0D+00 ,temperature)
!
!....           CH3Cl + OH = Cl + H2O + HO2
!
      rcarr(101,:) = skarr(  1.960D-12 ,1200.0D+00 ,temperature)
!
!....           CH3CCl3 + OH = 3 Cl + H2O
!
      rcarr(102,:) = skarr(  1.640D-12 ,1520.0D+00 ,temperature)
!
!....           HCFC22 + OH = Cl + H2O
!
      rcarr(103,:) = skarr(  9.200D-13 ,1560.0D+00 ,temperature)
!
!....           HCFC141b + OH = 2 Cl + H2O
!
      rcarr(104,:) = skarr(  1.250D-12 ,1600.0D+00 ,temperature)
!
!....           HCFC142b + OH = Cl + H2O
!
      rcarr(105,:) = skarr(  1.300D-12 ,1770.0D+00 ,temperature)
!
!....           CH3Cl + Cl = CO + 2 HCl + HO2
!
      rcarr(106,:) = skarr(  2.030D-11 ,1110.0D+00 ,temperature)
!
!....           CH3Br + OH = Br + H2O + HO2
!
      rcarr(107,:) = skarr(  1.420D-12 ,1150.0D+00 ,temperature)
!
!....           ALD2 + OH = H2O + MCO3
!
      rcarr(108,:) = skarr(  4.630D-12 ,-350.0D+00 ,temperature)
!
!....           ALD2 + NO3 = HNO3 + MCO3
!
      rcarr(109,:) = skarr(  1.400D-12 ,1900.0D+00 ,temperature)
!
!....           MCO3 + NO2 = PAN
!
      rcarr(110,:) = sktroe(  7.300D-29 ,4.10D0 & 
     &                     , 9.500D-12 ,1.60D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           PAN = MCO3 + NO2
!
      rcarr(111,:) = skpanan (temperature,adcol)
!
!....           MCO3 + NO = MO2 + NO2
!
      rcarr(112,:) = skarr(  8.100D-12 ,-270.0D+00 ,temperature)
!
!....           C2H6 + OH = ETO2 + H2O
!
      rcarr(113,:) = skarr(  7.660D-12 ,1020.0D+00 ,temperature)
!
!....           C2H6 + Cl = ETO2 + HCl
!
      rcarr(114,:) = skarr(  7.200D-11 ,70.0D+00 ,temperature)
!
!....           ETO2 + NO = ALD2 + HO2 + NO2
!
      rcarr(115,:) = sketo2_no_b (temperature ,adcol)
!
!....           C3H8 + OH = B3O2
!
      rcarr(116,:) = skc3h8oh_2 (temperature)
!
!....           C3H8 + OH = A3O2
!
      rcarr(117,:) = skc3h8oh_1 (temperature)
!
!....           A3O2 + NO = HO2 + NO2 + RCHO
!
      rcarr(118,:) = ska3o2_no_b (temperature ,adcol)
!
!....           ATO2 + NO =  0.19 CH2O +  0.77 HO2 +  0.19 MCO3 +  0.96 NO2
!
      rcarr(119,:) = skarr(  2.900D-12 ,-300.0D+00 ,temperature)
!
!....           KO2 + NO =  0.93 ALD2 +  0.93 MCO3 +  0.93 NO2
!
      rcarr(120,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           NO + RIO2 =  0.69 CH2O +  0.86 HO2 + NO2 +  0.14 RIO1
!
      rcarr(121,:) = skro2noabs_2 (temperature ,adcol)
!
!....           NO + RIO2 = HNO3
!
      rcarr(122,:) = skro2noadd_2 (temperature ,adcol)
!
!....           NO + RIO1 =  0.75 CH2O + HO2 + NO2
!
      rcarr(123,:) = skro2noabs_2 (temperature ,adcol)
!
!....           NO + RIO1 = HNO3
!
      rcarr(124,:) = skrio1_no_hno3 (temperature ,adcol)
!
!....           IAO2 + NO =  0.35 CH2O +  0.27 CO +  0.08 HNO3 +  0.92 HO2 +  0.92 NO2
!
      rcarr(125,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           B3O2 + NO = ACET + HO2 + NO2
!
      rcarr(126,:) = skb3o2_no_b (temperature ,adcol)
!
!....           ACTA + OH = H2O + MO2
!
      rcarr(127,:) = skarr(  3.150D-14 ,-920.0D+00 ,temperature)
!
!....           OH + RCHO = H2O + RCO3
!
      rcarr(128,:) = skarr(  6.000D-12 ,-410.0D+00 ,temperature)
!
!....           NO + RCO3 = ETO2 + NO2
!
      rcarr(129,:) = skarr(  6.700D-12 ,-340.0D+00 ,temperature)
!
!....           NO3 + RCHO = HNO3 + RCO3
!
      rcarr(130,:) = skarr(  6.500D-15 ,0.0D+00 ,temperature)
!
!....           ACET + OH = ATO2 + H2O
!
      rcarr(131,:) = skacetoh (temperature)
!
!....           A3O2 + MO2 =  0.75 CH2O + HO2 +  0.75 RCHO
!
      rcarr(132,:) = skarr(  5.920D-13 ,0.0D+00 ,temperature)
!
!....           ATO2 + HO2 = MCO3 + MO2
!
      rcarr(133,:) = skarr(  8.600D-13 ,-700.0D+00 ,temperature)
!
!....           HO2 + RIO2 = RIP
!
      rcarr(134,:) = skino2_ho2 (temperature)
!
!....           HO2 + RIO1 = RIP
!
      rcarr(135,:) = skarr(  2.910D-13 ,-1300.0D+00 ,temperature)
!
!....           HO2 + IAO2 = IAP
!
      rcarr(136,:) = skarr(  2.380D-13 ,-1300.0D+00 ,temperature)
!
!....           B3O2 + HO2 = RB3P
!
      rcarr(137,:) = skb3o2_ho2 (temperature)
!
!....           MEK + OH = H2O + KO2
!
      rcarr(138,:) = skohmek (temperature)
!
!....           ETO2 + MO2 =  0.75 ALD2 +  0.75 CH2O + HO2
!
      rcarr(139,:) = skarr(  3.000D-13 ,0.0D+00 ,temperature)
!
!....           MEK + NO3 = HNO3 + KO2
!
      rcarr(140,:) = skarr(  8.000D-16 ,0.0D+00 ,temperature)
!
!....           ATO2 + MO2 =  0.85 CH2O +  0.90 HO2 +  0.10 MCO3 +  0.25 MEK
!
      rcarr(141,:) = skarr(  7.500D-13 ,-500.0D+00 ,temperature)
!
!....           KO2 + MO2 =  0.50 ALD2 +  0.75 CH2O +  0.50 HO2 +  0.50 MCO3 +  0.25 MEK
!
      rcarr(142,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MO2 + RIO2 =  1.10 CH2O +  0.93 HO2 +  0.25 MEK +  0.07 RIO1
!
      rcarr(143,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MO2 + RIO1 =  1.13 CH2O + HO2 +  0.25 MEK
!
      rcarr(144,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           IAO2 + MO2 =  0.95 CH2O +  0.15 CO + HO2 +  0.25 MEK
!
      rcarr(145,:) = skarr(  1.300D-12 ,0.0D+00 ,temperature)
!
!....           B3O2 + MO2 =  0.75 ACET +  0.75 CH2O + HO2
!
      rcarr(146,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           ETO2 + ETO2 = ALD2
!
      rcarr(147,:) = skarr(  6.800D-14 ,0.0D+00 ,temperature)
!
!....           ETO2 + ETO2 = 2 ALD2 + 2 HO2
!
      rcarr(148,:) = skarr(  6.800D-14 ,0.0D+00 ,temperature)
!
!....           A3O2 + HO2 = RA3P
!
      rcarr(149,:) = ska3o2_ho2 (temperature)
!
!....           HO2 + MCO3 = ACTA + O3
!
      rcarr(150,:) = skho2mco3_1 (temperature)
!
!....           ISOP + OH = RIO2
!
      rcarr(151,:) = skarr(  3.000D-11 ,-360.0D+00 ,temperature)
!
!....           A3O2 + MCO3 = HO2 + MO2 + RCHO
!
      rcarr(152,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           A3O2 + MCO3 = ACTA + RCHO
!
      rcarr(153,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           ISOP + O3 =  0.90 CH2O +  0.05 CO +  0.06 HO2 +  0.10 O3 +  0.27 OH
!
      rcarr(154,:) = skarr(  1.100D-14 ,2000.0D+00 ,temperature)
!
!....           MO2 + RCO3 = CH2O + ETO2 + HO2
!
      rcarr(155,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           MO2 + RCO3 = CH2O
!
      rcarr(156,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO
!
      rcarr(157,:) = skarr(  5.180D-12 ,-200.0D+00 ,temperature)
!
!....           OH + RB3P =  0.50 B3O2 +  0.50 OH +  0.50 RCHO
!
      rcarr(158,:) = skarr(  8.780D-12 ,-200.0D+00 ,temperature)
!
!....           OH + RIP =  0.50 IAO2 +  0.10 RIO1 +  0.40 RIO2
!
      rcarr(159,:) = skarr(  4.750D-12 ,-200.0D+00 ,temperature)
!
!....           IAP + OH =  0.50 IAO2 +  0.50 OH +  0.50 RCHO
!
      rcarr(160,:) = skarr(  5.310D-12 ,-200.0D+00 ,temperature)
!
!....           C2H6 + NO3 = ETO2 + HNO3
!
      rcarr(161,:) = skarr(  1.400D-18 ,0.0D+00 ,temperature)
!
!....           MCO3 + MCO3 = 2 MO2
!
      rcarr(162,:) = skarr(  2.900D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MO2 = CH2O + HO2 + MO2
!
      rcarr(163,:) = skarr(  2.000D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MO2 = ACTA + CH2O
!
      rcarr(164,:) = skarr(  2.000D-13 ,-500.0D+00 ,temperature)
!
!....           KO2 + MCO3 = ALD2 + MCO3 + MO2
!
      rcarr(165,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           ATO2 + MCO3 =  0.20 CH2O +  0.80 HO2 +  0.20 MCO3 + MO2
!
      rcarr(166,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.14 RIO1
!
      rcarr(167,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RIO1 =  0.75 CH2O + HO2 + MO2
!
      rcarr(168,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           IAO2 + MCO3 =  0.40 CH2O +  0.29 CO + HO2 + MO2
!
      rcarr(169,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           B3O2 + MCO3 = ACET + HO2 + MO2
!
      rcarr(170,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           ATO2 + MCO3 = ACTA + MEK
!
      rcarr(171,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           KO2 + MCO3 = ACTA + MEK
!
      rcarr(172,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RIO2 = ACTA + MEK
!
      rcarr(173,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RIO1 = ACTA + MEK
!
      rcarr(174,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           IAO2 + MCO3 = ACTA + MEK
!
      rcarr(175,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           B3O2 + MCO3 = ACET + ACTA
!
      rcarr(176,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           ETO2 + MCO3 = ALD2 + HO2 + MO2
!
      rcarr(177,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           ETO2 + MCO3 = ACTA + ALD2
!
      rcarr(178,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RCO3 = ETO2 + MO2
!
      rcarr(179,:) = skarr(  2.500D-12 ,-500.0D+00 ,temperature)
!
!....           N2O5 = 2 HNO3
!
!... function to calc SO4 percent weight and gammas for ClONO2 + LBs het reacs
      gammas_code = sk_clono2_gammas (temperature, adcol, pressure &
                    , specarr(iclono2,:), specarr(ihcl,:), specarr(ih2o,:) &
                    , FRH, reff_lbs, wt_h2so4, g_clono2, g_clono2_h2o, g_clono2_hcl, g_hocl_hcl)
!
      rcarr(180,:) = sklbs_n2o5 (temperature ,pressure ,sad_lbs ,wt_h2so4 ,ptrop)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(181,:) = sklbs_clono2_h2o (temperature  & 
     &           ,pressure ,sad_lbs ,g_clono2_h2o ,mw(iCLONO2) ,ptrop)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(182,:) = sklbs_brono2 (temperature ,pressure ,sad_lbs ,wt_h2so4 ,ptrop)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(183,:) = sklbs_clono2_hcl (temperature  & 
     &           ,pressure ,sad_lbs ,g_clono2_hcl ,mw(iCLONO2) ,specarr( 38,:) ,ptrop)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(184,:) = sklbs_hocl_hcl (temperature  & 
     &           ,pressure ,sad_lbs ,g_hocl_hcl ,mw(iHOCL) ,specarr( 38,:) ,ptrop)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(185,:) = sklbs_hobr_hcl (temperature  & 
     &           ,adcol ,pressure ,sad_lbs ,specarr(  44,:) ,specarr( 38,:) ,water  & 
     &           ,ptrop)
!
!....           N2O5 = 2 HNO3
!
      rcarr(186,:) = sksts_n2o5 (temperature ,pressure ,sad_sts ,ptrop)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(187,:) = sksts_clono2 (temperature  & 
     &           ,adcol ,pressure ,sad_sts ,specarr(  38,:) ,water ,ptrop)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(188,:) = sksts_brono2 (temperature ,pressure ,sad_sts ,ptrop)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(189,:) = sksts_clono2_hcl (temperature  & 
     &           ,adcol ,pressure ,sad_sts ,specarr(    27,:) ,specarr(  38,:) ,water  & 
     &           ,ptrop)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(190,:) = sksts_hocl_hcl (temperature  & 
     &           ,adcol ,pressure ,sad_sts ,specarr(  45,:) ,specarr( 38,:) ,water  & 
     &           ,ptrop)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(191,:) = sksts_hobr_hcl (temperature  & 
     &           ,adcol ,pressure ,sad_sts ,specarr(  44,:) ,specarr( 38,:) ,water  & 
     &           ,ptrop)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(192,:) = sknat_clono2 (temperature ,pressure ,sad_nat ,ptrop)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(193,:) = sknat_brono2 (temperature ,pressure ,sad_nat ,ptrop)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(194,:) = sknat_hcl_clono2 (temperature  & 
     &           ,pressure ,sad_nat ,specarr(  38,:) ,ptrop)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(195,:) = sknat_hcl_hocl (temperature  & 
     &           ,pressure ,sad_nat ,specarr(  38,:) ,ptrop)
!
!....           BrONO2 + HCl = BrCl + HNO3
!
      rcarr(196,:) = sknat_hcl_brono2 (temperature  & 
     &           ,pressure ,sad_nat ,specarr(  38,:) ,ptrop)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(197,:) = sknat_hcl_hobr (temperature  & 
     &           ,pressure ,sad_nat ,specarr(  38,:) ,ptrop)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(198,:) = skice_clono2 (temperature ,pressure ,sad_ice ,ptrop)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(199,:) = skice_brono2 (temperature ,pressure ,sad_ice ,ptrop)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(200,:) = skice_hcl_clono2 (temperature  & 
     &           ,pressure ,sad_ice ,specarr(  38,:) ,ptrop)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(201,:) = skice_hcl_hocl (temperature  & 
     &           ,pressure ,sad_ice ,specarr(  38,:) ,ptrop)
!
!....           BrONO2 + HCl = BrCl + HNO3
!
      rcarr(202,:) = skice_hcl_brono2 (temperature  & 
     &           ,pressure ,sad_ice ,specarr(  38,:) ,ptrop)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(203,:) = skice_hcl_hobr (temperature  & 
     &           ,pressure ,sad_ice ,specarr(  38,:) ,ptrop)
!
!....           HNO3 = NO2 + OH
!
      rcarr(204,:) = skpyro_hno3 (temperature ,sad_pyro)
!
!....           NO3 + NO3 = 2 NO2 + O2
!
      rcarr(205,:) = skarr(  8.500D-13 ,2450.0D+00 ,temperature)
!
!....           HO2 =  0.50 H2O2
!
      rcarr(206,:) = sktrs_ho2 (temperature, & 
     &            sadcol2, adcol, radA, NSADaer, NSADdust, cPBLcol, pressure)
!
!....           NO2 =  0.50 HNO2 +  0.50 HNO3
!
      rcarr(207,:) = sktrs_no2 (temperature, & 
     &            sadcol2, adcol, radA, NSADaer,NSADdust,ptrop, pressure)
!
!....           NO3 = HNO3
!
      rcarr(208,:) = sktrs_no3 (temperature, & 
     &            sadcol2, adcol, radA,NSADaer,NSADdust,ptrop, pressure)
!
!....           N2O5 = 2 HNO3
!
      rcarr(209,:) = sktrs_n2o5 (temperature, & 
     &           sadcol2,adcol,radA,FRH,NSADaer,NSADdust, ptrop, pressure)
!
!....          End thermal rate constants
!
      CONTAINS
        FUNCTION skarr (af,ae,tk)
          real*8 &
     &      af ,ae ,tk(:)
          real*8, dimension(size(tk)) :: skarr
          skarr(:) = af * exp(-ae / tk(:))
        END FUNCTION skarr
        FUNCTION sklp (af,npwr,tk,ad)
          real*8 af ,npwr ,tk(:) ,ad(:)
          real*8, PARAMETER :: TSTD=298.0d0
          real*8, dimension(size(tk)) :: sklp
          sklp(:) = ad(:) * af * (TSTD/tk(:))**npwr
        END FUNCTION sklp
        FUNCTION skhp (ai,mpwr,tk)
          real*8 ai ,mpwr ,tk(:)
          real*8, PARAMETER :: TSTD=298.0d0
          real*8, dimension(size(tk)) :: skhp
          skhp(:) = ai * (TSTD/tk(:))**mpwr
        END FUNCTION skhp
        FUNCTION skfo (af,npwr,ai,mpwr,tk,ad)
          real*8 &
     &      af ,npwr ,ai ,mpwr ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: skfo
          skfo(:) = sklp(af,npwr,tk,ad) / skhp(ai,mpwr,tk)
        END FUNCTION skfo
        FUNCTION skterlp (af,npwr,ae,tk,ad)
          real*8 &
     &      af ,npwr ,ae ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: skterlp
          skterlp(:) = sklp(af,npwr,tk,ad) * exp(-ae / tk(:))
        END FUNCTION skterlp
        FUNCTION sktroe (af,npwr,ai,mpwr,ae,tk,ad,fc)
          real*8 &
     &      af ,npwr ,ai ,mpwr ,ae ,tk(:) ,ad(:)
          real*8, OPTIONAL :: fc
          real*8, dimension(size(tk)) :: sktroe
          real*8 fsubc
          real*8 skfo_local(size(tk))
          if (present (fc)) then; fsubc=fc; else; fsubc=0.6d0; end if
          skfo_local(:) = skfo(af,npwr,ai,mpwr,tk,ad)
          sktroe(:) = skterlp(af,npwr,ae,tk,ad) * fsubc** &
     &                (1.0d0/(1.0d0+log10(skfo_local(:))**2)) / &
     &                (1.0d0+skfo(af,npwr,ai,mpwr,tk,ad))
        END FUNCTION sktroe
!
!.... Harvard/GMI
        FUNCTION fyrno3(xcarbn,tk,ad)
          real*8  xcarbn
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: fyrno3
          real*8  aaa(size(tk)) ,rarb(size(tk)) ,xxyn(size(tk)) ,yyyn(size(tk)) ,zzyn(size(tk))
!
          xxyn(:)   = 1.94D-22 * exp(0.97d0 * xcarbn) * ad(:) * (300.0d0 / tk(:))**0.0d0
          yyyn(:)   = 0.826d0 * (300.0d0 / tk(:))**8.1d0
          aaa(:)    = log10(xxyn(:) / yyyn(:))
          zzyn(:)   = 1.0d0 / (1.0d0 + aaa(:) * aaa(:))
          rarb(:)   = (xxyn(:) / (1.0d0 + (xxyn(:) / yyyn(:)))) * 0.411d0**zzyn(:)
          fyrno3(:) = rarb(:) / (1.0d0 + rarb(:))
!
        END FUNCTION fyrno3
!.... Harvard/GMI from GEOSCHEM 14.3.1
!
      FUNCTION skro2_no_b (tk, ad, a0, c0, a1)
!
       real*8                      :: tk(:), ad(:)
       real*8, DIMENSION(size(tk)) :: skro2_no_b
       real*8                      :: a0, c0, a1
       real*8, DIMENSION(size(tk)) :: k0,  k, yyyn, xxyn
       real*8, DIMENSION(size(tk)) :: aaa, rarb, zzyn, fyrno3
!
!
! Reaction rate for the "B" branch of these RO2 + NO reactions:
!    ETO2 + NO = NO2 +     HO2 + ...
!    A3O2 + NO = NO2 +     HO2 + ...
!    R4O2 + NO = NO2 + 0.27HO2 + ...
!    B3O2 + NO = NO2 +     HO2 + ...
! in which the "a1" parameter is greater than 1.0.
!
!
       k0(:)     = a0 * EXP( c0 / tk(:) )
       xxyn   = 1.94d-22 * EXP(  0.97d0 * a1 ) * ad(:)
       yyyn(:)   = 0.826d0 * ( (300.d0/tk(:))**8.1d0 )
       aaa(:)    = LOG10( xxyn / yyyn(:) )
       zzyn(:)   = ( 1.0d0 / ( 1.0d0 + ( aaa(:)  * aaa(:)  ) ) )
       rarb(:)   = ( xxyn   / ( 1.0d0 + ( xxyn / yyyn(:) ) ) ) * ( 0.411d0**zzyn(:) )
       fyrno3(:) = ( rarb(:)   / ( 1.0d0 +   rarb(:)          ) )
       skro2_no_b(:) = k0(:) * ( 1.0d0 - fyrno3(:) )
!
      END FUNCTION skro2_no_b
!
!.... Harvard/GMI from GEOSCHEM 14.3.1
!
      FUNCTION skro2_no_a (tk, ad, a0, c0, a1)
!
       real*8                      :: tk(:), ad(:)
       real*8, DIMENSION(size(tk)) :: skro2_no_a
       real*8                      :: a0, c0, a1
       real*8, DIMENSION(size(tk)) :: k0,  k, yyyn, xxyn
       real*8, DIMENSION(size(tk)) :: aaa, rarb, zzyn, fyrno3
!
!
! Reaction rate for the "B" branch of these RO2 + NO reactions:
!    ETO2 + NO = NO2 +     HO2 + ...
!    A3O2 + NO = NO2 +     HO2 + ...
!    R4O2 + NO = NO2 + 0.27HO2 + ...
!    B3O2 + NO = NO2 +     HO2 + ...
! in which the "a1" parameter is greater than 1.0.
!
!
       k0(:)     = a0 * EXP( c0 / tk(:) )
       xxyn   = 1.94d-22 * EXP(  0.97d0 * a1 ) * ad(:)
       yyyn(:)   = 0.826d0 * ( (300.d0/tk(:))**8.1d0 )
       aaa(:)    = LOG10( xxyn / yyyn(:) )
       zzyn(:)   = ( 1.0d0 / ( 1.0d0 + ( aaa(:)  * aaa(:)  ) ) )
       rarb(:)   = ( xxyn   / ( 1.0d0 + ( xxyn / yyyn(:) ) ) ) * ( 0.411d0**zzyn(:) )
       fyrno3(:) = ( rarb(:)   / ( 1.0d0 +   rarb(:)          ) )
       skro2_no_a(:) = k0(:) * fyrno3(:)
!
      END FUNCTION skro2_no_a
!
!.... skho2dis (temperature ,adcol)
!
!_1_
!
!.... JPL 19-5
!
        FUNCTION skho2dis (tk,ad)
!
!....       HO2 + HO2 = H2O2 + O2
!....
!.... This routine returns a bimolecular rate constant that accounts
!.... for the pressure, but not the H2O, dependence of the reaction.
!.... The H2O dependence is treated separately.
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION (size(tk)) :: skho2dis
          skho2dis(:) = 3.0d-13 * exp(460.0d0 / tk(:)) +  &
     &                  2.1d-33 * ad(:) * exp(920.0d0 / tk(:))
        END FUNCTION skho2dis
!
!.... skho2h2o (temperature ,adcol)
!
!_2_
!
!.... JPL 19-5
!
        FUNCTION skho2h2o (tk,ad)
!
!....       HO2 + HO2 + H2O = H2O2 + O2 + H2O
! NOT USED : multiplication factor is given in B13
!            not right here as there is H2O dependence
!....
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION (size(tk)) :: skho2h2o
          skho2h2o(:) = skho2dis(tk ,ad) * 1.4d-21 * exp(2200.d0 / tk(:))
        END FUNCTION skho2h2o
!
!.... skohco (temperature ,adcol)
!
!_3_
!
!.... JPL 19-5 ; CO + OH is composed of two separate reactions
!....  it's density and temperature dependent.
!         M
! OH + CO -> HOCO , but HOCO + O2 -> HO2 + CO2 quickly ; termolecular
!         M
! OH + CO -> H + CO2, but H + O2 -> HO2 quickly ; chemical activation reaction
!
      FUNCTION skohco (tk,ad)
!
!                 M
!...   1) OH + CO = H + CO2; assume H+O2->HO2 quick
!...   2) OH + CO = HO2 + CO2
!...
!... Pressure in hPa
!
      real*8  tk(:), ad(:)
      real*8, DIMENSION (size(tk)) :: skohco
      real*8, DIMENSION (size(tk)) :: k0
      real*8, DIMENSION (size(tk)) :: kinf
      real*8, DIMENSION (size(tk)) :: kint
      real*8, DIMENSION (size(tk)) :: r
      real*8, DIMENSION (size(tk)) :: kf
!
      real*8 fsubc, k0_298 ,npwr ,kinf_298 ,mpwr ,kint_A, kint_B
!
!... start
      fsubc = 0.6d0
!
!... from JPL 19-5 table 2.2
      k0_298   = 6.9d-33
      npwr     = 2.1d0
      kinf_298 = 1.1d-12
      mpwr     = -1.3d0
      kint_A   = 1.85d-13
      kint_B   = 65.0d0
!
!... JPL 19 formulation
      k0(:)   =   k0_298*(298.0/tk(:))**(npwr)
      kinf(:) = kinf_298*(298.0/tk(:))**(mpwr)
      kint(:) = kint_A * exp(-kint_B/tk(:))
! 
      r(:)    = 1.0/(1.0+(log10(k0(:)*ad(:)/kinf(:)))**2)
      kf(:)   = ((kinf(:)*k0(:)*ad(:))/(kinf(:)+k0(:)*ad(:))) * fsubc**r(:)
!
      skohco(:) = kf(:) + kint(:) * (1.0-kf(:)/kinf(:))
!
      END FUNCTION skohco
!
!.... skohch4 (temperature)
! _1_
!
!.... JPL 19-5
!
        FUNCTION skohch4 (tk)
!
!....      OH + CH4 = MO2 + H2O
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skohch4
!
          skohch4(:) = 2.80D-14 * tk(:)**0.667d0 * exp(-1575.0d0 / tk(:))
!
        END FUNCTION skohch4
!
!.... skmo2dis_1 (temperature)
! _2_
!
!.... JPL 19-5
!
        FUNCTION skmo2dis_1 (tk)
!
!....      MO2 + MO2 = MOH + CH2O + O2
!....
!....     k: A=9.50D-14, E/R=-390
!....     k = ka + kb
!....     ka for: MO2 + MO2 = 2 CH2O + 2 HO2
!.... ==> kb for: MO2 + MO2 = MOH + CH2O + O2
!....     ka/kb = 26.2*exp(-1130/T)
!======================================================================
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skmo2dis_1
!
          skmo2dis_1(:) = 9.50D-14 * exp(390.0d0/tk(:)) &
                         / (1.0d0 + 26.2d0 * exp(-1130.0d0/tk(:)))
!
        END FUNCTION skmo2dis_1
!
!.... skmo2dis_2 (temperature)
! _3_
!
!.... JPL 19-5
!
        FUNCTION skmo2dis_2 (tk)
!
!....     k: A=9.50D-14, E/R=-390
!....     k = ka + kb
!.... ==> ka for: MO2 + MO2 = 2 CH2O + 2 HO2
!....     kb for: MO2 + MO2 = MOH + CH2O + O2
!....     ka/kb = 26.2*exp(-1130/T)
!======================================================================
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skmo2dis_2
!
          skmo2dis_2(:) = 9.50D-14 * exp(390.0d0 / tk(:))  &
     &                   / (1.0d0 + 1.0d0 / (26.2d0 * exp(-1130.0d0 / tk(:))))
!
        END FUNCTION skmo2dis_2
!
!.... skono2_d (temperature ,adcol)
!
!_2_
!
!.... JPL 19-5 ; O + NO2 is composed of two separate reactions
!....  it's density and temperature dependent.
!
!         M
! O + NO2 -> NO + O2 , dissociation
!
      FUNCTION skono2_d (tk,ad)
!
!              M
!...   O + NO2 = NO + O2 , dissociation
!...
!... Pressure in hPa
!
      real*8  tk(:), ad(:)
      real*8, DIMENSION (size(tk)) :: skono2_d
      real*8, DIMENSION (size(tk)) :: k0
      real*8, DIMENSION (size(tk)) :: kinf
      real*8, DIMENSION (size(tk)) :: kint
      real*8, DIMENSION (size(tk)) :: r
      real*8, DIMENSION (size(tk)) :: kf
!
      real*8 fsubc, k0_298 ,npwr ,kinf_298 ,mpwr ,kint_A, kint_B
!
!... start
      fsubc = 0.6d0
!
!... from JPL 19-5 table 2.2
      k0_298   = 3.4d-31
      npwr     = 1.6d0
      kinf_298 = 2.3d-11
      mpwr     = 0.2d0
      kint_A   = 5.3d-12
      kint_B   = -200.0d0
!
!... JPL 19 formulation
      k0(:)   =   k0_298*(298.0/tk(:))**(npwr)
      kinf(:) = kinf_298*(298.0/tk(:))**(mpwr)
      kint(:) = kint_A * exp(-kint_B/tk(:)) 
      r(:)    = 1.0/(1.0+(log10(k0(:)*ad(:)/kinf(:)))**2)
      kf(:)   = ((kinf(:)*k0(:)*ad(:))/(kinf(:)+k0(:)*ad(:))) * fsubc**r(:)
!
      skono2_d(:)   = kint(:) * (1.0-(kf(:)/kinf(:)))
!
      END FUNCTION skono2_d
!
!.... skohhno3_j19 (temperature ,adcol)
!
!_2_
!
!.... JPL 19-5 ; OH + HNO3 is composed of two separate reactions
!....  it's density and temperature dependent.
!           M
! OH + HNO3 -> OH.HONO2 goes to NO3 , association
!           M
! OH + HNO3 -> H2O + NO2            , dissociation
!
      FUNCTION skohhno3_j19 (tk,ad)
!
!...
!... Pressure in hPa
!
      real*8  tk(:), ad(:)
      real*8, DIMENSION (size(tk)) :: skohhno3_j19
      real*8, DIMENSION (size(tk)) :: k0
      real*8, DIMENSION (size(tk)) :: kinf
      real*8, DIMENSION (size(tk)) :: kint
      real*8, DIMENSION (size(tk)) :: r
      real*8, DIMENSION (size(tk)) :: kf
!
      real*8 fsubc, k0_298 ,npwr ,kinf_298 ,mpwr ,kint_A, kint_B
!
!... start
      fsubc = 0.6d0
!
!... from JPL 19-5 table 2.2
      k0_298   = 3.9d-31
      npwr     = 7.2d0
      kinf_298 = 1.5d-13
      mpwr     = 4.8d0
      kint_A   = 3.7d-14
      kint_B   = -240.0d0
!
!... JPL 19 formulation
      k0(:)   =   k0_298*(298.0/tk(:))**(npwr)
      kinf(:) = kinf_298*(298.0/tk(:))**(mpwr)
      kint(:) = kint_A * exp(-kint_B/tk(:))
      r(:)    = 1.0/(1.0+(log10(k0(:)*ad(:)/kinf(:)))**2)
      kf(:) = ((kinf(:)*k0(:)*ad(:))/(kinf(:)+k0(:)*ad(:))) * fsubc**r(:)
!
      skohhno3_j19(:) = kf(:) + kint(:) * (1.0-(kf(:)/kinf(:)))
!
      END FUNCTION skohhno3_j19
!
!.... skbrodis (temperature)
!
!_1_
!
!.... JPL 15-10 (unchanged from JPL 00-003)
!
        FUNCTION skbrodis (tk)
!
!....    Combined rate for:
!....       BrO + BrO = 2 Br + O2
!....       BrO + BrO = Br2 + O2
!....
!.... PSC - 8/23/2002
!....
!.... Combined two product channels into one.
!.... Appears as special function to avoid confusion with
!.... actual elementary reaction in reaction database.
!
          real*8  tk(:)
          real*8, DIMENSION (size(tk)) :: skbrodis
          skbrodis(:) = 1.5D-12 * exp(230.0d0 / tk(:))
        END FUNCTION skbrodis
!
!.... skpanan (temperature,adcol)
! _19_
!
!.... Harvard/GMI JPL 19-5
!
      FUNCTION skpanan (tk,ad)
!
!.... PAN = MCO3 + NO2
!....
      real*8 , INTENT(IN) :: tk(:), ad(:)
      real*8, DIMENSION(size(tk)) :: skpanan
      real*8                      :: a0, c0, a1, b1, a2, b2, fv
      real*8, DIMENSION(size(tk)) :: k0, k1, rlow, rhigh, xyrat, blog, fexp
!
!.old       skpanan(:) =  sktroe(9.700D-29,5.60D0,9.300D-12,1.50D0,0.0D0,tk(:),ad(:))  &
!.old     &  / skarr(9.000D-29, -14000.0D0,tk(:))
!... from Harvard 14.3
      a0 = 9.00d-29
      c0 = 14000.0d0
      a1 = 7.3d-29
      b1 = 4.1d0
      a2 = 9.5d-12
      b2 = 1.6d0
      fv = 0.6d0
!
      k0(:) = a0 * EXP( c0 / tk(:) )               ! backwards rxn rate
!
      rlow(:)  = a1 * ( (300.0/tk(:))**b1 ) * ad(:)
      rhigh(:) = a2 * ( (300.0/tk(:))**b2 )
      xyrat(:) = rlow(:) / rhigh(:)
      blog(:)  = LOG10( xyrat(:) )
      fexp(:)  = 1.0d0 / ( 1.0d0 + ( blog(:) * blog(:) ) )
      k1(:)    = rlow(:) * ( fv**fexp(:) ) / ( 1.0d0 + xyrat(:) )
      skpanan(:) = k1(:) / k0(:)
!
      END FUNCTION skpanan
!
!.... sketo2_no_b (temperature ,adcol)
! _33_
!
!.... Harvard/GMI from GEOSCHEM 14.3.1
!
      FUNCTION sketo2_no_b (tk, ad)
!
!.... ETO2 + NO = ETNO3
!
      real*8                      :: tk(:), ad(:)
      real*8, DIMENSION(size(tk)) :: sketo2_no_b
!      real*8, DIMENSION(size(tk)) :: k0,  k, yyyn, xxyn, aaa, rarb, zzyn, fyrno3
      real*8                      :: a0, c0, a1
!
!
      a0 = 2.60d-12
      c0 = 365.0d0
      a1 = 2.0d0
!
      sketo2_no_b(:) =  skro2_no_b (tk, ad, a0, c0, a1)
!
    END FUNCTION sketo2_no_b
!
!.... skc3h8oh_2 (temperature)
! _22_
!
!.... GMI JPL 19-5 - Secondary propane oxidation
!
        FUNCTION skc3h8oh_2 (tk)
!
!....      C3H8 + OH = B3O2
!....
!======================================================================
!
          real*8  tk(:)
          real*8, dimension(size(tk)) :: skc3h8oh_2
!
          skc3h8oh_2(:) = 1.97D-11 * (tk(:)/298.d0)**1.23 * exp(-675.d0/tk(:))
!
!.old          skc3h8ox_1(:) = 7.60D-12 * exp(-585.0d0 / tk(:)) * (1.0d0 /  &
!.old     &                    (1.0d0 + 5.87d0 * (300.0d0/tk(:))**(0.64d0) *  exp(-816.0d0 / tk(:))))
!
        END FUNCTION skc3h8oh_2
!
!.... skc3h8oh_1 (temperature)
! _23_
!
!.... GMI JPL 19-5 - Primary propane oxidation
!
        FUNCTION skc3h8oh_1 (tk)
!
!....      C3H8 + OH = A3O2
!....
!======================================================================
!
          real*8  tk(:)
          real*8, dimension(size(tk)) :: skc3h8oh_1
!
          skc3h8oh_1(:) = 8.54D-13 * (tk(:)/298.d0)**1.54 * exp(-19.d0/tk(:))
!
!.old          skc3h8ox_2(:) = 7.60D-12 * exp(-585.0d0 / tk(:)) * (1.0d0 /  &
!.old     &                    (1.0d0 + 0.17d0 * (300.0d0/tk(:))**(-0.64d0) * exp(816.0d0 / tk(:))))
!
        END FUNCTION skc3h8oh_1
!
!.... ska3o2_no_b (temperature ,adcol)
! _32_
!
!.... Harvard/GMI from GEOSCHEM 14.3.1
!
      FUNCTION ska3o2_no_b (tk, ad)
!
!.... A3O2 + NO = NO2 + HO2 + RCHO 
!
      real*8                      :: tk(:), ad(:)
      real*8, DIMENSION(size(tk)) :: ska3o2_no_b
!      real*8, DIMENSION(size(tk)) :: k0,  k, yyyn, xxyn, aaa, rarb, zzyn, fyrno3
      real*8                      :: a0, c0, a1
!
!
      a0 = 2.90d-12
      c0 = 350.0d0
      a1 = 3.0d0
!
      ska3o2_no_b(:) =  skro2_no_b (tk, ad, a0, c0, a1)
!
    END FUNCTION ska3o2_no_b
!
!.... skro2noabs_2 (temperature ,adcol)
! _14_
!
        FUNCTION skro2noabs_2 (tk,ad)
!
!.... RIO2 + NO = NO2 + 0.864 HO2 + 0.690 CH2O +
!....                   0.402 MVK + 0.288 MACR +
!....                   0.136 RIO1 + 0.127 IALD
!.... RIO1 + NO = NO2 + IALD + HO2 + 0.750 CH2O
!....
!.... PSC - 8/8/2002
!.... SDS - 7/6/2016 updated rate
!....
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noabs_2
!
          skro2noabs_2(:) = 2.70D-12 * exp(350.0d0 / tk(:)) * (1.0d0 - fyrno3(5.0d0,tk,ad))
!
        END FUNCTION skro2noabs_2
!
!.... skro2noadd_2 (temperature ,adcol)
! _17_
!
!.... Harvard/GMI
!
        FUNCTION skro2noadd_2 (tk,ad)
!
!.... RIO2 + NO = HNO3
!.... RIO1 + NO = HNO3
!....
!.... PSC - 8/8/2002
!.... SDS - 7/6/2016 updated rate
!....
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noadd_2
!
          skro2noadd_2(:) = 2.70D-12 * exp(350.0d0 / tk(:)) * fyrno3(5.0d0,tk,ad)
!
        END FUNCTION skro2noadd_2
!
!.... skrio1_no_hno3 (temperature ,adcol)
! _43_
!
!.... Harvard/GMI from GEOSCHEM 14.3.1
!
      FUNCTION skrio1_no_hno3 (tk,ad)
!
! Used to compute the rate for this reaction:
!    RIO1 + NO = NO2 + MVK  + HO2 + CH2O
!
      real*8, intent(in) :: tk(:), ad(:)
      REAL*8             :: a0, b0, c0, rn, x0, y0
      REAL*8, DIMENSION(size(tk)) :: k0, k1, k2, k3 ,k4, skrio1_no_hno3
!
!  FUNCTION GC_NIT(2.7d-12, 3.50d2, 1.190d0, 6.0d0,  1.1644d0,  7.05d-4) 
!
      a0 = 2.7d-12
      b0 = 350.0d0
      c0 = 1.190d0
      rn = 6.0d0
      x0 = 1.1644d0
      y0 = 7.05d-4
!
      k0(:) = 2.0d-22 * EXP( rn ) * ad(:)
      k1(:) = 4.3d-01 * ( tk(:) / 298.0d0)**(-8)
      k1(:) = k0(:) / k1(:)
      k2(:) = ( k0(:) / (1.0+k1(:)) ) * 4.1d-01**( 1.0d0 / ( 1.0d0+(LOG10(k1(:)))**2) )
      k3(:) = k2(:) / ( k2(:) + c0 )
      k4(:) = a0 * ( x0 - tk(:)*y0 )
      skrio1_no_hno3(:)  = k4(:) * EXP( b0 / tk(:) ) * k3(:)
      skrio1_no_hno3(:)  = MAX( skrio1_no_hno3(:), 0.0d0 )
      END FUNCTION skrio1_no_hno3
!
!.... skb3o2_no_b (temperature ,adcol)
! _31_
!
!.... Harvard/GMI from GEOSCHEM 14.3.1
!
      FUNCTION skb3o2_no_b (tk, ad)
!
!.... B3O2 + NO = NO2 + HO2 + ACET 
!
      real*8                      :: tk(:), ad(:)
      real*8, DIMENSION(size(tk)) :: skb3o2_no_b
      real*8                      :: a0, c0, a1
!
!
!.  RCONST(71) = (GC_RO2NO_B2_aca(2.70d-12,360.0d0,3.0d0))
      a0 = 2.70d-12
      c0 = 360.0d0
      a1 = 3.0d0
!
      skb3o2_no_b(:) =  skro2_no_b (tk, ad, a0, c0, a1)
!
    END FUNCTION skb3o2_no_b
!
!.... skacetoh (temperature)
! _24_
!
!.... Harvard/GMI JPL 19-5
!
        FUNCTION skacetoh (tk)
!
!....      ACET + OH = ATO2 + H2O
!....
!======================================================================
!
          real*8  tk(:)
          real*8, dimension(size(tk)) :: skacetoh
!
       skacetoh(:) = 1.330D-13 + skarr(3.8200D-11 ,2000.0D+00 ,tk(:))
!
        END FUNCTION skacetoh
!
!.... skino2_ho2 (temperature)
! _28_
!
!.... Harvard/GMI
!
        FUNCTION skino2_ho2 (tk)
!
!.... INO2 + HO2 = INPN
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skino2_ho2
!
! INO2 +  HO2 => INPN : 
!A  472 2.91E-13  0.0E+00   1300 1 HR  0.00     0.     0.         
!       5.00E+00  0.0E+00      0 0     0.00     0.     0.         
!
       skino2_ho2(:) = 2.91D-13 * exp(1300.0d0 / tk(:)) / (1.0D0+5.00E+00)
!
        END FUNCTION skino2_ho2
!
!.... skb3o2_ho2 (temperature)
! _26_
!
!.... Harvard/GMI from GEOSCHEM 14.3.1
!
      FUNCTION skb3o2_ho2 (tk)
!
!.... B3O2 + HO2 = RB3P
!
      real*8  tk(:)
      real*8, DIMENSION(size(tk)) :: skb3o2_ho2
!
! B3O2 +  HO2 => RB3P : 
!A  472 2.91E-13  0.0E+00   1300 1 HR  0.00     0.     0.         
!       3.00E+00  0.0E+00      0 0     0.00     0.     0.         
!
!.old       skb3o2_ho2(:) = 2.91D-13 * exp(1300.0d0 / tk(:)) / (1.0D0+3.00E+00)
      skb3o2_ho2(:) = 2.91D-13 * exp(1300.0d0 / tk(:)) * ( 1.0d0 - EXP( -0.245d0 * 3.0d0 ) )
!
      END FUNCTION skb3o2_ho2
!
!.... skohmek (temperature)
! _6_
!
!.... Harvard/GMI
!
        FUNCTION skohmek (tk)
!
!....      OH + MEK = KO2 + H2O
!....
!.... PSC - 8/8/2002
!....
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skohmek
!
!.old          skohmek(:) = 2.92D-12 * (300.0d0 / tk(:))**(-2.0d0) * exp(414.0d0 / tk(:))
!... JPL 19-5
          skohmek(:) = 1.33d-13 + 3.82d-11 * exp(-2000.0d0 / tk(:))
!
        END FUNCTION skohmek
!
!.... ska3o2_ho2 (temperature)
! _27_
!
!.... Harvard/GMI from GEOSCHEM 14.3.1
!
      FUNCTION ska3o2_ho2 (tk)
!
!.... A3O2 + HO2 = RA3P
!
      real*8  tk(:)
      real*8  a0, c0, a1
      real*8, DIMENSION(size(tk)) :: ska3o2_ho2
!
! A3O2 +  HO2 => RA3P : 
!A  486 2.91E-13  0.0E+00   1300 1 HR  0.00     0.     0.         
!       3.00E+00  0.0E+00      0 0     0.00     0.     0.         
!
!.older       ska3o2_ho2(:) = 2.91D-13 * exp(1300.0d0 / tk(:)) / (1.0D0+3.00E+00)
!
      a0 = 2.91d-13
      c0 = 1300.0d0
      a1 = 3.0d0
!
      ska3o2_ho2(:) = a0 * EXP( c0 / tk(:) )
      ska3o2_ho2(:) = ska3o2_ho2(:) * ( 1.0d0 - EXP( -0.245d0 * a1 ) )
      END FUNCTION ska3o2_ho2
!
!.... skho2mco3_1 (temperature)
! _7_
!
!.... Harvard/GMI JPL 10-6
!
        FUNCTION skho2mco3_1 (tk)
!
!....      HO2 + MCO3 = ACTA + O3
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skho2mco3_1
!
          skho2mco3_1(:) = 4.30D-13 * exp(1040.0d0 / tk(:)) * (1.0d0 /  &
     &                     (1.0d0 + 37.0d0 * exp(-660.0d0 / tk(:))))
!
        END FUNCTION skho2mco3_1
!
!.... sk_clono2_gammas (temperature ,adcol ,pressure ,specarr(ClONO2,:) ,specarr(HCl,:) ,water ,wt_h2so4, g_clono2, g_clono2_hcl, g_clono2_h2o, g_hocl_hcl)
!
!_13_
!
!.... JPL 19-5
!
      FUNCTION sk_clono2_gammas (tk, ad, pr, clono2, hcl, h2o, FRH, reff &
                  , wt_h2so4, g_clono2, g_clono2_h2o, g_clono2_hcl, g_hocl_hcl)
!
!... Following: Shi, Q., et al, JGR, V106, D20, pp24,259-24,274, OCTOBER 27, 2001.
!
  use ieee_arithmetic
!... return value
      real*8  sk_clono2_gammas
!... input variables
      real*8  tk(:), ad(:), pr(:), clono2(:), hcl(:), h2o(:), FRH(:), reff
!... output variables
      real*8, DIMENSION (size(tk)) :: wt_h2so4, g_clono2, g_clono2_hcl, g_clono2_h2o, g_hocl_hcl
!
!... local variables
      real*8, DIMENSION (size(tk)) :: p_clono2, p_hcl, aw, y1, y2, m, tmp, wt, p0_h2o, p_h2o
      real*8, DIMENSION (size(tk)) :: a1, b1, c1, d1, a2, b2, c2, d2
      real*8, DIMENSION (size(tk)) :: Z1, Z2, Z3, rho, M_h2so4, chi, T0, A, nu, alphaH
      real*8, DIMENSION (size(tk)) :: S_clono2, H_clono2, D_clono2, k_hydr, gamma_h2o
      real*8, DIMENSION (size(tk)) :: H_hcl, k_hcl, l_clono2, f_clono2, gamma_clono2_rxn, gamma_hcl
      real*8, DIMENSION (size(tk)) :: gamma_s, F_hcl, gamma_prime_s, gamma_prime_hcl, gamma_b
      real*8, DIMENSION (size(tk)) :: D_hocl, k_hocl_hcl, H_hocl, l_hocl, f_hocl, gamma_hocl_rxn
      integer :: l
      real*8  :: Rgas
!
!
      sk_clono2_gammas = 999.0
!
!... H2SO4 weight% calc from T and rel hum
      p_clono2(:) = (clono2(:)/ad(:)) * pr(:)/1013.25d0   ! hPa -> Atm
      p_hcl(:)    = (hcl(:)   /ad(:)) * pr(:)/1013.25d0   ! hPa -> Atm
      p_h2o(:)    = (h2o(:)   /ad(:)) * pr(:)             ! hPa
!... [HCL] and [ClONO2] must be gt 0.0, use floor of 1e-30
      where(p_clono2.lt.1.0d-20) p_clono2 = 1.0d-20
      where(p_hcl   .lt.1.0d-20) p_hcl    = 1.0d-20
!
!... H2O saturation partial pressure (hPa)
      p0_h2o(:) = exp(18.452406985 - 3505.1578807/tk(:) &
                   - 330918.55082/(tk(:)**2) &
                   + 12725068.262/(tk(:)**3) )
      aw(:) = p_h2o(:)/p0_h2o(:)
!      aw(:) = FRH(:)
!
      do l=1,size(tk)
        if(aw(l) .le. 0.05) then
           a1(l) = 12.37208932
           b1(l) = -0.16125516114
           c1(l) = -30.490657554
           d1(l) = -2.1133114241
           a2(l) = 13.455394705
           b2(l) = -0.1921312255
           c2(l) = -34.285174607
           d2(l) = -1.7620073078
         elseif(aw(l) .le. 0.85) then
!         elseif(aw(l) .gt. 0.05 .and. aw(l) .le. 0.85) then
           a1(l) = 11.820654354
           b1(l) = -0.20786404244
           c1(l) = -4.807306373
           d1(l) = -5.1727540348
           a2(l) = 12.891938068
           b2(l) = -0.23233847708
           c2(l) = -6.4261237757
           d2(l) = -4.9005471319
         else
           a1(l) = -180.06541028
           b1(l) = -0.38601102592
           c1(l) = -93.317846778
           d1(l) = 273.88132245
           a2(l) = -176.95814097
           b2(l) = -0.36257048154
           c2(l) = -90.469744201
           d2(l) = 267.45509988
         endif
       enddo

      y1(:) = a1(:)*aw(:)**b1(:) + c1(:)*aw(:) + d1(:)
      y2(:) = a2(:)*aw(:)**b2(:) + c2(:)*aw(:) + d2(:)
      m(:) = y1(:) + (tk(:)-190.0d0)*(y2(:)-y1(:))/70.0d0
!... keep m >= 0.0
      where(m(:).lt.0.0d0) m(:) = 0.0d0
!... H2SO4 weight percentage
      wt_h2so4(:) = (9800.0d0*m(:)) / (98.0d0*m(:) + 1000.0d0)
      where(wt_h2so4(:).lt.30.0d0) wt_h2so4(:) = 30.0d0
      where(wt_h2so4(:).gt.80.0d0) wt_h2so4(:) = 80.0d0
      wt(:) = wt_h2so4(:)
!
!... Parameters for the H2SO4 Solution
      Z1(:) = 0.12364 - 5.6d-7*tk(:)**2
      Z2(:) = -0.02954 + 1.814d-7*tk(:)**2
      Z3(:) = 2.343d-3 - 1.487d-6*tk(:) - 1.324d-8*tk(:)**2
!... H2SO4 solution density (g/cm^3)
      rho(:) = 1 + Z1(:)*m(:) + Z2(:)*m(:)**1.5 + Z3(:)*m(:)**2

!... H2SO4 Molarity  (mol/L)
      M_h2so4(:) = rho(:)*wt(:)/9.8
!... H2SO4 Mole Fraction
      chi(:) = wt(:)/(wt(:) + (100.0-wt(:))*98.0/18.0)
!... H2SO4 solution viscosity (cP)
      A(:)  = 169.5  + 5.18*wt(:)  - 0.0825*wt(:)**2 + 3.27d-3*wt(:)**3
      T0(:) = 144.11 + 0.166*wt(:) - 0.0150*wt(:)**2 + 2.18d-4*wt(:)**3
      tmp(:) = (448.0/(tk(:) - T0(:)))
!... limit tmp to prevent overflow
      where((tk(:)-T0(:)).lt.1.0)  tmp(:) = 1.0
      where( tmp(:)      .gt.15.0) tmp(:) = 15.0
!
      nu(:) = A(:) * tk(:)**(-1.43) * exp(tmp(:))
!... acid activity in molarity
      alphaH(:) = exp(60.51 - 0.095*wt(:) + 0.0077*wt(:)**2 - 1.61d-5*wt(:)**3 &
        - (1.76 + 2.52d-4*wt(:)**2)*sqrt(tk(:)) + (-805.89 + 253.05*wt(:)**0.076)/sqrt(tk(:)) )
!
!... Uptake Parameters
      S_clono2(:) = 0.306 + 24.0/tk(:)
      H_clono2(:) = 1.6d-6*exp(4710./tk(:))*exp(-S_clono2(:)*M_h2so4(:))
      D_clono2(:) = 5.0d-8*tk(:)/nu(:)
      k_hydr(:) = (1.95d10*exp(-2800.0/tk(:)))*aw(:) &
                  + (1.22d12*exp(-6200.0/tk(:)))*alphaH(:)*aw(:)
      Rgas = 0.082
      gamma_h2o(:) = 4.0*H_clono2(:)*Rgas*tk(:)*sqrt(D_clono2(:)*k_hydr(:)) &
                     / (1474.0*sqrt(tk(:)))

!...
      H_hcl(:) = (0.094 - 0.61*chi(:) + 1.2*chi(:)**2) &
                 * exp(-8.68 + (8515.0 - 10718.0*chi(:)**0.7)/tk(:))
      k_hcl(:) = 7.9d11*alphaH(:)*D_clono2(:)*H_hcl(:)*p_hcl(:)
!
      l_clono2(:) = sqrt(D_clono2(:)/(k_hydr(:)+k_hcl(:)))
      f_clono2(:) = 1.0/(TANH(reff/l_clono2(:))) - l_clono2(:)/reff
      gamma_clono2_rxn(:) = f_clono2(:)*gamma_h2o(:)*sqrt(1.0+k_hcl(:)/k_hydr(:))
!
      gamma_hcl(:) = gamma_clono2_rxn(:)*k_hcl(:)/(k_hcl(:)+k_hydr(:))
!
      gamma_s(:) = 66.12*exp(-1374.0/tk(:))*H_clono2(:)*H_hcl(:)*p_hcl(:)
      F_hcl(:) = 1.0/(1.0+0.612*(gamma_s(:)+gamma_hcl(:))*p_clono2/p_hcl)
      gamma_prime_s(:) = F_hcl(:)*gamma_s(:)
      gamma_prime_hcl(:) = F_hcl(:)*gamma_hcl(:)
      gamma_b(:) = gamma_prime_hcl(:) + gamma_clono2_rxn(:)*k_hydr(:)/(k_hcl(:)+k_hydr(:))
!
      g_clono2(:) = 1.0/(1.0+1.0/(gamma_prime_s(:)+gamma_b(:)))
!
      g_clono2_hcl(:) = g_clono2(:) &
                        * (gamma_prime_s(:)+gamma_prime_hcl(:)) &
                         / (gamma_prime_s(:)+gamma_b(:))
      g_clono2_h2o(:) = g_clono2(:)-g_clono2_hcl(:)
      where(g_clono2_h2o.lt.0.0d0) g_clono2_h2o = 0.0d0
!
!... now do HOCl+HCl uptake
      D_hocl(:) = 6.4d-8*tk(:)/nu(:)
      k_hocl_hcl(:) = 1.25d9*alphaH(:)*D_hocl(:)*H_hcl(:)*p_hcl(:)
      H_hocl(:) = 1.91d-6*exp(5862.4/tk(:))*exp(-(0.0776+59.18/tk(:))*M_h2so4(:))
      l_hocl(:) = sqrt(D_hocl(:)/k_hocl_hcl(:))
      where(l_hocl(:).gt.1.0d1) l_hocl = 1.0d1
      f_hocl(:) = 1.0/TANH(reff/l_hocl(:))-l_hocl(:)/reff
      gamma_hocl_rxn(:) = f_hocl(:)*4.0*H_hocl(:)*Rgas*tk(:)*sqrt(D_hocl(:)*k_hocl_hcl(:)) &
                          / (2009.0*sqrt(tk(:)))
!... HOCl+HCl gamma
      g_hocl_hcl(:) = 1.0/(1.0+1.0/(gamma_hocl_rxn(:)*F_hcl(:)))
!
      sk_clono2_gammas = 0.0
!
      return
!
      end FUNCTION sk_clono2_gammas
!
!.... sklbs_n2o5 (temperature ,pressure ,sad_lbs ,wt_h2so4 ,ptrop)
!
!_1_
!
!.... (1) JPL 15-10
!
      FUNCTION sklbs_n2o5 (tk,pr,sad,wt,ptrop)
!
      real*8, OPTIONAL :: ptrop
      real*8  tk(:) ,pr(:) ,sad(:), wt(:)
!
      real*8, DIMENSION (size(tk)) :: sklbs_n2o5
      real*8  pi
!.sds..      real*8  gamma
!... update
      real*8, DIMENSION (size(tk)) :: gamma, k0, k1, k2, avgvel
!
      pi = acos(-1.0d0)
!
!=======================================================================
!     N2O5 + stratospheric sulfate aerosol = 2 HNO3
!=======================================================================
!
!.sds..!.... First order reaction rate constant
!.sds..!.... PSC 3/30/99
!.sds..      gamma     = 0.10d0
!
!... updated gamma calc (JPL10-6 or older)
!
!... constant weight % of H2SO4
!.STS..      wt = 60.0d0
!.LBS..      wt = 75.0d0
!
!... JPL10-6 (or earlier?)
      k0(:) = -25.5265 - 0.133188*wt(:) + 0.00930846*wt(:)**2 - 9.0194e-5*wt(:)**3
      k1(:) =  9283.76 +  115.345*wt(:) -    5.19258*wt(:)**2 + 0.0483464*wt(:)**3
      k2(:) = -851801. -  22191.2*wt(:) +    766.916*wt(:)**2 -   6.85427*wt(:)**3
      gamma(:) = exp (k0(:) + k1(:)/tk(:) + k2(:)/(tk(:)**2))
!.end update
!
      avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 &
                   / (pi * mw(IN2O5)))**0.5d0
!
      where( sad > 0.0d0 )
        sklbs_n2o5(:) = 0.25d0 * gamma(:) * avgvel(:) * sad(:)
       elsewhere
        sklbs_n2o5(:) = 0.0d0
       endwhere
!
      if ( present(ptrop) ) then
        where( pr > ptrop ) sklbs_n2o5(:) = 0.0d0
      endif
!
      END FUNCTION sklbs_n2o5
!
!.... sklbs_clono2_h2o (temperature ,pressure ,sad_lbs ,g_clono2_h2o ,mw(iCLONO2) ,ptrop)
!
!_2_
!
!.... JPL 97-4
!
        FUNCTION sklbs_clono2_h2o (tk, pr, sad, gamma, spec_mw, ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:), pr(:), sad(:), gamma(:), spec_mw
          real*8, DIMENSION (size(tk)) :: sklbs_clono2_h2o
          real*8  pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!
          avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 &
                     / (pi * spec_mw))**0.5d0
!
          sklbs_clono2_h2o(:) = 0.25d0 * gamma(:) * avgvel(:) * sad(:)
!
          where( sad <= 0.0d0 ) sklbs_clono2_h2o = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sklbs_clono2_h2o = 0.0d0
          endif
!
        END FUNCTION sklbs_clono2_h2o
!
!.... sklbs_brono2 (temperature ,pressure ,sad_lbs ,wt_h2so4 ,ptrop)
!
!_3_
!
!.... (3) JPL 15-10
!
      FUNCTION sklbs_brono2 (tk,pr,sad,wt,ptrop)
!
      real*8, OPTIONAL :: ptrop
      real*8  tk(:) ,pr(:) ,sad(:)
      real*8, DIMENSION (size(tk)) :: sklbs_brono2
      real*8  pi
      real*8, DIMENSION (size(tk)) :: avgvel, gamma, wt
!... JPL15
!      real*8  wt
!
      pi = acos(-1.0d0)
!
!=======================================================================
!     BrONO2 + stratospheric sulfate aerosol = HOBr + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!.... David Hanson, personal communication, May 13, 1997
!
!.orig          gamma = 0.8d0
!
!... JPL15-6 formulation
!   Hanson has fit an empirical expression for measured gammas for BrONO2 + H2O in the form 
!    of:
!      1/gamma = 1/alpha + 1/gamma(rxn) 
!    where 
!      gamma(rxn) = exp(a+b*wt) 
!      alpha = 0.80,
!      a = 29.2,
!      b = -0.40.
!... assumed %wt for GMI
!      wt = 75.0d0
!
      gamma(:) = 1.0d0 / ( 1.0d0/0.80d0 + 1.0d0/(exp(29.2d0-0.40d0*wt(:))) )
!
      avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 &
               / (pi * mw(IBRONO2)))**0.5d0
!
      sklbs_brono2 = 0.25d0 * gamma(:) * avgvel(:) * sad(:)
!
      where( sad < 0.0d0 ) sklbs_brono2   = 0.0d0
!
      if ( present(ptrop) ) then
        where( pr > ptrop ) sklbs_brono2 = 0.0d0
      end if
!
      END FUNCTION sklbs_brono2
!
!.... sklbs_clono2_hcl (temperature ,pressure ,sad_lbs ,g_clono2_hcl ,mw(iCLONO2) ,specarr(HCl,:) ,ptrop)
!
!_4_
!
!.... JPL 19-5
!
        FUNCTION sklbs_clono2_hcl (tk, pr, sad, gamma, spec_mw, hcl, ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:), pr(:), sad(:), gamma(:), spec_mw, hcl(:) 
          real*8, DIMENSION (size(tk)) :: sklbs_clono2_hcl
          real*8  pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!
          avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 &
                     / (pi * spec_mw))**0.5d0
!
          sklbs_clono2_hcl(:) = 0.25d0 * gamma(:) * avgvel(:) * sad(:) / hcl(:)
!
          where( sad <= 0.0d0 ) sklbs_clono2_hcl = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sklbs_clono2_hcl = 0.0d0
          endif
!
        END FUNCTION sklbs_clono2_hcl
!
!.... sklbs_hocl_hcl (temperature ,pressure ,sad_lbs ,g_hocl_hcl ,mw(iHOCL) ,specarr(HCl,:) ,ptrop)
!
!_5_
!
!.... JPL 19-5
!
        FUNCTION sklbs_hocl_hcl (tk, pr, sad, gamma, spec_mw, hcl, ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:), pr(:), sad(:), gamma(:), spec_mw, hcl(:) 
          real*8, DIMENSION (size(tk)) :: sklbs_hocl_hcl
          real*8  pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!
          avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 &
                     / (pi * spec_mw))**0.5d0
!
          sklbs_hocl_hcl(:) = 0.25d0 * gamma(:) * avgvel(:) * sad(:) / hcl(:)
!
          where( sad <= 0.0d0 .and. hcl < 1.0d0 ) sklbs_hocl_hcl = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sklbs_hocl_hcl = 0.0d0
          endif
!
        END FUNCTION sklbs_hocl_hcl
!
!.... sklbs_hobr_hcl (temperature ,adcol ,pressure ,sad_lbs ,specarr(HOBr,:) ,specarr(HCl,:) ,water ,ptrop)
!
!_6_
!
!.... (6) JPL 97-4
!
        FUNCTION sklbs_hobr_hcl (tk,ad,pr,sad,hobr,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,ad(:) ,pr(:) ,sad(:) ,h2o(:) ,hcl(:) ,hobr(:)
          real*8, DIMENSION (size(tk)) :: sklbs_hobr_hcl
          real*8 adrop ,alpha ,d1 ,hsqrtd ,kii ,minconc ,pi,minadivl
          real*8 adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk)) &
     &     ,fterm(size(tk)) &
     &     ,gcalc(size(tk)) &
     &     ,gprob_tot(size(tk)) &
     &     ,hstar(size(tk)) &
     &     ,k(size(tk)) &
     &     ,ph2o(size(tk)) ,phcl(size(tk)) &
     &     ,tk_150(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HOBr + HCl on stratospheric sulfate aerosol = BrCl + H2O
!=======================================================================
!
!.... First order reaction rate constant
!
!.... Hanson and Ravi, JPC, 98, 5728, 1994
!.... DEK, 1/10/97
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc  = 1.0d0
          adrop    = 1.0D-05
          alpha    = 1.0d0
          d1       = 1.2D-08
          minadivl = 1.00D-15
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:) = (h2o(:) / ad(:)) * pr(:)
!
          ph2o(:) = ph2o(:) / 1013.25d0
!
          where( hcl(:) <= minconc )
            phcl = ((1.0d0 / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          elsewhere
            phcl = ((hcl(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:) = 1013.25d0 * ph2o(:) / 10.0d0**(9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))
!
          where( ah2o(:) > 1.1d0)
            ah2o = 1.1d0
          end where
!
          tk_150(:) = tk(:)
!
          where( tk(:) < 150.0d0 )
            tk_150 = 150.0d0
          endwhere
!
          hstar(:) = exp((6250.0d0 / tk_150(:)) - 10.414d0) * (ah2o(:)**3.49d0)
!
          kii      = 1.0D+05
!
          k(:)     = kii * hstar(:) * phcl(:)
!
          hsqrtd   = 110.0d0
!
          gcalc(:) = 2.2548D-05 * hsqrtd * sqrt(tk(:) * mw(IHOBR) * k(:))
!
          adivl(:) = adrop / sqrt(d1 / k(:))
!
           where( adivl(:) < minadivl)
              adivl = minadivl
           endwhere
!
          fterm(:) = ((exp(adivl(:)) + exp(-adivl(:))) / &
     &                (exp(adivl(:)) - exp(-adivl(:)))) - &
     &               (1.0d0 / adivl(:))
!
!....   NOTE: gprob_tot is the overall uptake coeff for HOCl
!
          where( fterm > 0.0d0 )
            gprob_tot = 1.0d0 / (1.0d0 / (fterm(:) * gcalc(:)) + 1.0d0 / alpha)
          elsewhere
            gprob_tot = 0.0d0
          end where
!
          avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 &
                     / (pi * mw(IHOBR)))**0.5d0
!
          where( hcl > minconc )
            sklbs_hobr_hcl = 0.25d0 * gprob_tot * avgvel * sad / hcl
         elsewhere
!.old            sklbs_hobr_hcl = 0.25d0 * gprob_tot * avgvel * sad
            sklbs_hobr_hcl = 0.0d0
         end where
!
          where( sad < 0.0d0 ) sklbs_hobr_hcl = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sklbs_hobr_hcl = 0.0d0
          endif
!
!... JPL02 has lots of caveats and uncertainity - ignore for now
!!!!!!          sklbs_hobr_hcl = 0.0d0
!
        END FUNCTION sklbs_hobr_hcl
!
!.... sksts_n2o5 (temperature ,pressure ,sad_sts ,ptrop)
!
!_1_
!
!.... (1) JPL 15-10
!
      FUNCTION sksts_n2o5 (tk,pr,sad,ptrop)
!
      real*8, OPTIONAL :: ptrop
      real*8  tk(:) ,pr(:) ,sad(:)
      real*8, DIMENSION (size(tk)) :: sksts_n2o5
      real*8, DIMENSION (size(tk)) :: gamma, avgvel
      real*8  pi
!.sds..      real*8  gamma
!... update
      real*8  wt, k0, k1, k2
!
      pi = acos(-1.0d0)
!
!=======================================================================
!     N2O5 + stratospheric sulfate aerosol = 2 HNO3
!=======================================================================
!
!.sds..!.... First order reaction rate constant
!.sds..!.... PSC 3/30/99
!.sds..      gamma     = 0.10d0
!
!... update gamma calc (older than JPL15)
!
!... weight % of H2SO4
      wt = 60.0d0
!.LBS..      wt = 75.0d0
!
!... JPL10-6 (or earlier?)
      k0 = -25.5265 - 0.133188*wt + 0.00930846*wt**2 - 9.0194e-5*wt**3
      k1 =  9283.76 +  115.345*wt -    5.19258*wt**2 + 0.0483464*wt**3
      k2 = -851801. -  22191.2*wt +    766.916*wt**2 -   6.85427*wt**3
      gamma(:) = exp (k0 + k1/tk(:) + k2/(tk(:)*tk(:)))
!.end update
!
      avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 &
                  / (pi * mw(IN2O5)))**0.5d0
!
      where( sad > 0.0d0 )
        sksts_n2o5 = 0.25d0 * gamma * avgvel * sad
       elsewhere
        sksts_n2o5 = 0.0d0
       endwhere
!
      if ( present(ptrop) ) then
        where( pr > ptrop ) sksts_n2o5 = 0.0d0
       endif
!
      END FUNCTION sksts_n2o5
!
!.... sksts_clono2 (temperature ,adcol ,pressure ,sad_sts ,specarr( HCl,:) ,water ,ptrop)
!
!_2_
!
!.... (2) JPL 97-4
!
        FUNCTION sksts_clono2 (tk,ad,pr,sad,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,ad(:) ,pr(:) ,sad(:) ,h2o(:) ,hcl(:)
          real*8, DIMENSION (size(tk)) :: sksts_clono2
          real*8  adrop ,alpha ,ksur ,minconc ,pi ,ro
          real*8  adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gamma(size(tk)) ,gam0(size(tk)) ,gcalc(size(tk))  &
     &     ,gprob_hcl(size(tk)) ,gprob_tot(size(tk)) ,gsurf(size(tk))  &
     &     ,hstar(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk)) ,prate(size(tk))  &
     &     ,tk_150(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     ClONO2 + stratospheric sulfate aerosol = HOCl + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... PSC 3/30/99
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc = 1.0d0
          alpha   = 1.0d0
          ksur    = 576.0d0
          ro      = 2000.0d0
          adrop   = 1.0d-05
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:) = ((h2o(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
!
          where( hcl(:) <= minconc )
            phcl = ((1.0d0 / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          elsewhere
            phcl = ((hcl(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:) = 1013.25d0 * ph2o(:) / 10.0d0**(9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))
!
          where( ah2o(:) > 1.1d0)
            ah2o = 1.1d0
          end where
!
          tk_150(:) = tk(:)
!
          where( tk(:) < 150.0d0 )
            tk_150 = 150.0d0
          endwhere
!
          hstar(:) = exp((6250.0d0 / tk_150(:)) - 10.414d0) * (ah2o(:)**3.49d0)
!
          gsurf(:) = ah2o(:) * ksur * hstar(:) * phcl(:)
!
          prate(:) = ro * hstar(:) * phcl(:) / ah2o(:)
!
          gam0(:)  = 1.18d-04 + (9.1d-03 * ah2o(:)) + (0.5d0 * ah2o(:)**2.0d0)
!
          gcalc(:) = gam0(:) * sqrt(1.0d0 + prate(:))
!
          adivl(:) = adrop / (1.4d-06 * sqrt(1.0d0 / ah2o(:)))
!
          fterm(:) = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                (exp(adivl(:)) - exp(-adivl(:)))) -  &
     &                (1.0d0 / adivl(:))
!
!....   NOTE: gprob_tot is the overall uptake coeff for ClONO2
!
          gprob_tot(:) = 1.0d0 / (1.0d0 / (gsurf(:) + fterm(:) * gcalc(:)) + 1.0d0 / alpha)
!
          gprob_hcl(:) = gprob_tot(:) * (gsurf(:) + fterm(:) * gcalc(:) * prate(:) /  &
     &                   (1.0d0 + prate(:))) / (gsurf(:) + fterm(:) * gcalc(:))
!
          gamma(:)     = gprob_tot(:) - gprob_hcl(:)
!
          avgvel(:)    = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                   (pi * mw(ICLONO2)))**0.5d0
!
          sksts_clono2 = 0.25d0 * gamma * avgvel * sad
!
          where( sad < 0.0d0 ) sksts_clono2   = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sksts_clono2 = 0.0d0
          end if
!
        END FUNCTION sksts_clono2
!
!.... sksts_brono2 (temperature ,pressure ,sad_sts ,ptrop)
!
!_2_
!
!.... (3) JPL 15-10
!
      FUNCTION sksts_brono2 (tk,pr,sad,ptrop)
!
      real*8, OPTIONAL :: ptrop
      real*8  tk(:) ,pr(:) ,sad(:)
      real*8, DIMENSION (size(tk)) :: sksts_brono2
      real*8  pi
      real*8  avgvel(size(tk))
      real*8  gamma
!... JPL15
      real*8  wt
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     BrONO2 + stratospheric sulfate aerosol = HOBr + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!.... David Hanson, personal communication, May 13, 1997
!
!.orig          gamma = 0.8d0
!
!... JPL15-6 formulation
!   Hanson has fit an empirical expression for measured gammas for BrONO2 + H2O in the form 
!    of:
!      1/gamma = 1/alpha + 1/gamma(rxn) 
!    where 
!      gamma(rxn) = exp(a+b*wt) 
!      alpha = 0.80,
!      a = 29.2,
!      b = -0.40.
!... assumed %wt for GMI STS
      wt = 60.0d0
!
      gamma = 1.0d0 / ( 1.0d0/0.80d0 + 1.0d0/(exp(29.2d0-0.40d0*wt)) )
!
      avgvel = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 &
               / (pi * mw(IBRONO2)))**0.5d0
!
      gamma  = 0.8d0
!
      sksts_brono2 = 0.25d0 * gamma * avgvel * sad
!
      where( sad < 0.0d0 ) sksts_brono2 = 0.0d0
!
      if ( present(ptrop) ) then
        where( pr > ptrop ) sksts_brono2 = 0.0d0
      end if
!
      END FUNCTION sksts_brono2
!
!.... sksts_clono2_hcl (temperature ,adcol ,pressure ,sad_sts ,specarr(ClONO2,:) ,specarr( HCl,:) ,water ,ptrop)
!
!_4_
!
!.... (4) JPL 97-4
!
        FUNCTION sksts_clono2_hcl (tk,ad,pr,sad,clono2,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,ad(:) ,pr(:) ,sad(:) ,h2o(:) ,hcl(:) ,clono2(:)
          real*8, DIMENSION (size(tk)) :: sksts_clono2_hcl
          real*8  adrop ,alpha ,ksur ,minconc ,pi ,ro
          real*8  adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gam0(size(tk)) ,gcalc(size(tk))  &
     &     ,gprob_hcl(size(tk)) ,gprob_tot(size(tk)) ,gsurf(size(tk))  &
     &     ,hstar(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk)) ,prate(size(tk))  &
     &     ,tk_150(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     ClONO2 + HCl on stratospheric sulfate aerosol = Cl2 + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!.... Hanson and Ravi, JPC, 98, 5728, 1994
!.... DEK, 1/10/97
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc = 1.0d0
          alpha   = 1.0d0
          ksur    = 576.0d0
          ro      = 2000.0d0
          adrop   = 1.0d-05
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:) = ((h2o(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
!
          where( hcl(:) <= minconc )
            phcl = ((1.0d0 / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          elsewhere
            phcl = ((hcl(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:) = 1013.25d0 * ph2o(:) / 10.0d0**(9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))
!
          where( ah2o(:) > 1.1d0)
            ah2o = 1.1d0
          end where
!
          tk_150(:) = tk(:)
!
          where( tk(:) < 150.0d0 )
            tk_150 = 150.0d0
          endwhere
!
          hstar(:) = exp((6250.0d0 / tk_150(:)) - 10.414d0) * (ah2o(:)**3.49d0)
!
          gsurf(:) = ah2o(:) * ksur * hstar(:) * phcl(:)
!
          prate(:) = ro * hstar(:) * phcl(:) / ah2o(:)
!
          gam0(:)  = 1.18d-04 + (9.1d-03 * ah2o(:)) + (0.5d0 * ah2o(:)**2.0d0)
!
          gcalc(:) = gam0(:) * sqrt(1.0d0 + prate(:))
!
          adivl(:) = adrop / (1.4d-06 * sqrt(1.0d0 / ah2o(:)))
!
          fterm(:) = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                (exp(adivl(:)) - exp(-adivl(:)))) -  &
     &               (1.0d0 / adivl(:))
!
!....   NOTE: gprob_tot is the overall uptake coeff for ClONO2
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          gprob_tot(:) = 1.0d0 / (1.0d0 /  &
     &                   (gsurf(:) + fterm(:) * gcalc(:)) + 1.0d0 / alpha)
!
          gprob_hcl(:) = gprob_tot(:) * (gsurf(:) +  &
     &                    fterm(:) * gcalc(:) * prate(:) /  &
     &                   (1.0d0 + prate(:))) / (gsurf(:) + fterm(:) * gcalc(:))
!
          avgvel(:)    = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                   (pi * mw(ICLONO2)))**0.5d0
!
          where( hcl > minconc )
            sksts_clono2_hcl   = 0.25d0 * gprob_hcl * avgvel * sad / hcl
          elsewhere
!.old            sksts_clono2_hcl   = 0.25d0 * gprob_hcl * avgvel * sad
            sksts_clono2_hcl   = 0.0d0
          end where
!
          where( sad < 0.0d0 ) sksts_clono2_hcl   = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sksts_clono2_hcl = 0.0d0
          end if
!
        END FUNCTION sksts_clono2_hcl
!
!.... sksts_hocl_hcl (temperature ,adcol ,pressure ,sad_sts ,specarr(HOCl,:) ,specarr(HCl,:) ,water ,ptrop)
!
!_5_
!
!.... (5) JPL 97-4
!
        FUNCTION sksts_hocl_hcl (tk,ad,pr,sad,hocl,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,ad(:) ,pr(:) ,sad(:) ,h2o(:) ,hcl(:) ,hocl(:)
          real*8, DIMENSION (size(tk)) :: sksts_hocl_hcl
          real*8  adrop ,alpha ,d1 ,minconc ,pi,minadivl
          real*8  adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,c1(size(tk)) ,c2(size(tk)) ,c3(size(tk)) ,conv(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gcalc(size(tk))  &
     &     ,gprob_tot(size(tk))  &
     &     ,hhuth(size(tk)) ,hm(size(tk))  &
     &     ,hsqrtd(size(tk)) ,hstar(size(tk)) ,hstar_hocl(size(tk))  &
     &     ,k(size(tk)) ,kii(size(tk))  &
     &     ,mterm(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk))  &
     &     ,rho(size(tk))  &
     &     ,tk_150(size(tk))  &
     &     ,wtper(size(tk))  &
     &     ,z(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HOCl + HCl on stratospheric sulfate aerosol = Cl2 + H2O
!=======================================================================
!
!.... First order reaction rate constant
!
!.... Hanson and Ravi, JPC, 98, 5728, 1994
!.... DEK, 1/10/97
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc  = 1.0d0
          adrop    = 1.0D-05
          alpha    = 1.0d0
          d1       = 9.0D-09
          minadivl = 1.00D-15
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:) = (h2o(:) / ad(:)) * pr(:)
!
          z(:) = log(ph2o(:))
!
          ph2o(:) = ph2o(:) / 1013.25d0
!
!
          where( hcl(:) <= minconc )
            phcl = ((1.0d0 / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          elsewhere
            phcl = ((hcl(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:) = 1013.25d0 * ph2o(:) / 10.0d0**(9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))
!
          where( ah2o(:) > 1.1d0)
            ah2o = 1.1d0
          end where
!
          wtper(:) = ((-14.0508d0 + 0.708928d0 * z(:)) * tk(:) + 3578.6d0) /  &
     &               (45.5374d0 + 1.55981d0 * z(:) - 0.197298d0 * tk(:))
!
          where( wtper < 40.0d0 )
            wtper = 40.0d0
          end where
!
          where( wtper > 80.0d0 )
            wtper = 80.0d0
          end where
!
          kii(:) = exp(2.303d0 * (6.08d0 - 1050.0d0 / tk(:) + 0.0747d0 * wtper(:)))
!
          tk_150(:) = tk(:)
!
          where( tk(:) < 150.0d0 )
            tk_150 = 150.0d0
          endwhere
!
          hstar(:) = exp((6250.0d0 / tk_150(:)) - 10.414d0) * (ah2o(:)**3.49d0)
!
          k(:) = kii(:) * hstar(:) * phcl(:)
!
          adivl(:) = adrop / sqrt(d1 / k(:))
!
           where( adivl(:) < minadivl)
              adivl = minadivl
           endwhere
!
          fterm(:) = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                (exp(adivl(:)) - exp(-adivl(:)))) -  &
     &               (1.0d0 / adivl(:))
!
          mterm(:) = 10.196d0 * wtper(:) / (100.0d0 - wtper(:))
!
          c1(:) = 123.64d0 - 5.6D-04 * tk(:)**2.0d0
!
          c2(:) = -29.54d0 + 1.814D-04 * tk(:)**2.0d0
!
          c3(:) = 2.243d0 - 1.487D-03 * tk(:) + 1.324D-05 * tk(:)**2.0d0
!
          rho(:) = 1000.0d0 + c1(:) * mterm(:) +  &
     &                        c2(:) * mterm(:)**1.5d0 +  &
     &                        c3(:) * mterm(:)**2.0d0
!
          conv(:) = (rho(:) / 1000.0d0) / (1.0d0 + mterm(:) * 0.09808d0)
!
          hhuth(:) = exp(6.4946d0 - mterm(:) *  &
     &                   (-0.04107d0 + 54.56d0 / tk(:)) -  &
     &                   5862.0d0 * (1.0d0 / 298.15d0 - 1.0d0 / tk(:)))
!
          hm(:) = hhuth(:) * conv(:)
!
          hstar_hocl(:) = hm(:) * (1.0d0 + 1.052d0 * exp(0.273d0 * (wtper(:) - 65.66d0)))
!
          hsqrtd(:) = hstar_hocl(:) * sqrt(d1)
!
          gcalc(:) = 2.2548D-05 * hsqrtd(:) * sqrt(tk(:) * mw(IHOCL) * k(:))
!
!
!....   NOTE: gprob_tot is the overall uptake coeff for HOCl
!
          where( fterm > 0.0d0 )
            gprob_tot = 1.0d0 / (1.0d0 / (fterm(:) * gcalc(:)) + 1.0d0 / alpha)
          elsewhere
            gprob_tot = 0.0d0
          end where
!
          avgvel(:)   = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 &
                       / (pi * mw(IHOCL)))**0.5d0
!
          where( hcl > minconc )
            sksts_hocl_hcl = 0.25d0 * gprob_tot * avgvel * sad / hcl
          elsewhere
!.old            sksts_hocl_hcl = 0.25d0 * gprob_tot * avgvel * sad
            sksts_hocl_hcl = 0.0d0
          end where
!
          where( sad < 0.0d0 ) sksts_hocl_hcl   = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sksts_hocl_hcl = 0.0d0
          end if
!
        END FUNCTION sksts_hocl_hcl
!
!.... sksts_hobr_hcl (temperature ,adcol ,pressure ,sad_sts ,specarr(HOBr,:) ,specarr(HCl,:) ,water ,ptrop)
!
!_6_
!
!.... (6) JPL 97-4
!
        FUNCTION sksts_hobr_hcl (tk,ad,pr,sad,hobr,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,ad(:) ,pr(:) ,sad(:) ,h2o(:) ,hcl(:) ,hobr(:)
          real*8, DIMENSION (size(tk)) :: sksts_hobr_hcl
          real*8  adrop ,alpha ,d1 ,hsqrtd ,kii ,minconc ,pi,minadivl
          real*8  adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gcalc(size(tk))  &
     &     ,gprob_tot(size(tk))  &
     &     ,hstar(size(tk))  &
     &     ,k(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk))  &
     &     ,tk_150(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HOBr + HCl on stratospheric sulfate aerosol = BrCl + H2O
!=======================================================================
!
!.... First order reaction rate constant
!
!.... Hanson and Ravi, JPC, 98, 5728, 1994
!.... DEK, 1/10/97
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc  = 1.0d0
          adrop    = 1.0D-05
          alpha    = 1.0d0
          d1       = 1.2D-08
          minadivl = 1.00D-15
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:) = (h2o(:) / ad(:)) * pr(:)
!
          ph2o(:) = ph2o(:) / 1013.25d0
!
          where( hcl(:) <= minconc )
            phcl = ((1.0d0 / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          elsewhere
            phcl = ((hcl(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:) = 1013.25d0 * ph2o(:) / 10.0d0**(9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))
!
          where( ah2o(:) > 1.1d0)
            ah2o = 1.1d0
          end where
!
          tk_150(:) = tk(:)
!
          where( tk(:) < 150.0d0 )
            tk_150 = 150.0d0
          endwhere
!
          hstar(:) = exp((6250.0d0 / tk_150(:)) - 10.414d0) * (ah2o(:)**3.49d0)
!
          kii      = 1.0D+05
!
          k(:)     = kii * hstar(:) * phcl(:)
!
          hsqrtd   = 110.0d0
!
          gcalc(:) = 2.2548D-05 * hsqrtd * sqrt(tk(:) * mw(IHOBR) * k(:))
!
          adivl(:) = adrop / sqrt(d1 / k(:))
!
          where( adivl(:) < minadivl)
             adivl = minadivl
          endwhere
!
          fterm(:) = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                (exp(adivl(:)) - exp(-adivl(:)))) -  &
     &               (1.0d0 / adivl(:))
!
!....   NOTE: gprob_tot is the overall uptake coeff for HOCl
!
          where( fterm > 0.0d0 )
            gprob_tot = 1.0d0 / (1.0d0 / (fterm(:) * gcalc(:)) + 1.0d0 / alpha)
          elsewhere
            gprob_tot = 0.0d0
          end where
!
          avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 &
                     / (pi * mw(IHOBR)))**0.5d0
!
          where( hcl > minconc )
            sksts_hobr_hcl = 0.25d0 * gprob_tot * avgvel * sad / hcl
          elsewhere
!.old            sksts_hobr_hcl = 0.25d0 * gprob_tot * avgvel * sad
            sksts_hobr_hcl = 0.0d0
          end where
!
          where( sad < 0.0d0 ) sksts_hobr_hcl = 0.0d0
!
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sksts_hobr_hcl = 0.0d0
          endif
!
        END FUNCTION sksts_hobr_hcl
!
!.... sknat_clono2 (temperature ,pressure ,sad_nat ,ptrop)
!
!_1_
!
!.... (1) JPL 15-10
!
      FUNCTION sknat_clono2 (tk,pr,sad,ptrop)
!
      real*8, OPTIONAL :: ptrop
      real*8  tk(:) ,pr(:) ,sad(:)
      real*8, DIMENSION(size(tk)) :: sknat_clono2
      real*8  gprob ,pi
      real*8  avgvel(size(tk))
!
      pi = acos(-1.0d0)
!
!=======================================================================
!     ClONO2 + PSC Type I NAT particles = HOCl + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!.... JPL 15-10
      gprob     = 0.004d0
!
      avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
                 (pi * mw(ICLONO2)))**0.5d0
!
      sknat_clono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
      if ( present(ptrop) ) then
        where( pr > ptrop ) sknat_clono2 = 0.0d0
      end if
!
      END FUNCTION sknat_clono2
!
!.... sknat_brono2 (temperature ,pressure ,sad_nat ,ptrop)
!
!_2_
!
!.... (2) JPL 15-10
!
      FUNCTION sknat_brono2 (tk,pr,sad,ptrop)
!
      real*8, OPTIONAL :: ptrop
      real*8  tk(:) ,pr(:) ,sad(:)
      real*8, DIMENSION(size(tk)) :: sknat_brono2
      real*8  gprob ,pi
      real*8  avgvel(size(tk))
!
      pi = acos(-1.0d0)
!
!=======================================================================
!     BrONO2 + PSC Type I NAT particles = HOBr + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!.orig          gprob     = 0.004d0
!... JPL15-10
      gprob     = 0.0d0
!
      avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
                  (pi * mw(IBRONO2)))**0.5d0
!
      sknat_brono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
      if ( present(ptrop) ) then
        where( pr > ptrop ) sknat_brono2 = 0.0d0
      end if
!
      END FUNCTION sknat_brono2
!
!.... sknat_hcl_clono2 (temperature ,pressure ,sad_nat ,specarr( HCl,:) ,ptrop)
!
!_3_
!
!.... (3) JPL 15-10
!
      FUNCTION sknat_hcl_clono2 (tk,pr,sad,hcl,ptrop)
!
      real*8, OPTIONAL :: ptrop
      real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
      real*8, DIMENSION(size(tk)) :: sknat_hcl_clono2
      real*8  gprob ,minconc ,pi
      real*8  avgvel(size(tk))
!
      pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + ClONO2 on PSC Type I NAT particles = Cl2 + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!
      minconc   = 1.0d0
!
!.... JPL 15-10
      gprob     = 0.20d0
!
      avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
                 (pi * mw(ICLONO2)))**0.5d0
!
      sknat_hcl_clono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
      where (hcl < minconc)
!.old        sknat_hcl_clono2  = sknat_hcl_clono2 / minconc
        sknat_hcl_clono2  = 0.0d0
      elsewhere
        sknat_hcl_clono2  = sknat_hcl_clono2 / hcl
      end where
      if ( present(ptrop) ) then
        where( pr > ptrop ) sknat_hcl_clono2 = 0.0d0
      end if
!
      END FUNCTION sknat_hcl_clono2
!
!.... sknat_hcl_hocl (temperature ,pressure ,sad_nat ,specarr( HCl,:) ,ptrop)
!
!_4_
!
!.... (4) JPL 15-10
!
      FUNCTION sknat_hcl_hocl (tk,pr,sad,hcl,ptrop)
!
      real*8, OPTIONAL :: ptrop
      real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
      real*8, DIMENSION(size(tk)) :: sknat_hcl_hocl
      real*8  gprob ,minconc ,pi
      real*8  avgvel(size(tk))
!
      pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + HOCl on PSC Type I NAT particles = Cl2 + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!
      minconc   = 1.0d0
!.... JPL 15-10
      gprob     = 0.10d0
!
      avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
                 (pi * mw(IHOCL)))**0.5d0
!
      sknat_hcl_hocl(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
      where (hcl < minconc)
!.old        sknat_hcl_hocl  = sknat_hcl_hocl / minconc
        sknat_hcl_hocl  = 0.0d0
      elsewhere
        sknat_hcl_hocl  = sknat_hcl_hocl / hcl
      end where
      if ( present(ptrop) ) then
        where( pr > ptrop ) sknat_hcl_hocl = 0.0d0
      end if
!
      END FUNCTION sknat_hcl_hocl
!
!.... sknat_hcl_brono2 (temperature ,pressure ,sad_nat ,specarr( HCl,:) ,ptrop)
!
!_5_
!
!.... (5) JPL 15-10
!
      FUNCTION sknat_hcl_brono2 (tk,pr,sad,hcl,ptrop)
!
      real*8, OPTIONAL :: ptrop
      real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
      real*8, DIMENSION(size(tk)) :: sknat_hcl_brono2
      real*8  gprob ,minconc ,pi
      real*8  avgvel(size(tk))
!
      pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + BrONO2 on PSC Type I NAT particles = BrCl + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!
      minconc   = 1.0d0
!
!... orig          gprob     = 0.20d0
!.... JPL 15-10
      gprob     = 0.0d0
!
      avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
                 (pi * mw(IBRONO2)))**0.5d0
!
      sknat_hcl_brono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
      where (hcl < minconc)
!.old        sknat_hcl_brono2  = sknat_hcl_brono2 / minconc
        sknat_hcl_brono2  = 0.0d0
      elsewhere
        sknat_hcl_brono2  = sknat_hcl_brono2 / hcl
      end where
      if ( present(ptrop) ) then
        where( pr > ptrop ) sknat_hcl_brono2 = 0.0d0
      end if
!
      END FUNCTION sknat_hcl_brono2
!
!.... sknat_hcl_hobr (temperature ,pressure ,sad_nat ,specarr( HCl,:) ,ptrop)
!
!_6_
!
!.... (6) JPL 15-10
!
      FUNCTION sknat_hcl_hobr (tk,pr,sad,hcl,ptrop)
!
      real*8, OPTIONAL :: ptrop
      real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
      real*8, DIMENSION(size(tk)) :: sknat_hcl_hobr
      real*8  gprob ,minconc ,pi
      real*8  avgvel(size(tk))
!
      pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + HOBr on PSC Type I NAT particles = BrCl + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!
      minconc   = 1.0d0
!... orig          gprob     = 0.10d0
!.... JPL 15-10
      gprob     = 0.0d0
!
      avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
                  (pi * mw(IHOBR)))**0.5d0
!
      sknat_hcl_hobr(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
      where (hcl < minconc)
!.old        sknat_hcl_hobr  = sknat_hcl_hobr / minconc
        sknat_hcl_hobr  = 0.0d0
      elsewhere
        sknat_hcl_hobr  = sknat_hcl_hobr / hcl
      end where
      if ( present(ptrop) ) then
        where( pr > ptrop ) sknat_hcl_hobr = 0.0d0
      end if
!
      END FUNCTION sknat_hcl_hobr
!
!.... skice_clono2 (temperature ,pressure ,sad_ice ,ptrop)
!
!_1_
!
!.... (1) JPL 15-10
!
        FUNCTION skice_clono2 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:)
          real*8, DIMENSION(size(tk)) :: skice_clono2
          real*8  gprob ,pi
          real*8  avgvel(size(tk))
!
          pi               = acos(-1.0d0)
!
!=======================================================================
!     ClONO2 + PSC Type II ice particles = HOCl + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!... JPL 15-10
          gprob            = 0.30d0
!
          avgvel(:)        = 100.0d0 *  &
                            (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
                            (pi * mw(ICLONO2)))**0.5d0
!
          skice_clono2(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
          where( sad < 0.0d0 ) skice_clono2   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_clono2 = 0.0d0
          end if
!
        END FUNCTION skice_clono2
!
!.... skice_brono2 (temperature ,pressure ,sad_ice ,ptrop)
!
!_2_
!
!.... (2) JPL 15-10
!
        FUNCTION skice_brono2 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:)
          real*8, DIMENSION(size(tk)) :: skice_brono2
          real*8  gprob ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     BrONO2 + PSC Type II ice particles = HOBr + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!... orig          gprob            = 0.30d0
!... JPL 15-10
          gprob            = 0.26d0
!
          avgvel(:)        = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
                            (pi * mw(IBRONO2)))**0.5d0
!
          skice_brono2(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
          where( sad < 0.0d0 ) skice_brono2   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_brono2 = 0.0d0
          end if
!
        END FUNCTION skice_brono2
!
!.... skice_hcl_clono2 (temperature ,pressure ,sad_ice ,specarr( HCl,:) ,ptrop)
!
!_3_
!
!.... (3) JPL 15-10
!
        FUNCTION skice_hcl_clono2 (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: skice_hcl_clono2
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + ClONO2 on PSC Type II ice particles = Cl2 + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!
          minconc              = 1.0d0
!... JPL 15-10
          gprob                = 0.30d0
!
          avgvel(:)            = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) *  &
                                 1000.0d0 / (pi * mw(ICLONO2)))**0.5d0
!
          skice_hcl_clono2(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
!.old            skice_hcl_clono2   = skice_hcl_clono2 / minconc
            skice_hcl_clono2   = 0.0d0
          elsewhere
            skice_hcl_clono2   = skice_hcl_clono2 / hcl
          end where
          where( sad < 0.0d0 ) skice_hcl_clono2   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_hcl_clono2 = 0.0d0
          end if
!
        END FUNCTION skice_hcl_clono2
!
!.... skice_hcl_hocl (temperature ,pressure ,sad_ice ,specarr( HCl,:) ,ptrop)
!
!_4_
!
!.... (4) JPL 15-10
!
        FUNCTION skice_hcl_hocl (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: skice_hcl_hocl
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + HOCl on PSC Type II ice particles = Cl2 + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!
          minconc            = 1.0d0
!... JPL 15-10
          gprob              = 0.20d0
!
          avgvel(:)          = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
                              (pi * mw(IHOCL)))**0.5d0
!
          skice_hcl_hocl(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
!.old            skice_hcl_hocl   = skice_hcl_hocl / minconc
            skice_hcl_hocl   = 0.0d0
          elsewhere
            skice_hcl_hocl   = skice_hcl_hocl / hcl
          end where
          where( sad < 0.0d0 ) skice_hcl_hocl   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_hcl_hocl = 0.0d0
          end if
!
        END FUNCTION skice_hcl_hocl
!
!.... skice_hcl_brono2 (temperature ,pressure ,sad_ice ,specarr( HCl,:) ,ptrop)
!
!_5_
!
!.... (5) JPL 15-10
!
        FUNCTION skice_hcl_brono2 (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: skice_hcl_brono2
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + BrONO2 on PSC Type II ice particles = BrCl + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!
          minconc              = 1.0d0
!
!.orig          gprob                = 0.30d0
!... JPL 15-10
          gprob                = 0.26d0
!
          avgvel(:)            = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) *  &
                                 1000.0d0 /  &
                                (pi * mw(IBRONO2)))**0.5d0
!
          skice_hcl_brono2(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
!.old            skice_hcl_brono2   = skice_hcl_brono2 / minconc
            skice_hcl_brono2   = 0.0d0
          elsewhere
            skice_hcl_brono2   = skice_hcl_brono2 / hcl
          end where
          where( sad < 0.0d0 ) skice_hcl_brono2   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_hcl_brono2 = 0.0d0
          end if
!
        END FUNCTION skice_hcl_brono2
!
!.... skice_hcl_hobr (temperature ,pressure ,sad_ice ,specarr( HCl,:) ,ptrop)
!
!_6_
!
!.... (6) JPL 15-10
!
        FUNCTION skice_hcl_hobr (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: skice_hcl_hobr
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + HOBr on PSC Type II ice particles = BrCl + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!
          minconc            = 1.0d0
!
!.orig          gprob              = 0.20d0
!... JPL 15-10
          gprob              = 0.30d0
!
          avgvel(:)          = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
                              (pi * mw(IHOBR)))**0.5d0
!
          skice_hcl_hobr(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
!.old            skice_hcl_hobr   = skice_hcl_hobr / minconc
            skice_hcl_hobr   = 0.0d0
          elsewhere
            skice_hcl_hobr   = skice_hcl_hobr / hcl
          end where
          where( sad < 0.0d0 ) skice_hcl_hobr   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_hcl_hobr = 0.0d0
          end if
!
        END FUNCTION skice_hcl_hobr
!
!.... skpyro_hno3 (temperature ,sad_pyro)
!
!.... (1)
!
        FUNCTION skpyro_hno3 (tk,sad)
          real*8  tk(:) ,sad(:)
          real*8, DIMENSION(size(tk)) :: skpyro_hno3
          real*8  gprob ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HNO3 + pyro particles = OH + NO2
!=======================================================================
!
!.... First order reaction rate constant
!.... PSC 1/16/2002
!
          gprob           = 0.0d0
          avgvel(:)       = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                      (pi * mw(IHNO3)))**0.5d0
!
          skpyro_hno3(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
        END FUNCTION skpyro_hno3
!
!.... sktrs_ho2 (temperature, sadcol2, adcol, radA, NSADaer, NSADdust, cPBLcol, pressure)
!
!_5_
!
!
!-------------------------------------------------------------------
!
      FUNCTION HO2( RADIUS, TEMP, DENAIR, SQM, HO2DENS, &
                    AEROTYPE, CONTINENTAL_PBL ) RESULT( GAMMA )
!
!=================================================================
! Internal function HO2 computes the GAMMA reaction probability
! for HO2 loss in aerosols based on the recommendation of 
! Thornton, Jaegle, and McNeill, 
! "Assessing Known Pathways For HO2 Loss in Aqueous Atmospheric
!  Aerosols: Regional and Global Impacts on Tropospheric Oxidants"
!  J. Geophys. Res.,  doi:10.1029/2007JD009236, 2008  
!
! gamma(HO2) is a function of aerosol type, radius, temperature
!
! jaegle 01/22/2008
! 
! Arguments as Input:
! ----------------------------------------------------------------
! (1 ) RADIUS   (REAL*8 ) : Aerosol radius [cm]
! (2 ) TEMP     (REAL*8 ) : Temperature [K]
! (3 ) DENAIR   (REAL*8 ) : Air Density [molec/cm3]
! (4 ) HO2DENS  (REAL*8 ) : HO2 Number Density [molec/cm3]
! (5 ) SQM      (REAL*8 ) : Square root of molecular weight [g/mole]
! (6 ) AEROTYPE (INTEGER) : # denoting aerosol type (cf FAST-J)
! (7 ) CONTINENTAL_PBL (INTEGER)  : Flag set to 1 if the
!         box is located in the continenal boundary layer,
!         otherwise it is zero. Also check for ICE/SNOW (to
!         disable this at high latitudes)
!
! NOTES:
!=================================================================
!
!... Arguments
      REAL*8,  INTENT(IN) :: RADIUS, TEMP, DENAIR, HO2DENS, SQM
      INTEGER, INTENT(IN) :: AEROTYPE, CONTINENTAL_PBL
!... Local variables
      REAL*8              :: ALPHA
      REAL*8              :: delG, Keq, w, H_eff
      REAL*8              :: A1, B1, k1, k2, A, B, C
      REAL*8              :: kaq, kmt, o2_ss, fluxrxn, DFKG
      REAL*8              :: TEST
!... Avogadro's number
      REAL*8,  PARAMETER   :: Na = 6.022d23
!... Ideal gas constant [atm cm3/mol/K], Raq
      REAL*8,  PARAMETER   :: Raq=82.d0

! Function return value
      REAL*8              :: GAMMA
!
!=================================================================
! FUNCTION HO2 begins here!
!=================================================================
!
!... Default value
      GAMMA = 0.0d0
!
!... Special handling for various aerosols
      SELECT CASE ( AEROTYPE )
!
!----------------
! Dust 
!----------------
        CASE ( 1, 2, 3, 4, 5, 6, 7 )      
!----------------
! Assume default gamma=0.1 on dust aerosols
!  This is tentative as no lab measurements presently exist
!  for gamma(HO2) on dust aerosols. We assume the rate to
!  be fast on dust aerosols as transition metal ion induced
!  chemistry is likely to occur in a thin aqueous surface layer.
          GAMMA = 0.1d0

!----------------
! For Sulfate(8), Black Carbon (9), Organic Carbon (10),
! Sea-salt accum & coarse (11,12) calculate the 
! reaction probability due to self reaction 
! by using the algebraic expression in Thornton et al.  (2008)
! (equation 7) which is a function of temperature, aerosol radius,
! air density and HO2 concentration. 
!
! Transition metal ions (such as copper and iron) in sea-salt and 
! carbonaceous aerosols are complexed to ligands and/or exist at 
! a concentration too low to catalyze HO2 loss efficiently, so we 
! apply the HO2 self reaction expression directly for these aerosols.
! 
! In the case of sulfate aerosol, the aerosols likely
! contain copper in the continental boundary layer and
! HO2 uptake proceeds rapidly. To account for the metal catalyzed
! uptake, we assume gamma(HO2)=0.07 (in the mid-range of the recommended
! 0.04-0.1 by Thornton et al, based on observed copper concentrations
! in the US boundary layer). Outside the continental boundary layer, we
! use the HO2-only algebraic expression.
!
!----------------
        CASE ( 8, 9, 10, 11, 12, 13)  
!
!... Mean molecular speed [cm/s]
          w = 14550.5d0 * sqrt(TEMP/(SQM*SQM))
!
!... DFKG = Gas phase diffusion coeff [cm2/s]
          DFKG  = 9.45D17/DENAIR * SQRT(TEMP) * SQRT(3.472D-2 + 1.D0/(SQM*SQM))
!... 
!... calculate T-dependent solubility and aq. reaction rate constants
!...  hydronium ion concentration
!...  A1 = 1.+(Keq/hplus) 
!...  with Keq = 2.1d-5 [M], Equilibrium constant for 
!...  HO2aq = H+ + O2- (Jacob, 2000)
!...       hplus=10.d0^(-pH), with pH = 5
!...  B1 = Req * TEMP
!...  with Req = 1.987d-3 [kcal/K/mol], Ideal gas constant
!...  Note that we assume a constant pH of 5.
          A1 = 1.+ (2.1d-5 / (10.d0**(-5) ) )
          B1 = 1.987d-3 * TEMP
!
!... Free energy change for solvation of HO2 (Hanson 1992, Golden 1991)
!...  in [kcal/mol]:
!...  delG = -4.9-(TEMP-298d0)*delS
!...  with delS=-0.023  [kcal/mol/K],  Entropy change for solvation of HO2
          delG  = -4.9d0 - (TEMP-298.d0) * (-0.023)
          H_eff = exp( -delG / B1 ) * A1
!
!... Estimated temp dependent value for HO2 + O2- (k1) and 
!...  HO2+HO2 (see Jacob 1989)
          k1  =   1.58d10 * exp( -3. / B1 )
          k2  =   2.4d9   * exp( -4.7 / B1 )
          kaq = ( k1 * (A1 - 1.d0) + k2) / (A1**2)
!
!... Calculate the mass transfer rate constant and s.s. conc. of 
!...  total HO2 in the aqueous phase:
!...  kmt = (RADIUS/DFKG + 4d0/w/alpha)^(-1)
!...  with alpha = mass accomodation coefficient, assumed 
!...  to be 1 (Thornton et al.)
          kmt = 1.d0/( RADIUS/DFKG + 4d0/w/1. )
!
!... use quadratic formula to obtain [O2-] in particle of radius RADIUS
          A = -2d0 * kaq
          B = -3d0 * kmt / RADIUS / (H_eff * 0.082 * TEMP)
          C =  3d0 * kmt * HO2DENS * 1000d0 / RADIUS / Na
!
!... Error check that B^2-(4d0*A*C) is not negative
          TEST= B**2-(4d0*A*C)
          IF ( TEST < 0d0 ) THEN
            GAMMA = 0d0
          ELSE
!... Calculate the concentration of O2- in the aerosol
            o2_ss = ( -B  -sqrt(B**2-(4d0*A*C)) )/(2d0*A)
!... Calculate the reactive flux
            fluxrxn = kmt*HO2DENS - o2_ss*Na*kmt/H_eff/Raq/TEMP
            IF ( fluxrxn <= 0d0 ) THEN
              GAMMA = 0d0
            ELSE
!... Gamma for HO2 at TEMP, ho2, and RADIUS given
              GAMMA = 1./( ( ( HO2DENS/fluxrxn ) - ( RADIUS/DFKG ) ) * w / 4.d0 )
            ENDIF
          ENDIF
!...  For sulfate aerosols, check whether we are in
!...  the continental boundary layer, in which case
!...  copper catalyzed HO2 uptake likely dominates and
!...  speeds up the reaction: we assume gamma=0.07,
!...  which is in the middle of the 0.04-0.1 range recommended
!...  by Thornton et al. (2008)
!... 
          IF ( AEROTYPE ==  8 .and. CONTINENTAL_PBL == 1) THEN
            GAMMA = 0.07
          ENDIF 
          IF ( AEROTYPE == 13 .and. CONTINENTAL_PBL == 1) THEN
            GAMMA = 0.07
          ENDIF 
!
!----------------
! Default
!----------------
        CASE DEFAULT
          WRITE (6,*) 'Not a suitable aerosol surface '
          WRITE (6,*) 'for HO2 uptake'
          WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
          STOP
!
      END SELECT
!    
!... If negative value is calculated, set it to zero
      GAMMA = max(GAMMA,0d0)
!
!... Return to sktrs_ho2
      END FUNCTION HO2
!
!.... Harvard/GMI Tropospheric Chemistry
!
      FUNCTION sktrs_ho2 (tk, sad, ad, radA, NSADaer, NSADdust, continental_pbl, pr)
!
      real*8  tk(:) ,sad(:,:), ad(:), radA(:,:), pr(:)
      real*8, DIMENSION (size(tk)) :: sktrs_ho2, gamma
      real*8  pi, sqm
      real*8  avgvel(size(tk)), dfkg(size(tk))
      integer jj,k,ntotA,NSADaer,NSADdust
      integer, DIMENSION (size(tk)) :: continental_pbl
!
      ntotA = NSADaer+NSADdust
!
!=======================================================================
!     HO2 + tropospheric aerosol = 0.5 H2O2
!=======================================================================
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)
! conPBLFlag = 1 if in continental PBL, 0 if not.
!
! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
      sktrs_ho2(:) = 0.d0
      pi           = acos(-1.0d0)
      sqm          = SQRT(mw(IHO2))

!... calculate gas phase diffusion coefficient (cm2/s)
      dfkg (:) = 9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &               + 1.D0/mw(IHO2) )**0.5d0
!
!... calculate mean molecular speed (cm/s)
      avgvel(:) = 100.0d0 *  &
     &                (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                (pi * mw(IHO2)))**0.5d0
!
!... loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA
!
        gamma(:) = 0.00D+00
!
!... obtain gamma when aerosol or dust is present
        do k=1,size(tk)
          if(sad(jj,k) > 0.0d0) &
     &       gamma(k) = HO2(radA(jj,k),tk(k),ad(k),sqm,specarr(IHO2,k),jj,continental_pbl(k))
        enddo
!
!... Bound gamma to avoid division by zero immediately below
        where(gamma(:) < 1.00D-05) gamma(:) = 1.00D-05
!
        where( sad(jj,:) > 0.0d0 )
          sktrs_ho2(:) = sktrs_ho2(:) +  &
     &         sad(jj,:) * ( 4.0d0 / ( gamma(:) * avgvel(:) )+  &
     &         radA(jj,:) / dfkg(:) )**(-1.0d0)
        endwhere
      enddo
!
      END FUNCTION sktrs_ho2
!
!.... sktrs_no2 (temperature, sadcol2, adcol, radA, NSADaer,NSADdust,ptrop, pressure)
!
!_2_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_no2 (tk, sad, ad, radA, NSADaer,NSADdust, ptrop, pr)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,sad(:,:), ad(:), radA(:,:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_no2
          real*8  gamma ,pi
          real*8  avgvel(size(tk)), dfkg(size(tk))
          integer ntotA,jj,NSADaer,NSADdust
!
           ntotA=NSADaer+NSADdust
!
!=======================================================================
!     NO2 + tropospheric aerosol = 0.5 HNO3 + 0.5 HONO
!=======================================================================
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)
!
! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_no2(:)= 0.d0
!
!.sds JPL15?          gamma       = 1.0d-02
          gamma       = 1.0d-04
          pi          = acos(-1.0d0)
!
! calculate gas phase diffusion coefficient (cm2/s)
          dfkg (:) = 9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &              + 1.D0/mw(INO2) )**0.5d0
!
! calculate mean molecular speed (cm/s)
          avgvel(:) = 100.0d0 *  &
     &               (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &               (pi * mw(INO2)))**0.5d0
!
! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA
!
          where( sad(jj,:) > 0.0d0 )
            sktrs_no2(:) = sktrs_no2(:) +  &
     &            sad(jj,:) * ( 4.0d0 / ( gamma * avgvel(:) )+  &
     &            radA(jj,:) / dfkg(:) )**(-1.0d0)
                end where
!
       enddo
!... reaction only used in troposphere if ptrop is passed in
          if ( present(ptrop) ) then
            where( pr <= ptrop ) sktrs_no2 = 0.0d0
          end if
        END FUNCTION sktrs_no2
!
!.... sktrs_no3 (temperature, sadcol2, adcol, radA,NSADaer,NSADdust,ptrop, pressure)
!
!_3_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_no3 (tk, sad, ad, radA, NSADaer,NSADdust, ptrop, pr)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,sad(:,:), ad(:), radA(:,:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_no3
          real*8  gamma ,pi
          real*8  avgvel(size(tk)), dfkg(size(tk))
          integer NSADaer,NSADdust,ntotA,jj
!
          ntotA=NSADaer+NSADdust
!
!=======================================================================
!     NO3 + tropospheric aerosol = HNO3
!=======================================================================
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)
!
! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_no3(:)= 0.d0
!
          gamma       = 1.0d-03
          pi          = acos(-1.0d0)
!
! calculate gas phase diffusion coefficient (cm2/s)
          dfkg(:) = 9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &              + 1.D0/mw(INO3) )**0.5d0
!
! calculate mean molecular speed (cm/s)
          avgvel(:) = 100.0d0 *  &
     &               (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &               (pi * mw(INO3)))**0.5d0
!
! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA
!
          where( sad(jj,:) > 0.0d0 )
            sktrs_no3(:) =  sktrs_no3(:) +  &
     &           sad(jj,:) * ( 4.0d0 / ( gamma * avgvel(:) )+  &
     &           radA(jj,:) / dfkg(:) )**(-1.0d0)
          endwhere
!
       enddo
!... reaction only used in troposphere if ptrop is passed in
          if ( present(ptrop) ) then
            where( pr <= ptrop ) sktrs_no3 = 0.0d0
          end if
        END FUNCTION sktrs_no3
!
!.... sktrs_n2o5 (temperature,sadcol2,adcol,radA,FRH,NSADaer,NSADdust, ptrop, pressure)
!
!_4_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_n2o5 (tk, sad, ad, radA,FRH, NSADaer,NSADdust,ptrop, pr)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,sad(:,:), ad(:), radA(:,:), FRH(:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_n2o5
          real*8  pi, FRH_P(size(tk)),ttk(size(tk)),fact(size(tk))
          real*8  avgvel(size(tk)),dfkg(size(tk)),gamma(size(tk))
          integer NSADaer,NSADdust,ntotA,jj
!
          ntotA=NSADaer+NSADdust
!=======================================================================
! N2O5 + tropospheric aerosol = 2 HNO3                  Branch 4
!=======================================================================
! FRH = relative humidity fraction (0-1)
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)
!
! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_n2o5(:)= 0.d0
!
          pi           = acos(-1.0d0)
!
! calculate gas phase diffusion coefficient (cm2/s)
          dfkg (:) =  9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &                + 1.D0/mw(IN2O5) )**0.5d0
!
! calculate mean molecular speed (cm/s)
          avgvel(:)    = 100.0d0 *  &
     &                   (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                   (pi * mw(IN2O5)))**0.5d0
!
! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA
!
!***********************************************************
! calculate gamma which is function of aerosol type, T, & RH
! following Evans and Jacob, "Impact of new laboratory studies of N2O5
! hydrolysis on global model budgets of tropospheric nitrogen oxides,
! ozone, and OH"
!***********************************************************
      ! Convert RH to % (max = 100%)
      FRH_P(:)  = FRH(:) * 100.d0
      where( FRH_P(:) > 100.d0) FRH_P(:) = 100.d0
!
        gamma(:) = 0d0
!
! DUST
!    Based on unpublished Crowley work
      if(jj.le.7) gamma(:) = 0.01d0
! SULFATE
      if(jj.eq.8) then
!===========================================================
! RH dependence from Kane et al., Heterogenous uptake of
! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
! J. Phys. Chem. A , 2001, 105, 6465-6470
!===========================================================
            gamma(:) = 2.79d-4 + FRH_P(:)*(  1.30d-4 +  &
     &                        FRH_P(:)*( -3.43d-6 +  &
     &                        FRH_P(:)*(  7.52d-8 ) ) )
!
!===========================================================
! Temperature dependence factor (Cox et al, Cambridge UK)
! is of the form:
!
!          10^( LOG10( G294 ) - 0.04 * ( TTEMP - 294 ) )
! FACT = -------------------------------------------------
!                     10^( LOG10( G294 ) )
!
! Where G294 = 1e-2 and TTEMP is MAX( TEMP, 282 ).
!
! For computational speed, replace LOG10( 1e-2 ) with -2
! and replace 10^( LOG10( G294 ) ) with G294
!===========================================================
            ttk(:) = tk(:)
            where( ttk(:) < 282d0) ttk(:) = 282d0
            fact(:) = 10d0**( -2d0 - 4d-2*(ttk(:) - 294.d0))/1d-2
!
            ! Apply temperature dependence
            gamma(:) = gamma(:) * fact(:)
       endif
! BLACK CARBON
!     from IUPAC
      if(jj.eq.9) gamma(:) = 0.005d0
! ORGANIC CARBON
      if(jj.eq.10) then
!===========================================================
! Based on Thornton, Braban and Abbatt, 2003
! N2O5 hydrolysis on sub-micron organic aerosol: the effect
! of relative humidity, particle phase and particle size
!===========================================================
            where ( FRH_P(:) >= 57d0 ) gamma(:) = 0.03d0
            where ( FRH_P(:) <  57d0 ) gamma(:) = FRH_P(:) * 5.2d-4
!
! Bryan & Jules 01/14/05
! Set gamma to very samll number when gamma=0, which occurs
! when relative humidity is 0 (e.g., UT/LS polar night).
!
            where ( gamma(:) == 0.0d0) gamma(:) = 0.03d-2
!
       endif
! SEA SALT
!     Based on IUPAC recomendation
      if(jj.ge.11) then
            where ( FRH_P(:) >= 62 ) gamma(:) = 0.03d0
            where ( FRH_P(:) <  62 ) gamma(:) = 0.005d0
      endif
!
! end calculation of gamma
!***********************************************************
! calculate loss rate
!***********************************************************
          where( sad(jj,:) > 0.0d0 )
            sktrs_n2o5(:) = sktrs_n2o5(:) +  &
     &            sad(jj,:) * ( 4.0d0 / ( gamma(:) * avgvel(:) )+  &
     &            radA(jj,:) / dfkg(:) )**(-1.0d0)
                end where
!
       enddo
!... reaction only used in troposphere if ptrop is passed in
          if ( present(ptrop) ) then
            where( pr <= ptrop ) sktrs_n2o5 = 0.0d0
          end if
!
        END FUNCTION sktrs_n2o5
      END
