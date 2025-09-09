!=======================================================================
!
! $Id: $
!
! ROUTINE
!   Calc_rate_Setkin - GMI version (setkin_ratecalc.F)
!   13 JUN 02 - PSC
!
! DESCRIPTION
!   Calculates and returns rates of kinetic processes for block
!   of boxes
!
! ARGUMENTS
!   INPUT
!    nblock  : number of blocks
!    numj    : number of photolysis rate constants
!    numk    : number of thermal rate constants
!    numspc  ; number of species
!    qk     : thermal reaction rate constants (cm^3 molecule^-1 s^-1)
!    qj     : photolysis rate constants (s^-1)
!    y      : concentrations of transported species (molecules cm^-3)
!   OUTPUT
!    qqk    : rates of thermal processes (molecules cm^-3 s^-1)
!    qqj    : rates of photolytic processes (molecules cm^-3 s^-1)
!
!  Input mechanism:        GeosCCM_Combo_Minimal2_Mechanism.txt
!  Reaction dictionary:    GMI_reactions_JPL19.db
!  Setkin files generated: Wed Mar  5 20:39:38 2025
!
!=======================================================================
      subroutine Calc_rate_Setkin &
     &  (nblock, numj, numk, numspc, &
     &   qk, qj, y, qqk, qqj)

      implicit none

      integer &
     &  nblock, &
     &  numj, numk, numspc

      real*8 &
     &  qj    (nblock, numj), &
     &  qk    (nblock, numk), &
     &  qqj   (nblock, numj), &
     &  qqk   (nblock, numk), &
     &  y     (nblock, numspc)

      integer &
     &  kloop

      do kloop = 1, nblock
!
!                  start thermal reactions
!
!....         O + O2 = O3
!
      qqk(kloop,1)=qk(kloop,1)*y(kloop,60)*y(kloop,75)
!
!....         O + O3 = 2 O2
!
      qqk(kloop,2)=qk(kloop,2)*y(kloop,60)*y(kloop,62)
!
!....         N2 + O1D = N2 + O
!
      qqk(kloop,3)=qk(kloop,3)*y(kloop,74)*y(kloop,61)
!
!....         O1D + O2 = O + O2
!
      qqk(kloop,4)=qk(kloop,4)*y(kloop,61)*y(kloop,75)
!
!....         O1D + O3 = 2 O2
!
      qqk(kloop,5)=qk(kloop,5)*y(kloop,61)*y(kloop,62)
!
!....         O1D + O3 = 2 O + O2
!
      qqk(kloop,6)=qk(kloop,6)*y(kloop,61)*y(kloop,62)
!
!....         H2O + O1D = 2 OH
!
      qqk(kloop,7)=qk(kloop,7)*y(kloop,33)*y(kloop,61)
!
!....         H2 + O1D = H + OH
!
      qqk(kloop,8)=qk(kloop,8)*y(kloop,31)*y(kloop,61)
!
!....         N2O + O1D = N2 + O2
!
      qqk(kloop,9)=qk(kloop,9)*y(kloop,55)*y(kloop,61)
!
!....         N2O + O1D = 2 NO
!
      qqk(kloop,10)=qk(kloop,10)*y(kloop,55)*y(kloop,61)
!
!....         CH4 + O1D = MO2 + OH
!
      qqk(kloop,11)=qk(kloop,11)*y(kloop,22)*y(kloop,61)
!
!....         CH4 + O1D = CH2O + H + HO2
!
      qqk(kloop,12)=qk(kloop,12)*y(kloop,22)*y(kloop,61)
!
!....         CH4 + O1D = CH2O + H2
!
      qqk(kloop,13)=qk(kloop,13)*y(kloop,22)*y(kloop,61)
!
!....         CFC12 + O1D = 2 Cl
!
      qqk(kloop,14)=qk(kloop,14)*y(kloop,17)*y(kloop,61)
!
!....         CFC113 + O1D = 3 Cl
!
      qqk(kloop,15)=qk(kloop,15)*y(kloop,15)*y(kloop,61)
!
!....         HCFC22 + O1D = Cl
!
      qqk(kloop,16)=qk(kloop,16)*y(kloop,37)*y(kloop,61)
!
!....         HCFC141b + O1D = 2 Cl
!
      qqk(kloop,17)=qk(kloop,17)*y(kloop,35)*y(kloop,61)
!
!....         HCFC142b + O1D = Cl
!
      qqk(kloop,18)=qk(kloop,18)*y(kloop,36)*y(kloop,61)
!
!....         H + O2 = HO2
!
      qqk(kloop,19)=qk(kloop,19)*y(kloop,30)*y(kloop,75)
!
!....         H + O3 = O2 + OH
!
      qqk(kloop,20)=qk(kloop,20)*y(kloop,30)*y(kloop,62)
!
!....         O + OH = H + O2
!
      qqk(kloop,21)=qk(kloop,21)*y(kloop,60)*y(kloop,64)
!
!....         HO2 + O = O2 + OH
!
      qqk(kloop,22)=qk(kloop,22)*y(kloop,43)*y(kloop,60)
!
!....         H + HO2 = 2 OH
!
      qqk(kloop,23)=qk(kloop,23)*y(kloop,30)*y(kloop,43)
!
!....         NO + O3 = NO2 + O2
!
      qqk(kloop,24)=qk(kloop,24)*y(kloop,57)*y(kloop,62)
!
!....         O3 + OH = HO2 + O2
!
      qqk(kloop,25)=qk(kloop,25)*y(kloop,62)*y(kloop,64)
!
!....         HO2 + O3 = 2 O2 + OH
!
      qqk(kloop,26)=qk(kloop,26)*y(kloop,43)*y(kloop,62)
!
!....         NO2 + O3 = NO3 + O2
!
      qqk(kloop,27)=qk(kloop,27)*y(kloop,58)*y(kloop,62)
!
!....         OH + OH = H2O + O
!
      qqk(kloop,28)=qk(kloop,28)*y(kloop,64)*y(kloop,64)
!
!....         OH + OH = H2O2
!
      qqk(kloop,29)=qk(kloop,29)*y(kloop,64)*y(kloop,64)
!
!....         HO2 + OH = H2O + O2
!
      qqk(kloop,30)=qk(kloop,30)*y(kloop,43)*y(kloop,64)
!
!....         H2O2 + OH = H2O + HO2
!
      qqk(kloop,31)=qk(kloop,31)*y(kloop,32)*y(kloop,64)
!
!....         HO2 + NO = NO2 + OH
!
      qqk(kloop,32)=qk(kloop,32)*y(kloop,43)*y(kloop,57)
!
!....         HO2 + HO2 = H2O2 + O2
!
      qqk(kloop,33)=qk(kloop,33)*y(kloop,43)*y(kloop,43)
!
!....         H2O + HO2 + HO2 = H2O + H2O2 + O2
!
      qqk(kloop,34)=qk(kloop,34)*y(kloop,33)*y(kloop,43)*y(kloop,43)
!
!....         H2 + OH = H + H2O
!
      qqk(kloop,35)=qk(kloop,35)*y(kloop,31)*y(kloop,64)
!
!....         CO + OH = H
!
      qqk(kloop,36)=qk(kloop,36)*y(kloop,28)*y(kloop,64)
!
!....         CH4 + OH = H2O + MO2
!
      qqk(kloop,37)=qk(kloop,37)*y(kloop,22)*y(kloop,64)
!
!....         ClO + MO2 = CH2O + Cl + HO2 + O2
!
      qqk(kloop,38)=qk(kloop,38)*y(kloop,26)*y(kloop,52)
!
!....         MO2 + NO = CH2O + HO2 + NO2
!
      qqk(kloop,39)=qk(kloop,39)*y(kloop,52)*y(kloop,57)
!
!....         HO2 + MO2 = MP + O2
!
      qqk(kloop,40)=qk(kloop,40)*y(kloop,43)*y(kloop,52)
!
!....         MO2 + MO2 = CH2O + O2
!
      qqk(kloop,41)=qk(kloop,41)*y(kloop,52)*y(kloop,52)
!
!....         MO2 + MO2 = 2 CH2O + 2 HO2
!
      qqk(kloop,42)=qk(kloop,42)*y(kloop,52)*y(kloop,52)
!
!....         MP + OH = H2O + MO2
!
      qqk(kloop,43)=qk(kloop,43)*y(kloop,53)*y(kloop,64)
!
!....         MP + OH = CH2O + H2O + OH
!
      qqk(kloop,44)=qk(kloop,44)*y(kloop,53)*y(kloop,64)
!
!....         CH2O + OH = CO + H2O + HO2
!
      qqk(kloop,45)=qk(kloop,45)*y(kloop,18)*y(kloop,64)
!
!....         N + O2 = NO + O
!
      qqk(kloop,46)=qk(kloop,46)*y(kloop,54)*y(kloop,75)
!
!....         N + NO = N2 + O
!
      qqk(kloop,47)=qk(kloop,47)*y(kloop,54)*y(kloop,57)
!
!....         NO2 + O = NO + O2
!
      qqk(kloop,48)=qk(kloop,48)*y(kloop,58)*y(kloop,60)
!
!....         NO3 + O = NO2 + O2
!
      qqk(kloop,49)=qk(kloop,49)*y(kloop,59)*y(kloop,60)
!
!....         NO2 + OH = HNO3
!
      qqk(kloop,50)=qk(kloop,50)*y(kloop,58)*y(kloop,64)
!
!....         HNO3 + OH = H2O + NO3
!
      qqk(kloop,51)=qk(kloop,51)*y(kloop,41)*y(kloop,64)
!
!....         NO + OH = HNO2
!
      qqk(kloop,52)=qk(kloop,52)*y(kloop,57)*y(kloop,64)
!
!....         HNO2 + OH = H2O + NO2
!
      qqk(kloop,53)=qk(kloop,53)*y(kloop,40)*y(kloop,64)
!
!....         HO2 + NO2 = HNO4
!
      qqk(kloop,54)=qk(kloop,54)*y(kloop,43)*y(kloop,58)
!
!....         HNO4 = HO2 + NO2
!
      qqk(kloop,55)=qk(kloop,55)*y(kloop,42)
!
!....         HNO4 + OH = H2O + NO2 + O2
!
      qqk(kloop,56)=qk(kloop,56)*y(kloop,42)*y(kloop,64)
!
!....         HO2 + NO3 = NO2 + O2 + OH
!
      qqk(kloop,57)=qk(kloop,57)*y(kloop,43)*y(kloop,59)
!
!....         NO + NO3 = 2 NO2
!
      qqk(kloop,58)=qk(kloop,58)*y(kloop,57)*y(kloop,59)
!
!....         NO3 + OH = HO2 + NO2
!
      qqk(kloop,59)=qk(kloop,59)*y(kloop,59)*y(kloop,64)
!
!....         NO2 + NO3 = N2O5
!
      qqk(kloop,60)=qk(kloop,60)*y(kloop,58)*y(kloop,59)
!
!....         N2O5 = NO2 + NO3
!
      qqk(kloop,61)=qk(kloop,61)*y(kloop,56)
!
!....         HCOOH + OH = H2O + HO2
!
      qqk(kloop,62)=qk(kloop,62)*y(kloop,39)*y(kloop,64)
!
!....         NO2 + NO3 = NO + NO2 + O2
!
      qqk(kloop,63)=qk(kloop,63)*y(kloop,58)*y(kloop,59)
!
!....         CH2O + NO3 = CO + HNO3 + HO2
!
      qqk(kloop,64)=qk(kloop,64)*y(kloop,18)*y(kloop,59)
!
!....         Cl + O3 = ClO + O2
!
      qqk(kloop,65)=qk(kloop,65)*y(kloop,25)*y(kloop,62)
!
!....         Cl + H2 = H + HCl
!
      qqk(kloop,66)=qk(kloop,66)*y(kloop,25)*y(kloop,31)
!
!....         Cl + H2O2 = HCl + HO2
!
      qqk(kloop,67)=qk(kloop,67)*y(kloop,25)*y(kloop,32)
!
!....         Cl + HO2 = HCl + O2
!
      qqk(kloop,68)=qk(kloop,68)*y(kloop,25)*y(kloop,43)
!
!....         Cl + HO2 = ClO + OH
!
      qqk(kloop,69)=qk(kloop,69)*y(kloop,25)*y(kloop,43)
!
!....         ClO + O = Cl + O2
!
      qqk(kloop,70)=qk(kloop,70)*y(kloop,26)*y(kloop,60)
!
!....         ClO + OH = Cl + HO2
!
      qqk(kloop,71)=qk(kloop,71)*y(kloop,26)*y(kloop,64)
!
!....         ClO + OH = HCl + O2
!
      qqk(kloop,72)=qk(kloop,72)*y(kloop,26)*y(kloop,64)
!
!....         ClO + HO2 = HOCl + O2
!
      qqk(kloop,73)=qk(kloop,73)*y(kloop,26)*y(kloop,43)
!
!....         ClO + NO = Cl + NO2
!
      qqk(kloop,74)=qk(kloop,74)*y(kloop,26)*y(kloop,57)
!
!....         ClO + NO2 = ClONO2
!
      qqk(kloop,75)=qk(kloop,75)*y(kloop,26)*y(kloop,58)
!
!....         ClO + ClO = 2 Cl + O2
!
      qqk(kloop,76)=qk(kloop,76)*y(kloop,26)*y(kloop,26)
!
!....         ClO + ClO = Cl2 + O2
!
      qqk(kloop,77)=qk(kloop,77)*y(kloop,26)*y(kloop,26)
!
!....         ClO + ClO = Cl + OClO
!
      qqk(kloop,78)=qk(kloop,78)*y(kloop,26)*y(kloop,26)
!
!....         ClO + ClO = Cl2O2
!
      qqk(kloop,79)=qk(kloop,79)*y(kloop,26)*y(kloop,26)
!
!....         Cl2O2 = 2 ClO
!
      qqk(kloop,80)=qk(kloop,80)*y(kloop,24)
!
!....         HCl + OH = Cl + H2O
!
      qqk(kloop,81)=qk(kloop,81)*y(kloop,38)*y(kloop,64)
!
!....         HOCl + OH = ClO + H2O
!
      qqk(kloop,82)=qk(kloop,82)*y(kloop,45)*y(kloop,64)
!
!....         ClONO2 + O = ClO + NO3
!
      qqk(kloop,83)=qk(kloop,83)*y(kloop,27)*y(kloop,60)
!
!....         ClONO2 + OH = HOCl + NO3
!
      qqk(kloop,84)=qk(kloop,84)*y(kloop,27)*y(kloop,64)
!
!....         Cl + ClONO2 = Cl2 + NO3
!
      qqk(kloop,85)=qk(kloop,85)*y(kloop,25)*y(kloop,27)
!
!....         Br + O3 = BrO + O2
!
      qqk(kloop,86)=qk(kloop,86)*y(kloop,6)*y(kloop,62)
!
!....         Br + HO2 = HBr + O2
!
      qqk(kloop,87)=qk(kloop,87)*y(kloop,6)*y(kloop,43)
!
!....         Br + CH2O = CO + HBr + HO2
!
      qqk(kloop,88)=qk(kloop,88)*y(kloop,6)*y(kloop,18)
!
!....         BrO + O = Br + O2
!
      qqk(kloop,89)=qk(kloop,89)*y(kloop,8)*y(kloop,60)
!
!....         BrO + HO2 = HOBr + O2
!
      qqk(kloop,90)=qk(kloop,90)*y(kloop,8)*y(kloop,43)
!
!....         BrO + NO = Br + NO2
!
      qqk(kloop,91)=qk(kloop,91)*y(kloop,8)*y(kloop,57)
!
!....         BrO + NO2 = BrONO2
!
      qqk(kloop,92)=qk(kloop,92)*y(kloop,8)*y(kloop,58)
!
!....         BrO + ClO = Br + OClO
!
      qqk(kloop,93)=qk(kloop,93)*y(kloop,8)*y(kloop,26)
!
!....         BrO + ClO = Br + Cl + O2
!
      qqk(kloop,94)=qk(kloop,94)*y(kloop,8)*y(kloop,26)
!
!....         BrO + ClO = BrCl + O2
!
      qqk(kloop,95)=qk(kloop,95)*y(kloop,8)*y(kloop,26)
!
!....         BrO + BrO = 2 Br + O2
!
      qqk(kloop,96)=qk(kloop,96)*y(kloop,8)*y(kloop,8)
!
!....         HBr + OH = Br + H2O
!
      qqk(kloop,97)=qk(kloop,97)*y(kloop,34)*y(kloop,64)
!
!....         CH2O + O = CO + HO2 + OH
!
      qqk(kloop,98)=qk(kloop,98)*y(kloop,18)*y(kloop,60)
!
!....         CH4 + Cl = HCl + MO2
!
      qqk(kloop,99)=qk(kloop,99)*y(kloop,22)*y(kloop,25)
!
!....         CH2O + Cl = CO + HCl + HO2
!
      qqk(kloop,100)=qk(kloop,100)*y(kloop,18)*y(kloop,25)
!
!....         CH3Cl + OH = Cl + H2O + HO2
!
      qqk(kloop,101)=qk(kloop,101)*y(kloop,21)*y(kloop,64)
!
!....         CH3CCl3 + OH = 3 Cl + H2O
!
      qqk(kloop,102)=qk(kloop,102)*y(kloop,20)*y(kloop,64)
!
!....         HCFC22 + OH = Cl + H2O
!
      qqk(kloop,103)=qk(kloop,103)*y(kloop,37)*y(kloop,64)
!
!....         HCFC141b + OH = 2 Cl + H2O
!
      qqk(kloop,104)=qk(kloop,104)*y(kloop,35)*y(kloop,64)
!
!....         HCFC142b + OH = Cl + H2O
!
      qqk(kloop,105)=qk(kloop,105)*y(kloop,36)*y(kloop,64)
!
!....         CH3Cl + Cl = CO + 2 HCl + HO2
!
      qqk(kloop,106)=qk(kloop,106)*y(kloop,21)*y(kloop,25)
!
!....         CH3Br + OH = Br + H2O + HO2
!
      qqk(kloop,107)=qk(kloop,107)*y(kloop,19)*y(kloop,64)
!
!....         ALD2 + OH = H2O + MCO3
!
      qqk(kloop,108)=qk(kloop,108)*y(kloop,3)*y(kloop,64)
!
!....         ALD2 + NO3 = HNO3 + MCO3
!
      qqk(kloop,109)=qk(kloop,109)*y(kloop,3)*y(kloop,59)
!
!....         MCO3 + NO2 = PAN
!
      qqk(kloop,110)=qk(kloop,110)*y(kloop,50)*y(kloop,58)
!
!....         PAN = MCO3 + NO2
!
      qqk(kloop,111)=qk(kloop,111)*y(kloop,65)
!
!....         MCO3 + NO = MO2 + NO2
!
      qqk(kloop,112)=qk(kloop,112)*y(kloop,50)*y(kloop,57)
!
!....         C2H6 + OH = ETO2 + H2O
!
      qqk(kloop,113)=qk(kloop,113)*y(kloop,10)*y(kloop,64)
!
!....         C2H6 + Cl = ETO2 + HCl
!
      qqk(kloop,114)=qk(kloop,114)*y(kloop,10)*y(kloop,25)
!
!....         ETO2 + NO = ALD2 + HO2 + NO2
!
      qqk(kloop,115)=qk(kloop,115)*y(kloop,29)*y(kloop,57)
!
!....         C3H8 + OH = B3O2
!
      qqk(kloop,116)=qk(kloop,116)*y(kloop,11)*y(kloop,64)
!
!....         C3H8 + OH = A3O2
!
      qqk(kloop,117)=qk(kloop,117)*y(kloop,11)*y(kloop,64)
!
!....         A3O2 + NO = HO2 + NO2 + RCHO
!
      qqk(kloop,118)=qk(kloop,118)*y(kloop,1)*y(kloop,57)
!
!....         ATO2 + NO =  0.19 CH2O +  0.77 HO2 +  0.19 MCO3 +  0.96 NO2
!
      qqk(kloop,119)=qk(kloop,119)*y(kloop,4)*y(kloop,57)
!
!....         KO2 + NO =  0.93 ALD2 +  0.93 MCO3 +  0.93 NO2
!
      qqk(kloop,120)=qk(kloop,120)*y(kloop,49)*y(kloop,57)
!
!....         NO + RIO2 =  0.69 CH2O +  0.86 HO2 + NO2 +  0.14 RIO1
!
      qqk(kloop,121)=qk(kloop,121)*y(kloop,57)*y(kloop,71)
!
!....         NO + RIO2 = HNO3
!
      qqk(kloop,122)=qk(kloop,122)*y(kloop,57)*y(kloop,71)
!
!....         NO + RIO1 =  0.75 CH2O + HO2 + NO2
!
      qqk(kloop,123)=qk(kloop,123)*y(kloop,57)*y(kloop,70)
!
!....         NO + RIO1 = HNO3
!
      qqk(kloop,124)=qk(kloop,124)*y(kloop,57)*y(kloop,70)
!
!....         IAO2 + NO =  0.35 CH2O +  0.27 CO +  0.08 HNO3 +  0.92 HO2 +  0.92 NO2
!
      qqk(kloop,125)=qk(kloop,125)*y(kloop,46)*y(kloop,57)
!
!....         B3O2 + NO = ACET + HO2 + NO2
!
      qqk(kloop,126)=qk(kloop,126)*y(kloop,5)*y(kloop,57)
!
!....         ACTA + OH = H2O + MO2
!
      qqk(kloop,127)=qk(kloop,127)*y(kloop,2)*y(kloop,64)
!
!....         OH + RCHO = H2O + RCO3
!
      qqk(kloop,128)=qk(kloop,128)*y(kloop,64)*y(kloop,68)
!
!....         NO + RCO3 = ETO2 + NO2
!
      qqk(kloop,129)=qk(kloop,129)*y(kloop,57)*y(kloop,69)
!
!....         NO3 + RCHO = HNO3 + RCO3
!
      qqk(kloop,130)=qk(kloop,130)*y(kloop,59)*y(kloop,68)
!
!....         ACET + OH = ATO2 + H2O
!
      qqk(kloop,131)=qk(kloop,131)*y(kloop,73)*y(kloop,64)
!
!....         A3O2 + MO2 =  0.75 CH2O + HO2 +  0.75 RCHO
!
      qqk(kloop,132)=qk(kloop,132)*y(kloop,1)*y(kloop,52)
!
!....         ATO2 + HO2 = MCO3 + MO2
!
      qqk(kloop,133)=qk(kloop,133)*y(kloop,4)*y(kloop,43)
!
!....         HO2 + RIO2 = RIP
!
      qqk(kloop,134)=qk(kloop,134)*y(kloop,43)*y(kloop,71)
!
!....         HO2 + RIO1 = RIP
!
      qqk(kloop,135)=qk(kloop,135)*y(kloop,43)*y(kloop,70)
!
!....         HO2 + IAO2 = IAP
!
      qqk(kloop,136)=qk(kloop,136)*y(kloop,43)*y(kloop,46)
!
!....         B3O2 + HO2 = RB3P
!
      qqk(kloop,137)=qk(kloop,137)*y(kloop,5)*y(kloop,43)
!
!....         MEK + OH = H2O + KO2
!
      qqk(kloop,138)=qk(kloop,138)*y(kloop,51)*y(kloop,64)
!
!....         ETO2 + MO2 =  0.75 ALD2 +  0.75 CH2O + HO2
!
      qqk(kloop,139)=qk(kloop,139)*y(kloop,29)*y(kloop,52)
!
!....         MEK + NO3 = HNO3 + KO2
!
      qqk(kloop,140)=qk(kloop,140)*y(kloop,51)*y(kloop,59)
!
!....         ATO2 + MO2 =  0.85 CH2O +  0.90 HO2 +  0.10 MCO3 +  0.25 MEK
!
      qqk(kloop,141)=qk(kloop,141)*y(kloop,4)*y(kloop,52)
!
!....         KO2 + MO2 =  0.50 ALD2 +  0.75 CH2O +  0.50 HO2 +  0.50 MCO3 +  0.25 MEK
!
      qqk(kloop,142)=qk(kloop,142)*y(kloop,49)*y(kloop,52)
!
!....         MO2 + RIO2 =  1.10 CH2O +  0.93 HO2 +  0.25 MEK +  0.07 RIO1
!
      qqk(kloop,143)=qk(kloop,143)*y(kloop,52)*y(kloop,71)
!
!....         MO2 + RIO1 =  1.13 CH2O + HO2 +  0.25 MEK
!
      qqk(kloop,144)=qk(kloop,144)*y(kloop,52)*y(kloop,70)
!
!....         IAO2 + MO2 =  0.95 CH2O +  0.15 CO + HO2 +  0.25 MEK
!
      qqk(kloop,145)=qk(kloop,145)*y(kloop,46)*y(kloop,52)
!
!....         B3O2 + MO2 =  0.75 ACET +  0.75 CH2O + HO2
!
      qqk(kloop,146)=qk(kloop,146)*y(kloop,5)*y(kloop,52)
!
!....         ETO2 + ETO2 = ALD2
!
      qqk(kloop,147)=qk(kloop,147)*y(kloop,29)*y(kloop,29)
!
!....         ETO2 + ETO2 = 2 ALD2 + 2 HO2
!
      qqk(kloop,148)=qk(kloop,148)*y(kloop,29)*y(kloop,29)
!
!....         A3O2 + HO2 = RA3P
!
      qqk(kloop,149)=qk(kloop,149)*y(kloop,1)*y(kloop,43)
!
!....         HO2 + MCO3 = ACTA + O3
!
      qqk(kloop,150)=qk(kloop,150)*y(kloop,43)*y(kloop,50)
!
!....         ISOP + OH = RIO2
!
      qqk(kloop,151)=qk(kloop,151)*y(kloop,48)*y(kloop,64)
!
!....         A3O2 + MCO3 = HO2 + MO2 + RCHO
!
      qqk(kloop,152)=qk(kloop,152)*y(kloop,1)*y(kloop,50)
!
!....         A3O2 + MCO3 = ACTA + RCHO
!
      qqk(kloop,153)=qk(kloop,153)*y(kloop,1)*y(kloop,50)
!
!....         ISOP + O3 =  0.90 CH2O +  0.05 CO +  0.06 HO2 +  0.10 O3 +  0.27 OH
!
      qqk(kloop,154)=qk(kloop,154)*y(kloop,48)*y(kloop,62)
!
!....         MO2 + RCO3 = CH2O + ETO2 + HO2
!
      qqk(kloop,155)=qk(kloop,155)*y(kloop,52)*y(kloop,69)
!
!....         MO2 + RCO3 = CH2O
!
      qqk(kloop,156)=qk(kloop,156)*y(kloop,52)*y(kloop,69)
!
!....         OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,157)=qk(kloop,157)*y(kloop,64)*y(kloop,66)
!
!....         OH + RB3P =  0.50 B3O2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,158)=qk(kloop,158)*y(kloop,64)*y(kloop,67)
!
!....         OH + RIP =  0.50 IAO2 +  0.10 RIO1 +  0.40 RIO2
!
      qqk(kloop,159)=qk(kloop,159)*y(kloop,64)*y(kloop,72)
!
!....         IAP + OH =  0.50 IAO2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,160)=qk(kloop,160)*y(kloop,47)*y(kloop,64)
!
!....         C2H6 + NO3 = ETO2 + HNO3
!
      qqk(kloop,161)=qk(kloop,161)*y(kloop,10)*y(kloop,59)
!
!....         MCO3 + MCO3 = 2 MO2
!
      qqk(kloop,162)=qk(kloop,162)*y(kloop,50)*y(kloop,50)
!
!....         MCO3 + MO2 = CH2O + HO2 + MO2
!
      qqk(kloop,163)=qk(kloop,163)*y(kloop,50)*y(kloop,52)
!
!....         MCO3 + MO2 = ACTA + CH2O
!
      qqk(kloop,164)=qk(kloop,164)*y(kloop,50)*y(kloop,52)
!
!....         KO2 + MCO3 = ALD2 + MCO3 + MO2
!
      qqk(kloop,165)=qk(kloop,165)*y(kloop,49)*y(kloop,50)
!
!....         ATO2 + MCO3 =  0.20 CH2O +  0.80 HO2 +  0.20 MCO3 + MO2
!
      qqk(kloop,166)=qk(kloop,166)*y(kloop,4)*y(kloop,50)
!
!....         MCO3 + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.14 RIO1
!
      qqk(kloop,167)=qk(kloop,167)*y(kloop,50)*y(kloop,71)
!
!....         MCO3 + RIO1 =  0.75 CH2O + HO2 + MO2
!
      qqk(kloop,168)=qk(kloop,168)*y(kloop,50)*y(kloop,70)
!
!....         IAO2 + MCO3 =  0.40 CH2O +  0.29 CO + HO2 + MO2
!
      qqk(kloop,169)=qk(kloop,169)*y(kloop,46)*y(kloop,50)
!
!....         B3O2 + MCO3 = ACET + HO2 + MO2
!
      qqk(kloop,170)=qk(kloop,170)*y(kloop,5)*y(kloop,50)
!
!....         ATO2 + MCO3 = ACTA + MEK
!
      qqk(kloop,171)=qk(kloop,171)*y(kloop,4)*y(kloop,50)
!
!....         KO2 + MCO3 = ACTA + MEK
!
      qqk(kloop,172)=qk(kloop,172)*y(kloop,49)*y(kloop,50)
!
!....         MCO3 + RIO2 = ACTA + MEK
!
      qqk(kloop,173)=qk(kloop,173)*y(kloop,50)*y(kloop,71)
!
!....         MCO3 + RIO1 = ACTA + MEK
!
      qqk(kloop,174)=qk(kloop,174)*y(kloop,50)*y(kloop,70)
!
!....         IAO2 + MCO3 = ACTA + MEK
!
      qqk(kloop,175)=qk(kloop,175)*y(kloop,46)*y(kloop,50)
!
!....         B3O2 + MCO3 = ACET + ACTA
!
      qqk(kloop,176)=qk(kloop,176)*y(kloop,5)*y(kloop,50)
!
!....         ETO2 + MCO3 = ALD2 + HO2 + MO2
!
      qqk(kloop,177)=qk(kloop,177)*y(kloop,29)*y(kloop,50)
!
!....         ETO2 + MCO3 = ACTA + ALD2
!
      qqk(kloop,178)=qk(kloop,178)*y(kloop,29)*y(kloop,50)
!
!....         MCO3 + RCO3 = ETO2 + MO2
!
      qqk(kloop,179)=qk(kloop,179)*y(kloop,50)*y(kloop,69)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,180)=qk(kloop,180)*y(kloop,56)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,181)=qk(kloop,181)*y(kloop,27)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,182)=qk(kloop,182)*y(kloop,9)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,183)=qk(kloop,183)*y(kloop,27)*y(kloop,38)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,184)=qk(kloop,184)*y(kloop,38)*y(kloop,45)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,185)=qk(kloop,185)*y(kloop,38)*y(kloop,44)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,186)=qk(kloop,186)*y(kloop,56)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,187)=qk(kloop,187)*y(kloop,27)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,188)=qk(kloop,188)*y(kloop,9)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,189)=qk(kloop,189)*y(kloop,27)*y(kloop,38)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,190)=qk(kloop,190)*y(kloop,38)*y(kloop,45)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,191)=qk(kloop,191)*y(kloop,38)*y(kloop,44)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,192)=qk(kloop,192)*y(kloop,27)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,193)=qk(kloop,193)*y(kloop,9)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,194)=qk(kloop,194)*y(kloop,27)*y(kloop,38)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,195)=qk(kloop,195)*y(kloop,38)*y(kloop,45)
!
!....         BrONO2 + HCl = BrCl + HNO3
!
      qqk(kloop,196)=qk(kloop,196)*y(kloop,9)*y(kloop,38)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,197)=qk(kloop,197)*y(kloop,38)*y(kloop,44)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,198)=qk(kloop,198)*y(kloop,27)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,199)=qk(kloop,199)*y(kloop,9)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,200)=qk(kloop,200)*y(kloop,27)*y(kloop,38)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,201)=qk(kloop,201)*y(kloop,38)*y(kloop,45)
!
!....         BrONO2 + HCl = BrCl + HNO3
!
      qqk(kloop,202)=qk(kloop,202)*y(kloop,9)*y(kloop,38)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,203)=qk(kloop,203)*y(kloop,38)*y(kloop,44)
!
!....         HNO3 = NO2 + OH
!
      qqk(kloop,204)=qk(kloop,204)*y(kloop,41)
!
!....         NO3 + NO3 = 2 NO2 + O2
!
      qqk(kloop,205)=qk(kloop,205)*y(kloop,59)*y(kloop,59)
!
!....         HO2 =  0.50 H2O2
!
      qqk(kloop,206)=qk(kloop,206)*y(kloop,43)
!
!....         NO2 =  0.50 HNO2 +  0.50 HNO3
!
      qqk(kloop,207)=qk(kloop,207)*y(kloop,58)
!
!....         NO3 = HNO3
!
      qqk(kloop,208)=qk(kloop,208)*y(kloop,59)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,209)=qk(kloop,209)*y(kloop,56)
!
!                  start photolytic reactions
!
!....  O2 + hv = 2 O
!
      qqj(kloop,1)=qj(kloop,1)*y(kloop,75)
!
!....  O3 + hv = O + O2
!
      qqj(kloop,2)=qj(kloop,2)*y(kloop,62)
!
!....  O3 + hv = O1D + O2
!
      qqj(kloop,3)=qj(kloop,3)*y(kloop,62)
!
!....  NO + hv = N + O
!
      qqj(kloop,4)=qj(kloop,4)*y(kloop,57)
!
!....  N2O + hv = N2 + O1D
!
      qqj(kloop,5)=qj(kloop,5)*y(kloop,55)
!
!....  NO2 + hv = NO + O
!
      qqj(kloop,6)=qj(kloop,6)*y(kloop,58)
!
!....  H2O2 + hv = 2 OH
!
      qqj(kloop,7)=qj(kloop,7)*y(kloop,32)
!
!....  MP + hv = CH2O + HO2 + OH
!
      qqj(kloop,8)=qj(kloop,8)*y(kloop,53)
!
!....  CH2O + hv = CO + H + HO2
!
      qqj(kloop,9)=qj(kloop,9)*y(kloop,18)
!
!....  CH2O + hv = CO + H2
!
      qqj(kloop,10)=qj(kloop,10)*y(kloop,18)
!
!....  HNO3 + hv = NO2 + OH
!
      qqj(kloop,11)=qj(kloop,11)*y(kloop,41)
!
!....  HNO2 + hv = NO + OH
!
      qqj(kloop,12)=qj(kloop,12)*y(kloop,40)
!
!....  HNO4 + hv = NO3 + OH
!
      qqj(kloop,13)=qj(kloop,13)*y(kloop,42)
!
!....  HNO4 + hv = HO2 + NO2
!
      qqj(kloop,14)=qj(kloop,14)*y(kloop,42)
!
!....  NO3 + hv = NO2 + O3
!
      qqj(kloop,15)=qj(kloop,15)*y(kloop,59)
!
!....  NO3 + hv = NO + O2
!
      qqj(kloop,16)=qj(kloop,16)*y(kloop,59)
!
!....  N2O5 + hv = NO2 + NO3
!
      qqj(kloop,17)=qj(kloop,17)*y(kloop,56)
!
!....  N2O5 + hv = NO + NO3 + O3
!
      qqj(kloop,18)=qj(kloop,18)*y(kloop,56)
!
!....  Cl2 + hv = 2 Cl
!
      qqj(kloop,19)=qj(kloop,19)*y(kloop,23)
!
!....  OClO + hv = ClO + O
!
      qqj(kloop,20)=qj(kloop,20)*y(kloop,63)
!
!....  Cl2O2 + hv = 2 Cl + O2
!
      qqj(kloop,21)=qj(kloop,21)*y(kloop,24)
!
!....  HOCl + hv = Cl + OH
!
      qqj(kloop,22)=qj(kloop,22)*y(kloop,45)
!
!....  ClONO2 + hv = Cl + NO3
!
      qqj(kloop,23)=qj(kloop,23)*y(kloop,27)
!
!....  ClONO2 + hv = ClO + NO2
!
      qqj(kloop,24)=qj(kloop,24)*y(kloop,27)
!
!....  BrCl + hv = Br + Cl
!
      qqj(kloop,25)=qj(kloop,25)*y(kloop,7)
!
!....  BrO + hv = Br + O
!
      qqj(kloop,26)=qj(kloop,26)*y(kloop,8)
!
!....  HOBr + hv = Br + OH
!
      qqj(kloop,27)=qj(kloop,27)*y(kloop,44)
!
!....  BrONO2 + hv = Br + NO3
!
      qqj(kloop,28)=qj(kloop,28)*y(kloop,9)
!
!....  BrONO2 + hv = BrO + NO2
!
      qqj(kloop,29)=qj(kloop,29)*y(kloop,9)
!
!....  CH3Cl + hv = Cl + MO2
!
      qqj(kloop,30)=qj(kloop,30)*y(kloop,21)
!
!....  CCl4 + hv = 4 Cl
!
      qqj(kloop,31)=qj(kloop,31)*y(kloop,12)
!
!....  CH3CCl3 + hv = 3 Cl
!
      qqj(kloop,32)=qj(kloop,32)*y(kloop,20)
!
!....  CFC11 + hv = 3 Cl
!
      qqj(kloop,33)=qj(kloop,33)*y(kloop,16)
!
!....  CFC12 + hv = 2 Cl
!
      qqj(kloop,34)=qj(kloop,34)*y(kloop,17)
!
!....  CFC113 + hv = 3 Cl
!
      qqj(kloop,35)=qj(kloop,35)*y(kloop,15)
!
!....  HCFC141b + hv = 2 Cl
!
      qqj(kloop,36)=qj(kloop,36)*y(kloop,35)
!
!....  HCFC142b + hv = Cl
!
      qqj(kloop,37)=qj(kloop,37)*y(kloop,36)
!
!....  CH3Br + hv = Br + MO2
!
      qqj(kloop,38)=qj(kloop,38)*y(kloop,19)
!
!....  CF3Br + hv = Br
!
      qqj(kloop,39)=qj(kloop,39)*y(kloop,14)
!
!....  CF2ClBr + hv = Br + Cl
!
      qqj(kloop,40)=qj(kloop,40)*y(kloop,13)
!
!....  ALD2 + hv = CO + HO2 + MO2
!
      qqj(kloop,41)=qj(kloop,41)*y(kloop,3)
!
!....  ALD2 + hv = CH4 + CO
!
      qqj(kloop,42)=qj(kloop,42)*y(kloop,3)
!
!....  PAN + hv = MCO3 + NO2
!
      qqj(kloop,43)=qj(kloop,43)*y(kloop,65)
!
!....  RCHO + hv = CO + ETO2 + HO2
!
      qqj(kloop,44)=qj(kloop,44)*y(kloop,68)
!
!....  ACET + hv = MCO3 + MO2
!
      qqj(kloop,45)=qj(kloop,45)*y(kloop,73)
!
!....  RA3P + hv = HO2 + OH + RCHO
!
      qqj(kloop,46)=qj(kloop,46)*y(kloop,66)
!
!....  RB3P + hv = HO2 + OH + RCHO
!
      qqj(kloop,47)=qj(kloop,47)*y(kloop,67)
!
!....  RIP + hv =  0.69 CH2O +  0.86 HO2 + OH +  0.14 RIO1
!
      qqj(kloop,48)=qj(kloop,48)*y(kloop,72)
!
!....  IAP + hv =  0.67 CO + HO2 + OH
!
      qqj(kloop,49)=qj(kloop,49)*y(kloop,47)
      enddo
!
!.... End of kloop
!
      return
      end
