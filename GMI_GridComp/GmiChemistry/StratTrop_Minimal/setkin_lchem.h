!=======================================================================
!
! $Id: $
!
! FILE
!   setkin_lchem.h - character labels for species and reactions
!             (setkin_lchem.h)
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   Include file that provides ascii strings identifying reactions
!   and species
!
!  Input mechanism:        GeosCCM_Combo_Minimal2_Mechanism.txt
!  Reaction dictionary:    GMI_reactions_JPL19.db
!  Setkin files generated: Wed Mar  5 20:39:38 2025
!
!=======================================================================


      integer kmg_i

      character*16 lchemvar(NSP)
      character*180 lqkchem(NUM_K),  lqjchem(NUM_J)
!
!.... All species labels
!
      data lchemvar(1) /"A3O2"/
      data lchemvar(2) /"ACTA"/
      data lchemvar(3) /"ALD2"/
      data lchemvar(4) /"ATO2"/
      data lchemvar(5) /"B3O2"/
      data lchemvar(6) /"Br"/
      data lchemvar(7) /"BrCl"/
      data lchemvar(8) /"BrO"/
      data lchemvar(9) /"BrONO2"/
      data lchemvar(10) /"C2H6"/
      data lchemvar(11) /"C3H8"/
      data lchemvar(12) /"CCl4"/
      data lchemvar(13) /"CF2ClBr"/
      data lchemvar(14) /"CF3Br"/
      data lchemvar(15) /"CFC113"/
      data lchemvar(16) /"CFC11"/
      data lchemvar(17) /"CFC12"/
      data lchemvar(18) /"CH2O"/
      data lchemvar(19) /"CH3Br"/
      data lchemvar(20) /"CH3CCl3"/
      data lchemvar(21) /"CH3Cl"/
      data lchemvar(22) /"CH4"/
      data lchemvar(23) /"Cl2"/
      data lchemvar(24) /"Cl2O2"/
      data lchemvar(25) /"Cl"/
      data lchemvar(26) /"ClO"/
      data lchemvar(27) /"ClONO2"/
      data lchemvar(28) /"CO"/
      data lchemvar(29) /"ETO2"/
      data lchemvar(30) /"H"/
      data lchemvar(31) /"H2"/
      data lchemvar(32) /"H2O2"/
      data lchemvar(33) /"H2O"/
      data lchemvar(34) /"HBr"/
      data lchemvar(35) /"HCFC141b"/
      data lchemvar(36) /"HCFC142b"/
      data lchemvar(37) /"HCFC22"/
      data lchemvar(38) /"HCl"/
      data lchemvar(39) /"HCOOH"/
      data lchemvar(40) /"HNO2"/
      data lchemvar(41) /"HNO3"/
      data lchemvar(42) /"HNO4"/
      data lchemvar(43) /"HO2"/
      data lchemvar(44) /"HOBr"/
      data lchemvar(45) /"HOCl"/
      data lchemvar(46) /"IAO2"/
      data lchemvar(47) /"IAP"/
      data lchemvar(48) /"ISOP"/
      data lchemvar(49) /"KO2"/
      data lchemvar(50) /"MCO3"/
      data lchemvar(51) /"MEK"/
      data lchemvar(52) /"MO2"/
      data lchemvar(53) /"MP"/
      data lchemvar(54) /"N"/
      data lchemvar(55) /"N2O"/
      data lchemvar(56) /"N2O5"/
      data lchemvar(57) /"NO"/
      data lchemvar(58) /"NO2"/
      data lchemvar(59) /"NO3"/
      data lchemvar(60) /"O"/
      data lchemvar(61) /"O1D"/
      data lchemvar(62) /"O3"/
      data lchemvar(63) /"OClO"/
      data lchemvar(64) /"OH"/
      data lchemvar(65) /"PAN"/
      data lchemvar(66) /"RA3P"/
      data lchemvar(67) /"RB3P"/
      data lchemvar(68) /"RCHO"/
      data lchemvar(69) /"RCO3"/
      data lchemvar(70) /"RIO1"/
      data lchemvar(71) /"RIO2"/
      data lchemvar(72) /"RIP"/
      data lchemvar(73) /"ACET"/
      data lchemvar(74) /"N2"/
      data lchemvar(75) /"O2"/
      data lchemvar(76) /"NUMDENS"/
      data lchemvar(77) /"HNO3COND"/
!
!.... Thermal reaction labels
!
      data (lqkchem(kmg_i), kmg_i=1,10) / &
     & 'O + O2 = O3', &
     & 'O + O3 = 2 O2', &
     & 'N2 + O1D = N2 + O', &
     & 'O1D + O2 = O + O2', &
     & 'O1D + O3 = 2 O2', &
     & 'O1D + O3 = 2 O + O2', &
     & 'H2O + O1D = 2 OH', &
     & 'H2 + O1D = H + OH', &
     & 'N2O + O1D = N2 + O2', &
     & 'N2O + O1D = 2 NO' /

      data (lqkchem(kmg_i), kmg_i=11,20) / &
     & 'CH4 + O1D = MO2 + OH', &
     & 'CH4 + O1D = CH2O + H + HO2', &
     & 'CH4 + O1D = CH2O + H2', &
     & 'CFC12 + O1D = 2 Cl', &
     & 'CFC113 + O1D = 3 Cl', &
     & 'HCFC22 + O1D = Cl', &
     & 'HCFC141b + O1D = 2 Cl', &
     & 'HCFC142b + O1D = Cl', &
     & 'H + O2 = HO2', &
     & 'H + O3 = O2 + OH' /

      data (lqkchem(kmg_i), kmg_i=21,30) / &
     & 'O + OH = H + O2', &
     & 'HO2 + O = O2 + OH', &
     & 'H + HO2 = 2 OH', &
     & 'NO + O3 = NO2 + O2', &
     & 'O3 + OH = HO2 + O2', &
     & 'HO2 + O3 = 2 O2 + OH', &
     & 'NO2 + O3 = NO3 + O2', &
     & 'OH + OH = H2O + O', &
     & 'OH + OH = H2O2', &
     & 'HO2 + OH = H2O + O2' /

      data (lqkchem(kmg_i), kmg_i=31,40) / &
     & 'H2O2 + OH = H2O + HO2', &
     & 'HO2 + NO = NO2 + OH', &
     & 'HO2 + HO2 = H2O2 + O2', &
     & 'H2O + HO2 + HO2 = H2O + H2O2 + O2', &
     & 'H2 + OH = H + H2O', &
     & 'CO + OH = H', &
     & 'CH4 + OH = H2O + MO2', &
     & 'ClO + MO2 = CH2O + Cl + HO2 + O2', &
     & 'MO2 + NO = CH2O + HO2 + NO2', &
     & 'HO2 + MO2 = MP + O2' /

      data (lqkchem(kmg_i), kmg_i=41,50) / &
     & 'MO2 + MO2 = CH2O + O2', &
     & 'MO2 + MO2 = 2 CH2O + 2 HO2', &
     & 'MP + OH = H2O + MO2', &
     & 'MP + OH = CH2O + H2O + OH', &
     & 'CH2O + OH = CO + H2O + HO2', &
     & 'N + O2 = NO + O', &
     & 'N + NO = N2 + O', &
     & 'NO2 + O = NO + O2', &
     & 'NO3 + O = NO2 + O2', &
     & 'NO2 + OH = HNO3' /

      data (lqkchem(kmg_i), kmg_i=51,60) / &
     & 'HNO3 + OH = H2O + NO3', &
     & 'NO + OH = HNO2', &
     & 'HNO2 + OH = H2O + NO2', &
     & 'HO2 + NO2 = HNO4', &
     & 'HNO4 = HO2 + NO2', &
     & 'HNO4 + OH = H2O + NO2 + O2', &
     & 'HO2 + NO3 = NO2 + O2 + OH', &
     & 'NO + NO3 = 2 NO2', &
     & 'NO3 + OH = HO2 + NO2', &
     & 'NO2 + NO3 = N2O5' /

      data (lqkchem(kmg_i), kmg_i=61,70) / &
     & 'N2O5 = NO2 + NO3', &
     & 'HCOOH + OH = H2O + HO2', &
     & 'NO2 + NO3 = NO + NO2 + O2', &
     & 'CH2O + NO3 = CO + HNO3 + HO2', &
     & 'Cl + O3 = ClO + O2', &
     & 'Cl + H2 = H + HCl', &
     & 'Cl + H2O2 = HCl + HO2', &
     & 'Cl + HO2 = HCl + O2', &
     & 'Cl + HO2 = ClO + OH', &
     & 'ClO + O = Cl + O2' /

      data (lqkchem(kmg_i), kmg_i=71,80) / &
     & 'ClO + OH = Cl + HO2', &
     & 'ClO + OH = HCl + O2', &
     & 'ClO + HO2 = HOCl + O2', &
     & 'ClO + NO = Cl + NO2', &
     & 'ClO + NO2 = ClONO2', &
     & 'ClO + ClO = 2 Cl + O2', &
     & 'ClO + ClO = Cl2 + O2', &
     & 'ClO + ClO = Cl + OClO', &
     & 'ClO + ClO = Cl2O2', &
     & 'Cl2O2 = 2 ClO' /

      data (lqkchem(kmg_i), kmg_i=81,90) / &
     & 'HCl + OH = Cl + H2O', &
     & 'HOCl + OH = ClO + H2O', &
     & 'ClONO2 + O = ClO + NO3', &
     & 'ClONO2 + OH = HOCl + NO3', &
     & 'Cl + ClONO2 = Cl2 + NO3', &
     & 'Br + O3 = BrO + O2', &
     & 'Br + HO2 = HBr + O2', &
     & 'Br + CH2O = CO + HBr + HO2', &
     & 'BrO + O = Br + O2', &
     & 'BrO + HO2 = HOBr + O2' /

      data (lqkchem(kmg_i), kmg_i=91,100) / &
     & 'BrO + NO = Br + NO2', &
     & 'BrO + NO2 = BrONO2', &
     & 'BrO + ClO = Br + OClO', &
     & 'BrO + ClO = Br + Cl + O2', &
     & 'BrO + ClO = BrCl + O2', &
     & 'BrO + BrO = 2 Br + O2', &
     & 'HBr + OH = Br + H2O', &
     & 'CH2O + O = CO + HO2 + OH', &
     & 'CH4 + Cl = HCl + MO2', &
     & 'CH2O + Cl = CO + HCl + HO2' /

      data (lqkchem(kmg_i), kmg_i=101,110) / &
     & 'CH3Cl + OH = Cl + H2O + HO2', &
     & 'CH3CCl3 + OH = 3 Cl + H2O', &
     & 'HCFC22 + OH = Cl + H2O', &
     & 'HCFC141b + OH = 2 Cl + H2O', &
     & 'HCFC142b + OH = Cl + H2O', &
     & 'CH3Cl + Cl = CO + 2 HCl + HO2', &
     & 'CH3Br + OH = Br + H2O + HO2', &
     & 'ALD2 + OH = H2O + MCO3', &
     & 'ALD2 + NO3 = HNO3 + MCO3', &
     & 'MCO3 + NO2 = PAN' /

      data (lqkchem(kmg_i), kmg_i=111,120) / &
     & 'PAN = MCO3 + NO2', &
     & 'MCO3 + NO = MO2 + NO2', &
     & 'C2H6 + OH = ETO2 + H2O', &
     & 'C2H6 + Cl = ETO2 + HCl', &
     & 'ETO2 + NO = ALD2 + HO2 + NO2', &
     & 'C3H8 + OH = B3O2', &
     & 'C3H8 + OH = A3O2', &
     & 'A3O2 + NO = HO2 + NO2 + RCHO', &
     & 'ATO2 + NO =  0.19 CH2O +  0.77 HO2 +  0.19 MCO3 +  0.96 NO2', &
     & 'KO2 + NO =  0.93 ALD2 +  0.93 MCO3 +  0.93 NO2' /

      data (lqkchem(kmg_i), kmg_i=121,130) / &
     & 'NO + RIO2 =  0.69 CH2O +  0.86 HO2 + NO2 +  0.14 RIO1', &
     & 'NO + RIO2 = HNO3', &
     & 'NO + RIO1 =  0.75 CH2O + HO2 + NO2', &
     & 'NO + RIO1 = HNO3', &
     & 'IAO2 + NO =  0.35 CH2O +  0.27 CO +  0.08 HNO3 +  0.92 HO2 +  0.92 NO2', &
     & 'B3O2 + NO = ACET + HO2 + NO2', &
     & 'ACTA + OH = H2O + MO2', &
     & 'OH + RCHO = H2O + RCO3', &
     & 'NO + RCO3 = ETO2 + NO2', &
     & 'NO3 + RCHO = HNO3 + RCO3' /

      data (lqkchem(kmg_i), kmg_i=131,140) / &
     & 'ACET + OH = ATO2 + H2O', &
     & 'A3O2 + MO2 =  0.75 CH2O + HO2 +  0.75 RCHO', &
     & 'ATO2 + HO2 = MCO3 + MO2', &
     & 'HO2 + RIO2 = RIP', &
     & 'HO2 + RIO1 = RIP', &
     & 'HO2 + IAO2 = IAP', &
     & 'B3O2 + HO2 = RB3P', &
     & 'MEK + OH = H2O + KO2', &
     & 'ETO2 + MO2 =  0.75 ALD2 +  0.75 CH2O + HO2', &
     & 'MEK + NO3 = HNO3 + KO2' /

      data (lqkchem(kmg_i), kmg_i=141,150) / &
     & 'ATO2 + MO2 =  0.85 CH2O +  0.90 HO2 +  0.10 MCO3 +  0.25 MEK', &
     & 'KO2 + MO2 =  0.50 ALD2 +  0.75 CH2O +  0.50 HO2 +  0.50 MCO3 +  0.25 MEK', &
     & 'MO2 + RIO2 =  1.10 CH2O +  0.93 HO2 +  0.25 MEK +  0.07 RIO1', &
     & 'MO2 + RIO1 =  1.13 CH2O + HO2 +  0.25 MEK', &
     & 'IAO2 + MO2 =  0.95 CH2O +  0.15 CO + HO2 +  0.25 MEK', &
     & 'B3O2 + MO2 =  0.75 ACET +  0.75 CH2O + HO2', &
     & 'ETO2 + ETO2 = ALD2', &
     & 'ETO2 + ETO2 = 2 ALD2 + 2 HO2', &
     & 'A3O2 + HO2 = RA3P', &
     & 'HO2 + MCO3 = ACTA + O3' /

      data (lqkchem(kmg_i), kmg_i=151,160) / &
     & 'ISOP + OH = RIO2', &
     & 'A3O2 + MCO3 = HO2 + MO2 + RCHO', &
     & 'A3O2 + MCO3 = ACTA + RCHO', &
     & 'ISOP + O3 =  0.90 CH2O +  0.05 CO +  0.06 HO2 +  0.10 O3 +  0.27 OH', &
     & 'MO2 + RCO3 = CH2O + ETO2 + HO2', &
     & 'MO2 + RCO3 = CH2O', &
     & 'OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO', &
     & 'OH + RB3P =  0.50 B3O2 +  0.50 OH +  0.50 RCHO', &
     & 'OH + RIP =  0.50 IAO2 +  0.10 RIO1 +  0.40 RIO2', &
     & 'IAP + OH =  0.50 IAO2 +  0.50 OH +  0.50 RCHO' /

      data (lqkchem(kmg_i), kmg_i=161,170) / &
     & 'C2H6 + NO3 = ETO2 + HNO3', &
     & 'MCO3 + MCO3 = 2 MO2', &
     & 'MCO3 + MO2 = CH2O + HO2 + MO2', &
     & 'MCO3 + MO2 = ACTA + CH2O', &
     & 'KO2 + MCO3 = ALD2 + MCO3 + MO2', &
     & 'ATO2 + MCO3 =  0.20 CH2O +  0.80 HO2 +  0.20 MCO3 + MO2', &
     & 'MCO3 + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.14 RIO1', &
     & 'MCO3 + RIO1 =  0.75 CH2O + HO2 + MO2', &
     & 'IAO2 + MCO3 =  0.40 CH2O +  0.29 CO + HO2 + MO2', &
     & 'B3O2 + MCO3 = ACET + HO2 + MO2' /

      data (lqkchem(kmg_i), kmg_i=171,180) / &
     & 'ATO2 + MCO3 = ACTA + MEK', &
     & 'KO2 + MCO3 = ACTA + MEK', &
     & 'MCO3 + RIO2 = ACTA + MEK', &
     & 'MCO3 + RIO1 = ACTA + MEK', &
     & 'IAO2 + MCO3 = ACTA + MEK', &
     & 'B3O2 + MCO3 = ACET + ACTA', &
     & 'ETO2 + MCO3 = ALD2 + HO2 + MO2', &
     & 'ETO2 + MCO3 = ACTA + ALD2', &
     & 'MCO3 + RCO3 = ETO2 + MO2', &
     & 'N2O5 = 2 HNO3' /

      data (lqkchem(kmg_i), kmg_i=181,190) / &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3', &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'N2O5 = 2 HNO3', &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3', &
     & 'HCl + HOCl = Cl2 + H2O' /

      data (lqkchem(kmg_i), kmg_i=191,200) / &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3', &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'BrONO2 + HCl = BrCl + HNO3', &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3' /

      data (lqkchem(kmg_i), kmg_i=201,209) / &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'BrONO2 + HCl = BrCl + HNO3', &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'HNO3 = NO2 + OH', &
     & 'NO3 + NO3 = 2 NO2 + O2', &
     & 'HO2 =  0.50 H2O2', &
     & 'NO2 =  0.50 HNO2 +  0.50 HNO3', &
     & 'NO3 = HNO3', &
     & 'N2O5 = 2 HNO3' /
!
!.... Photolytic reaction labels
!
      data (lqjchem(kmg_i), kmg_i=1,10) / &
     & 'O2 + hv = 2 O', &
     & 'O3 + hv = O + O2', &
     & 'O3 + hv = O1D + O2', &
     & 'NO + hv = N + O', &
     & 'N2O + hv = N2 + O1D', &
     & 'NO2 + hv = NO + O', &
     & 'H2O2 + hv = 2 OH', &
     & 'MP + hv = CH2O + HO2 + OH', &
     & 'CH2O + hv = CO + H + HO2', &
     & 'CH2O + hv = CO + H2' /

      data (lqjchem(kmg_i), kmg_i=11,20) / &
     & 'HNO3 + hv = NO2 + OH', &
     & 'HNO2 + hv = NO + OH', &
     & 'HNO4 + hv = NO3 + OH', &
     & 'HNO4 + hv = HO2 + NO2', &
     & 'NO3 + hv = NO2 + O3', &
     & 'NO3 + hv = NO + O2', &
     & 'N2O5 + hv = NO2 + NO3', &
     & 'N2O5 + hv = NO + NO3 + O3', &
     & 'Cl2 + hv = 2 Cl', &
     & 'OClO + hv = ClO + O' /

      data (lqjchem(kmg_i), kmg_i=21,30) / &
     & 'Cl2O2 + hv = 2 Cl + O2', &
     & 'HOCl + hv = Cl + OH', &
     & 'ClONO2 + hv = Cl + NO3', &
     & 'ClONO2 + hv = ClO + NO2', &
     & 'BrCl + hv = Br + Cl', &
     & 'BrO + hv = Br + O', &
     & 'HOBr + hv = Br + OH', &
     & 'BrONO2 + hv = Br + NO3', &
     & 'BrONO2 + hv = BrO + NO2', &
     & 'CH3Cl + hv = Cl + MO2' /

      data (lqjchem(kmg_i), kmg_i=31,40) / &
     & 'CCl4 + hv = 4 Cl', &
     & 'CH3CCl3 + hv = 3 Cl', &
     & 'CFC11 + hv = 3 Cl', &
     & 'CFC12 + hv = 2 Cl', &
     & 'CFC113 + hv = 3 Cl', &
     & 'HCFC141b + hv = 2 Cl', &
     & 'HCFC142b + hv = Cl', &
     & 'CH3Br + hv = Br + MO2', &
     & 'CF3Br + hv = Br', &
     & 'CF2ClBr + hv = Br + Cl' /

      data (lqjchem(kmg_i), kmg_i=41,49) / &
     & 'ALD2 + hv = CO + HO2 + MO2', &
     & 'ALD2 + hv = CH4 + CO', &
     & 'PAN + hv = MCO3 + NO2', &
     & 'RCHO + hv = CO + ETO2 + HO2', &
     & 'ACET + hv = MCO3 + MO2', &
     & 'RA3P + hv = HO2 + OH + RCHO', &
     & 'RB3P + hv = HO2 + OH + RCHO', &
     & 'RIP + hv =  0.69 CH2O +  0.86 HO2 + OH +  0.14 RIO1', &
     & 'IAP + hv =  0.67 CO + HO2 + OH' /

!                                  --^--

