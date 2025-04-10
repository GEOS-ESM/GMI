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
!  Input mechanism:        GeosCCM_Combo_2020_HFC_S_VSLCL_JPL19.txt
!  Reaction dictionary:    GMI_reactions_JPL19.db
!  Setkin files generated: Tue Oct 22 16:51:18 2024
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
      data lchemvar(4) /"ALK4"/
      data lchemvar(5) /"ATO2"/
      data lchemvar(6) /"B3O2"/
      data lchemvar(7) /"Br"/
      data lchemvar(8) /"BrCl"/
      data lchemvar(9) /"BrO"/
      data lchemvar(10) /"BrONO2"/
      data lchemvar(11) /"C2Cl4"/
      data lchemvar(12) /"C2H4Cl2"/
      data lchemvar(13) /"C2H6"/
      data lchemvar(14) /"C3H8"/
      data lchemvar(15) /"CCl4"/
      data lchemvar(16) /"CF2Br2"/
      data lchemvar(17) /"CF2ClBr"/
      data lchemvar(18) /"CF3Br"/
      data lchemvar(19) /"CFC113"/
      data lchemvar(20) /"CFC114"/
      data lchemvar(21) /"CFC115"/
      data lchemvar(22) /"CFC11"/
      data lchemvar(23) /"CFC12"/
      data lchemvar(24) /"CHCl3"/
      data lchemvar(25) /"CH2Br2"/
      data lchemvar(26) /"CH2Cl2"/
      data lchemvar(27) /"CH2O"/
      data lchemvar(28) /"CH3Br"/
      data lchemvar(29) /"CH3CCl3"/
      data lchemvar(30) /"CH3Cl"/
      data lchemvar(31) /"CH4"/
      data lchemvar(32) /"CHBr3"/
      data lchemvar(33) /"Cl2"/
      data lchemvar(34) /"Cl2O2"/
      data lchemvar(35) /"Cl"/
      data lchemvar(36) /"ClO"/
      data lchemvar(37) /"ClONO2"/
      data lchemvar(38) /"CO"/
      data lchemvar(39) /"DMS"/
      data lchemvar(40) /"EOH"/
      data lchemvar(41) /"ETO2"/
      data lchemvar(42) /"ETP"/
      data lchemvar(43) /"GLYC"/
      data lchemvar(44) /"GLYX"/
      data lchemvar(45) /"H"/
      data lchemvar(46) /"H2"/
      data lchemvar(47) /"H2402"/
      data lchemvar(48) /"H2O2"/
      data lchemvar(49) /"H2O"/
      data lchemvar(50) /"HAC"/
      data lchemvar(51) /"HBr"/
      data lchemvar(52) /"HCFC141b"/
      data lchemvar(53) /"HCFC142b"/
      data lchemvar(54) /"HCFC22"/
      data lchemvar(55) /"HCl"/
      data lchemvar(56) /"HCOOH"/
      data lchemvar(57) /"HFC125"/
      data lchemvar(58) /"HFC134a"/
      data lchemvar(59) /"HFC143a"/
      data lchemvar(60) /"HFC152a"/
      data lchemvar(61) /"HFC23"/
      data lchemvar(62) /"HFC32"/
      data lchemvar(63) /"HNO2"/
      data lchemvar(64) /"HNO3"/
      data lchemvar(65) /"HNO4"/
      data lchemvar(66) /"HO2"/
      data lchemvar(67) /"HOBr"/
      data lchemvar(68) /"HOCl"/
      data lchemvar(69) /"IALD"/
      data lchemvar(70) /"IAO2"/
      data lchemvar(71) /"INO2"/
      data lchemvar(72) /"INPN"/
      data lchemvar(73) /"ISOP"/
      data lchemvar(74) /"KO2"/
      data lchemvar(75) /"MACR"/
      data lchemvar(76) /"MAO3"/
      data lchemvar(77) /"MAOP"/
      data lchemvar(78) /"MAP"/
      data lchemvar(79) /"MCO3"/
      data lchemvar(80) /"MEK"/
      data lchemvar(81) /"MGLY"/
      data lchemvar(82) /"MO2"/
      data lchemvar(83) /"MOH"/
      data lchemvar(84) /"MP"/
      data lchemvar(85) /"MRP"/
      data lchemvar(86) /"MVK"/
      data lchemvar(87) /"N"/
      data lchemvar(88) /"N2O"/
      data lchemvar(89) /"N2O5"/
      data lchemvar(90) /"NO"/
      data lchemvar(91) /"NO2"/
      data lchemvar(92) /"NO3"/
      data lchemvar(93) /"O"/
      data lchemvar(94) /"O1D"/
      data lchemvar(95) /"O3"/
      data lchemvar(96) /"OClO"/
      data lchemvar(97) /"OCSg"/
      data lchemvar(98) /"OH"/
      data lchemvar(99) /"PAN"/
      data lchemvar(100) /"PMN"/
      data lchemvar(101) /"PO2"/
      data lchemvar(102) /"PP"/
      data lchemvar(103) /"PPN"/
      data lchemvar(104) /"PRN1"/
      data lchemvar(105) /"PRPE"/
      data lchemvar(106) /"PRPN"/
      data lchemvar(107) /"R4N1"/
      data lchemvar(108) /"R4N2"/
      data lchemvar(109) /"R4O2"/
      data lchemvar(110) /"R4P"/
      data lchemvar(111) /"RA3P"/
      data lchemvar(112) /"RB3P"/
      data lchemvar(113) /"RCHO"/
      data lchemvar(114) /"RCO3"/
      data lchemvar(115) /"RCOOH"/
      data lchemvar(116) /"RIO1"/
      data lchemvar(117) /"RIO2"/
      data lchemvar(118) /"RIPA"/
      data lchemvar(119) /"RIPB"/
      data lchemvar(120) /"ROH"/
      data lchemvar(121) /"RP"/
      data lchemvar(122) /"VRO2"/
      data lchemvar(123) /"VRP"/
      data lchemvar(124) /"SO2"/
      data lchemvar(125) /"H2SO4"/
      data lchemvar(126) /"ACET"/
      data lchemvar(127) /"N2"/
      data lchemvar(128) /"O2"/
      data lchemvar(129) /"NUMDENS"/
      data lchemvar(130) /"HNO3COND"/
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
     & 'CFC114 + O1D = 2 Cl', &
     & 'CFC115 + O1D = Cl', &
     & 'HCFC22 + O1D = Cl', &
     & 'HCFC141b + O1D = 2 Cl', &
     & 'HCFC142b + O1D = Cl' /

      data (lqkchem(kmg_i), kmg_i=21,30) / &
     & 'H + O2 = HO2', &
     & 'H + O3 = O2 + OH', &
     & 'O + OH = H + O2', &
     & 'HO2 + O = O2 + OH', &
     & 'H + HO2 = 2 OH', &
     & 'O3 + OH = HO2 + O2', &
     & 'HO2 + O3 = 2 O2 + OH', &
     & 'OH + OH = H2O + O', &
     & 'OH + OH = H2O2', &
     & 'HO2 + OH = H2O + O2' /

      data (lqkchem(kmg_i), kmg_i=31,40) / &
     & 'H2O2 + OH = H2O + HO2', &
     & 'HO2 + HO2 = H2O2 + O2', &
     & 'H2O + HO2 + HO2 = H2O + H2O2 + O2', &
     & 'H2 + OH = H + H2O', &
     & 'NO + O3 = NO2 + O2', &
     & 'NO2 + O3 = NO3 + O2', &
     & 'HO2 + NO = NO2 + OH', &
     & 'CO + OH = H', &
     & 'CH4 + OH = H2O + MO2', &
     & 'MO2 + NO = CH2O + HO2 + NO2' /

      data (lqkchem(kmg_i), kmg_i=41,50) / &
     & 'ClO + MO2 = CH2O + Cl + HO2 + O2', &
     & 'HO2 + MO2 = MP + O2', &
     & 'MO2 + MO2 = CH2O + MOH + O2', &
     & 'MO2 + MO2 = 2 CH2O + 2 HO2', &
     & 'MP + OH = H2O + MO2', &
     & 'MP + OH = CH2O + H2O + OH', &
     & 'CH2O + OH = CO + H2O + HO2', &
     & 'NO2 + O = NO + O2', &
     & 'NO2 + O = NO3', &
     & 'NO3 + O = NO2 + O2' /

      data (lqkchem(kmg_i), kmg_i=51,60) / &
     & 'N + O2 = NO + O', &
     & 'N + NO = N2 + O', &
     & 'N + NO2 = N2O + O', &
     & 'NO2 + OH = HNO3', &
     & 'HNO3 + OH = H2O + NO3', &
     & 'NO + OH = HNO2', &
     & 'HNO2 + OH = H2O + NO2', &
     & 'HO2 + NO2 = HNO4', &
     & 'HNO4 = HO2 + NO2', &
     & 'HNO4 + OH = H2O + NO2 + O2' /

      data (lqkchem(kmg_i), kmg_i=61,70) / &
     & 'HO2 + NO3 = NO2 + O2 + OH', &
     & 'NO + NO3 = 2 NO2', &
     & 'NO3 + OH = HO2 + NO2', &
     & 'NO2 + NO3 = N2O5', &
     & 'N2O5 = NO2 + NO3', &
     & 'HCOOH + OH = H2O + HO2', &
     & 'MOH + OH = CH2O + HO2', &
     & 'NO2 + NO3 = NO + NO2 + O2', &
     & 'CH2O + NO3 = CO + HNO3 + HO2', &
     & 'Cl + O3 = ClO + O2' /

      data (lqkchem(kmg_i), kmg_i=71,80) / &
     & 'Cl + H2 = H + HCl', &
     & 'Cl + H2O2 = HCl + HO2', &
     & 'Cl + HO2 = HCl + O2', &
     & 'Cl + HO2 = ClO + OH', &
     & 'ClO + O = Cl + O2', &
     & 'ClO + OH = Cl + HO2', &
     & 'ClO + OH = HCl + O2', &
     & 'ClO + HO2 = HOCl + O2', &
     & 'ClO + NO = Cl + NO2', &
     & 'ClO + NO2 = ClONO2' /

      data (lqkchem(kmg_i), kmg_i=81,90) / &
     & 'ClO + ClO = 2 Cl + O2', &
     & 'ClO + ClO = Cl2 + O2', &
     & 'ClO + ClO = Cl + OClO', &
     & 'ClO + ClO = Cl2O2', &
     & 'Cl2O2 = 2 ClO', &
     & 'HCl + OH = Cl + H2O', &
     & 'HOCl + OH = ClO + H2O', &
     & 'ClONO2 + O = ClO + NO3', &
     & 'ClONO2 + OH = HOCl + NO3', &
     & 'Cl + ClONO2 = Cl2 + NO3' /

      data (lqkchem(kmg_i), kmg_i=91,100) / &
     & 'Br + O3 = BrO + O2', &
     & 'Br + HO2 = HBr + O2', &
     & 'Br + CH2O = CO + HBr + HO2', &
     & 'BrO + O = Br + O2', &
     & 'BrO + HO2 = HOBr + O2', &
     & 'BrO + NO = Br + NO2', &
     & 'BrO + NO2 = BrONO2', &
     & 'BrO + ClO = Br + OClO', &
     & 'BrO + ClO = Br + Cl + O2', &
     & 'BrO + ClO = BrCl + O2' /

      data (lqkchem(kmg_i), kmg_i=101,110) / &
     & 'BrO + BrO = 2 Br + O2', &
     & 'HBr + OH = Br + H2O', &
     & 'CHBr3 + OH = 3 Br', &
     & 'CH2Br2 + OH = 2 Br', &
     & 'CH2O + O = CO + HO2 + OH', &
     & 'CH4 + Cl = HCl + MO2', &
     & 'CH2O + Cl = CO + HCl + HO2', &
     & 'CH3Cl + OH = Cl + H2O + HO2', &
     & 'CH3CCl3 + OH = 3 Cl + H2O', &
     & 'HCFC22 + OH = Cl + H2O' /

      data (lqkchem(kmg_i), kmg_i=111,120) / &
     & 'HCFC141b + OH = 2 Cl + H2O', &
     & 'HCFC142b + OH = Cl + H2O', &
     & 'CH3Cl + Cl = CO + 2 HCl + HO2', &
     & 'CH3Br + OH = Br + H2O + HO2', &
     & 'HFC23 + O1D =  0.25 H2O +  0.75 O', &
     & 'HFC32 + O1D =  0.30 H2O +  0.70 O', &
     & 'HFC125 + O1D =  0.15 H2O +  0.25 O +  0.60 OH', &
     & 'HFC134a + O1D =  0.11 H2O +  0.65 O +  0.24 OH', &
     & 'HFC143a + O1D =  0.27 H2O +  0.35 O +  0.38 OH', &
     & 'HFC152a + O1D =  0.40 H2O +  0.45 O +  0.15 OH' /

      data (lqkchem(kmg_i), kmg_i=121,130) / &
     & 'HFC23 + OH = H2O', &
     & 'HFC32 + OH = H2O', &
     & 'HFC125 + OH = H2O', &
     & 'HFC134a + OH = H2O', &
     & 'HFC143a + OH = H2O', &
     & 'HFC152a + OH = H2O', &
     & 'CHCl3 + OH = 3 Cl', &
     & 'CH2Cl2 + OH = 2 Cl', &
     & 'C2H4Cl2 + OH = 2 Cl', &
     & 'CHCl3 + Cl = 3 Cl + HCl' /

      data (lqkchem(kmg_i), kmg_i=131,140) / &
     & 'CH2Cl2 + Cl = 2 Cl + HCl', &
     & 'C2Cl4 + Cl = 2 Cl', &
     & 'A3O2 + HO2 = RA3P', &
     & 'A3O2 + MO2 =  0.75 CH2O + HO2 +  0.25 MOH +  0.75 RCHO +  0.25 ROH', &
     & 'A3O2 + NO = HO2 + NO2 + RCHO', &
     & 'ACET + OH = ATO2 + H2O', &
     & 'ACTA + OH = H2O + MO2', &
     & 'ALD2 + NO3 = HNO3 + MCO3', &
     & 'ALD2 + OH =  0.05 CH2O +  0.05 CO + H2O +  0.05 HO2 +  0.95 MCO3', &
     & 'ALK4 + NO3 = HNO3 + R4O2' /

      data (lqkchem(kmg_i), kmg_i=141,150) / &
     & 'ALK4 + OH = R4O2', &
     & 'ATO2 + HO2 =  0.15 CH2O +  0.85 HCOOH +  0.15 MCO3 +  0.15 OH', &
     & 'ATO2 + MCO3 =  0.10 ACTA +  0.90 CH2O +  0.90 MCO3 +  0.10 MGLY +  0.90 MO2', &
     & 'ATO2 + MO2 =  0.50 CH2O +  0.20 HAC +  0.30 HO2 +  0.30 MCO3 +  0.50 MGLY +  0.50 MOH', &
     & 'ATO2 + NO = CH2O + MCO3 + NO2', &
     & 'B3O2 + HO2 = RB3P', &
     & 'B3O2 + MCO3 = ACET +  0.10 ACTA +  0.90 HO2 +  0.90 MO2', &
     & 'B3O2 + MO2 =  0.75 ACET +  0.75 CH2O + HO2 +  0.25 MOH +  0.25 ROH', &
     & 'B3O2 + NO = ACET + HO2 + NO2', &
     & 'C2H6 + NO3 = ETO2 + HNO3' /

      data (lqkchem(kmg_i), kmg_i=151,160) / &
     & 'C2H6 + OH = ETO2 + H2O', &
     & 'C2H6 + Cl = ETO2 + HCl', &
     & 'C3H8 + OH = A3O2', &
     & 'C3H8 + OH = B3O2', &
     & 'EOH + OH = ALD2 + HO2', &
     & 'ETO2 + ETO2 =  1.60 ALD2 +  0.40 EOH +  1.20 HO2', &
     & 'ETO2 + NO = ALD2 + HO2 + NO2', &
     & 'ETP + OH =  0.64 ALD2 +  0.36 ETO2 +  0.64 OH', &
     & 'GLYC + OH =  0.73 CH2O +  0.50 CO +  0.13 GLYX +  0.13 HCOOH +  0.77 HO2 +  0.23 OH', &
     & 'GLYC + OH = CO + HCOOH + OH' /

      data (lqkchem(kmg_i), kmg_i=161,170) / &
     & 'GLYX + NO3 = 2 CO + HNO3 + HO2', &
     & 'GLYX + OH = 2 CO + HO2', &
     & 'HAC + OH = HO2 + MGLY', &
     & 'HAC + OH =  0.50 ACTA +  0.50 CO +  0.50 HCOOH +  0.50 MO2 + OH', &
     & 'ETO2 + HO2 = ETP', &
     & 'HO2 + MCO3 =  0.13 ACTA +  0.37 MAP +  0.50 MO2 +  0.13 O3 +  0.50 OH', &
     & 'IALD + OH =  0.15 HO2 +  0.44 IAO2 +  0.41 MAO3', &
     & 'HO2 + IAO2 =  0.50 MACR +  0.50 MVK + 2 OH', &
     & 'IAO2 + NO =  0.50 MACR +  0.50 MVK + NO2 + OH', &
     & 'HO2 + INO2 = INPN' /

      data (lqkchem(kmg_i), kmg_i=171,180) / &
     & 'INO2 + MCO3 =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR + MO2 +  0.05 MVK +  0.15 NO2', &
     & 'INO2 + MCO3 = ACTA + NO2 + RCHO', &
     & 'INO2 + MO2 =  0.83 CH2O +  0.43 HNO3 +  0.90 HO2 +  0.05 MACR +  0.25 MOH +  0.03 MVK +  0.57 NO2 +  0.25 RCHO +  0.25 ROH', &
     & 'INO2 + NO =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR +  0.05 MVK +  1.15 NO2', &
     & 'INPN + OH = INO2', &
     & 'ISOP + NO3 = INO2', &
     & 'ISOP + O3 =  0.83 CH2O +  0.41 CO +  0.01 H2O2 +  0.58 HCOOH +  0.16 HO2 +  0.42 MACR +  0.41 MO2 +  0.18 MVK +  0.28 OH', &
     & 'ISOP + OH = RIO2', &
     & 'HO2 + KO2 =  0.15 ALD2 +  0.85 HCOOH +  0.15 MCO3 +  0.85 MO2 +  0.15 OH', &
     & 'KO2 + MCO3 =  0.10 ACTA +  0.90 ALD2 +  0.90 MCO3 +  0.10 MEK +  0.90 MO2' /

      data (lqkchem(kmg_i), kmg_i=181,190) / &
     & 'KO2 + MO2 =  0.50 ALD2 +  0.75 CH2O +  0.50 HO2 +  0.50 MCO3 +  0.25 MEK +  0.25 MO2 +  0.25 MOH +  0.25 ROH', &
     & 'KO2 + NO =  0.92 ALD2 +  0.92 MCO3 +  0.93 NO2 +  0.07 R4N2', &
     & 'MACR + NO3 =  0.68 CO +  0.32 HNO3 +  0.32 MAO3 +  0.68 MGLY +  0.68 NO2', &
     & 'MACR + OH = MAO3', &
     & 'MACR + O3 =  0.12 CH2O +  0.12 CO +  0.88 HCOOH +  0.12 MCO3 +  0.88 MGLY +  0.12 OH', &
     & 'HO2 + MAO3 =  0.50 CH2O +  0.32 CO +  0.50 MAOP +  0.17 MCO3 +  0.32 MO2 +  0.13 O3 +  0.50 OH', &
     & 'MAO3 + NO = CH2O +  0.65 CO +  0.35 MCO3 +  0.65 MO2 + NO2', &
     & 'MAO3 + NO2 = PMN', &
     & 'MAOP + OH =  0.25 CH2O +  0.49 CO +  0.49 HAC +  0.17 MAO3 +  0.10 MAOP +  0.09 MCO3 +  0.16 MO2 +  0.58 OH', &
     & 'A3O2 + MCO3 =  0.10 ACTA +  0.90 HO2 +  0.90 MO2 + RCHO' /

      data (lqkchem(kmg_i), kmg_i=191,200) / &
     & 'ETO2 + MCO3 =  0.10 ACTA + ALD2 +  0.90 HO2 +  0.90 MO2', &
     & 'MCO3 + MCO3 = 2 MO2', &
     & 'MCO3 + MO2 =  0.10 ACTA + CH2O +  0.90 HO2 +  0.90 MO2', &
     & 'MCO3 + NO2 = PAN', &
     & 'MCO3 + NO = MO2 + NO2', &
     & 'MCO3 + PO2 =  0.10 ACTA +  0.90 ALD2 +  0.90 CH2O +  0.06 HAC +  0.90 HO2 +  0.90 MO2 +  0.04 RCHO', &
     & 'MEK + NO3 = HNO3 + KO2', &
     & 'MEK + OH = H2O + KO2', &
     & 'MGLY + NO3 = CO + HNO3 + MCO3', &
     & 'MGLY + OH = CO + MCO3' /

      data (lqkchem(kmg_i), kmg_i=201,210) / &
     & 'MVK + O3 =  0.04 ALD2 +  0.80 CH2O +  0.05 CO +  0.06 HO2 +  0.82 MGLY +  0.20 O3 +  0.08 OH', &
     & 'MVK + OH = VRO2', &
     & 'MAP + OH =  0.22 CH2O +  0.78 MCO3 +  0.22 OH', &
     & 'OH + RCHO = H2O + RCO3', &
     & 'OH + RCOOH = ETO2', &
     & 'PAN = MCO3 + NO2', &
     & 'PMN = MAO3 + NO2', &
     & 'OH + PMN =  0.25 CO +  0.25 HAC +  0.75 MAOP + NO3', &
     & 'HO2 + PO2 = PP', &
     & 'MO2 + PO2 =  0.50 ALD2 +  1.25 CH2O +  0.16 HAC + HO2 +  0.25 MOH +  0.09 RCHO +  0.25 ROH' /

      data (lqkchem(kmg_i), kmg_i=211,220) / &
     & 'NO + PO2 = ALD2 + CH2O + HO2 + NO2', &
     & 'PPN = NO2 + RCO3', &
     & 'OH + PP =  0.79 HAC +  0.79 OH +  0.21 PO2', &
     & 'HO2 + PRN1 = PRPN', &
     & 'MCO3 + PRN1 =  0.10 ACTA +  0.90 ALD2 +  0.90 CH2O +  0.90 MO2 + NO2 +  0.10 RCHO', &
     & 'MO2 + PRN1 =  0.50 ALD2 +  1.25 CH2O +  0.50 HO2 +  0.25 MOH + NO2 +  0.25 RCHO +  0.25 ROH', &
     & 'NO + PRN1 = ALD2 + CH2O + 2 NO2', &
     & 'NO3 + PRPE = PRN1', &
     & 'O3 + PRPE =  0.12 ACTA +  0.50 ALD2 +  0.50 CH2O +  0.10 CH4 +  0.56 CO +  0.22 HCOOH +  0.28 HO2 +  0.28 MO2 +  0.36 OH', &
     & 'OH + PRPE = PO2' /

      data (lqkchem(kmg_i), kmg_i=221,230) / &
     & 'OH + PRPN =  0.79 MGLY +  0.79 NO2 +  0.21 PRN1', &
     & 'HO2 + R4N1 = R4N2', &
     & 'MCO3 + R4N1 =  0.10 ACTA +  0.68 ALD2 +  0.35 CH2O +  0.90 MO2 + NO2 +  0.27 R4O2 +  0.61 RCHO', &
     & 'MO2 + R4N1 =  0.38 ALD2 +  0.95 CH2O +  0.50 HO2 +  0.25 MOH + NO2 +  0.15 R4O2 +  0.58 RCHO +  0.38 ROH', &
     & 'NO + R4N1 =  0.97 ALD2 +  0.64 CH2O + 2 NO2 +  0.64 RCHO', &
     & 'OH + R4N2 = H2O + R4N1', &
     & 'HO2 + R4O2 = R4P', &
     & 'MCO3 + R4O2 =  0.05 A3O2 +  0.29 ACET +  0.10 ACTA +  0.29 ALD2 +  0.16 B3O2 +  0.29 ETO2 +  0.24 HO2 +  0.27 MEK +  0.90 MO2 +  0.25 RCHO', &
     & 'MO2 + R4O2 =  0.03 A3O2 +  0.16 ACET +  0.16 ALD2 +  0.09 B3O2 +  0.75 CH2O +  0.16 ETO2 +  0.64 HO2 +  0.35 MEK +  0.09 MO2 +  0.25 MOH +  0.13 RCHO +  0.38 ROH', &
     & 'NO + R4O2 =  0.05 A3O2 +  0.34 ACET +  0.34 ALD2 +  0.19 B3O2 +  0.34 ETO2 +  0.27 HO2 +  0.19 MEK +  0.19 MO2 + NO2 +  0.15 RCHO' /

      data (lqkchem(kmg_i), kmg_i=231,240) / &
     & 'NO + R4O2 = R4N2', &
     & 'OH + R4P =  0.79 OH +  0.21 R4O2 +  1.18 RCHO', &
     & 'OH + RA3P =  0.36 A3O2 +  0.64 OH +  0.64 RCHO', &
     & 'OH + RB3P =  0.79 ACET +  0.21 B3O2 +  0.79 OH', &
     & 'NO3 + RCHO = HNO3 + RCO3', &
     & 'HO2 + RCO3 =  0.03 A3O2 +  0.12 B3O2 +  0.22 ETO2 +  0.15 O3 +  0.44 OH +  0.15 RCOOH +  0.41 RP', &
     & 'MCO3 + RCO3 =  0.07 A3O2 +  0.27 B3O2 +  0.49 ETO2 + MO2', &
     & 'MO2 + RCO3 =  0.07 A3O2 +  0.27 B3O2 + CH2O +  0.49 ETO2 + HO2', &
     & 'NO2 + RCO3 = PPN', &
     & 'NO + RCO3 =  0.07 A3O2 +  0.27 B3O2 +  0.49 ETO2 + NO2' /

      data (lqkchem(kmg_i), kmg_i=241,250) / &
     & 'HO2 + RIO1 =  0.06 CH2O +  0.06 HO2 +  0.06 MVK +  0.06 OH +  0.94 RIPA', &
     & 'HO2 + RIO2 =  0.06 CH2O +  0.06 HO2 +  0.06 MACR +  0.06 OH +  0.94 RIPB', &
     & 'MO2 + RIO1 = 2 CH2O + 2 HO2 + MVK', &
     & 'MO2 + RIO2 = 2 CH2O + 2 HO2 + MACR', &
     & 'RIO1 + RIO1 = 2 CH2O + 2 HO2 + 2 MVK', &
     & 'RIO2 + RIO2 = 2 CH2O + 2 HO2 + 2 MACR', &
     & 'RIO1 + RIO2 = 2 CH2O + 2 HO2 + MACR + MVK', &
     & 'NO + RIO1 = CH2O + HO2 + MVK + NO2', &
     & 'NO + RIO2 = CH2O + HO2 + MACR + NO2', &
     & 'NO + RIO1 = HNO3' /

      data (lqkchem(kmg_i), kmg_i=251,260) / &
     & 'NO + RIO2 = HNO3', &
     & 'RIO1 = CH2O + MVK + OH', &
     & 'RIO2 = CH2O + MACR + OH', &
     & 'RIO1 =  0.30 CH2O +  0.60 CO +  0.40 HO2 +  0.40 IALD +  0.30 MCO3 +  0.30 MGLY +  1.50 OH', &
     & 'RIO2 =  0.30 CH2O +  0.90 CO +  0.30 HCOOH +  0.70 HO2 +  0.40 IALD +  0.30 MGLY +  1.50 OH', &
     & 'OH + RIPA =  0.25 CO +  0.25 HO2 +  0.12 MRP +  0.12 MVK +  0.75 RIO1', &
     & 'OH + RIPB =  0.33 CO +  0.33 HO2 +  0.16 IAO2 +  0.17 MACR +  0.17 MRP +  0.51 RIO2', &
     & 'OH + ROH = HO2 + RCHO', &
     & 'OH + RP =  0.50 ALD2 +  0.50 OH +  0.50 RCO3', &
     & 'HO2 + VRO2 =  0.05 CH2O +  0.36 GLYC +  0.31 HO2 +  0.25 MAOP +  0.36 MCO3 +  0.05 MGLY +  0.67 OH +  0.34 VRP' /

      data (lqkchem(kmg_i), kmg_i=261,270) / &
     & 'NO + VRO2 = HNO3', &
     & 'NO + VRO2 =  0.24 CH2O +  0.76 GLYC +  0.24 HO2 +  0.76 MCO3 +  0.24 MGLY + NO2', &
     & 'OH + VRP =  1.19 CO +  0.53 HO2 +  0.53 MCO3 +  0.19 MGLY +  0.19 OH', &
     & 'N2O5 = 2 HNO3', &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3', &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'N2O5 = 2 HNO3' /

      data (lqkchem(kmg_i), kmg_i=271,280) / &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3', &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3', &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'BrONO2 + HCl = BrCl + HNO3' /

      data (lqkchem(kmg_i), kmg_i=281,290) / &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3', &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'BrONO2 + HCl = BrCl + HNO3', &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'HNO3 = NO2 + OH', &
     & 'NO3 + NO3 = 2 NO2 + O2', &
     & 'HO2 =  0.50 H2O' /

      data (lqkchem(kmg_i), kmg_i=291,300) / &
     & 'NO2 =  0.50 HNO2 +  0.50 HNO3', &
     & 'NO3 = HNO3', &
     & 'N2O5 = 2 HNO3', &
     & 'DMS + OH = O2 + SO2', &
     & 'DMS + NO3 = HNO3 + SO2', &
     & 'O + SO2 = H2SO4', &
     & 'OH + SO2 = H2SO4', &
     & 'H2O2 + SO2 = H2SO4', &
     & 'O3 + SO2 = H2SO4', &
     & 'O + OCSg = CO + SO2' /

      data (lqkchem(kmg_i), kmg_i=301,301) / &
     & 'OCSg + OH = SO2' /
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
     & 'Cl2 + hv = 2 Cl', &
     & 'OClO + hv = ClO + O', &
     & 'Cl2O2 + hv = 2 Cl + O2' /

      data (lqjchem(kmg_i), kmg_i=21,30) / &
     & 'HOCl + hv = Cl + OH', &
     & 'ClONO2 + hv = Cl + NO3', &
     & 'ClONO2 + hv = ClO + NO2', &
     & 'BrCl + hv = Br + Cl', &
     & 'BrO + hv = Br + O', &
     & 'HOBr + hv = Br + OH', &
     & 'BrONO2 + hv = Br + NO3', &
     & 'CHBr3 + hv = 3 Br', &
     & 'CH2Br2 + hv = 2 Br', &
     & 'CH3Cl + hv = Cl + MO2' /

      data (lqjchem(kmg_i), kmg_i=31,40) / &
     & 'CCl4 + hv = 4 Cl', &
     & 'CH3CCl3 + hv = 3 Cl', &
     & 'CFC11 + hv = 3 Cl', &
     & 'CFC12 + hv = 2 Cl', &
     & 'CFC113 + hv = 3 Cl', &
     & 'CFC114 + hv = 2 Cl', &
     & 'CFC115 + hv = Cl', &
     & 'HCFC141b + hv = 2 Cl', &
     & 'HCFC142b + hv = Cl', &
     & 'CH3Br + hv = Br + MO2' /

      data (lqjchem(kmg_i), kmg_i=41,50) / &
     & 'CF3Br + hv = Br', &
     & 'CF2Br2 + hv = 2 Br', &
     & 'H2402 + hv = 2 Br', &
     & 'CF2ClBr + hv = Br + Cl', &
     & 'ALD2 + hv = CO + HO2 + MO2', &
     & 'PAN + hv =  0.70 MCO3 +  0.30 MO2 +  0.70 NO2 +  0.30 NO3', &
     & 'RCHO + hv =  0.07 A3O2 +  0.27 B3O2 + CO +  0.49 ETO2 + HO2', &
     & 'ACET + hv = MCO3 + MO2', &
     & 'MEK + hv =  0.06 A3O2 +  0.23 B3O2 +  0.41 ETO2 +  0.85 MCO3 +  0.15 MO2 +  0.15 RCO3', &
     & 'GLYC + hv =  0.90 CH2O + CO +  1.73 HO2 +  0.10 MOH +  0.07 OH' /

      data (lqjchem(kmg_i), kmg_i=51,60) / &
     & 'GLYX + hv = 2 CO + H2', &
     & 'GLYX + hv = 2 CO + 2 HO2', &
     & 'GLYX + hv = CH2O + CO', &
     & 'MGLY + hv = CO + HO2 + MCO3', &
     & 'MVK + hv = CO + PRPE', &
     & 'MVK + hv = MO2 + RCO3', &
     & 'MACR + hv = CH2O + CO + HO2 + MCO3', &
     & 'HAC + hv = CH2O + HO2 + MCO3', &
     & 'INPN + hv = HO2 + NO2 + OH + RCHO', &
     & 'PRPN + hv = HO2 + NO2 + OH + RCHO' /

      data (lqjchem(kmg_i), kmg_i=61,70) / &
     & 'ETP + hv = ALD2 + HO2 + OH', &
     & 'RA3P + hv = HO2 + OH + RCHO', &
     & 'RB3P + hv = ACET + HO2 + OH', &
     & 'R4P + hv = HO2 + OH +  1.50 RCHO', &
     & 'PP + hv = ALD2 + CH2O + HO2 + OH', &
     & 'RP + hv = ALD2 + HO2 + MO2 + OH', &
     & 'RIPA + hv = CH2O + HO2 + MVK + OH', &
     & 'RIPB + hv = CH2O + HO2 + MACR + OH', &
     & 'VRP + hv =  0.53 CH2O +  0.76 CO +  0.23 GLYX +  1.24 HO2 +  0.76 MCO3 +  0.23 MGLY', &
     & 'MRP + hv =  0.35 CH2O +  0.38 CO +  0.26 GLYC +  0.38 HAC +  0.74 HO2 +  0.26 MCO3 +  0.35 MGLY + OH' /

      data (lqjchem(kmg_i), kmg_i=71,76) / &
     & 'MAOP + hv =  0.23 CH2O +  0.77 CO +  0.77 HAC + HO2 +  0.23 MGLY + OH', &
     & 'R4N2 + hv =  0.05 A3O2 +  0.34 ACET +  0.34 ALD2 +  0.19 B3O2 +  0.34 ETO2 +  0.27 HO2 +  0.19 MEK +  0.19 MO2 + NO2 +  0.15 RCHO', &
     & 'MAP + hv = MO2 + OH', &
     & 'CH2Cl2 + hv = 2 Cl', &
     & 'OCSg + hv = CO + SO2', &
     & 'H2SO4 + hv = 2 OH + SO2' /

!                                  --^--

