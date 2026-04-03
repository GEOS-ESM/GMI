!=======================================================================
!
! $Id: $
!
! FILE
!   setkin_par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Input mechanism:        Standard.txt
!  Reaction dictionary:    GMI_reactions_JPL19.db
!  Setkin files generated: Tue Mar 24 19:30:27 2026
!
!========1=========2=========3=========4=========5=========6=========7==


      integer &
     &  NSP &
     & ,NUM_J &
     & ,NUM_K

      parameter (NSP   = 126)
      parameter (NUM_J =  81)
      parameter (NUM_K = 303)

!.... Species affiliations and assignments

      integer &
     &  NACT &
     & ,NCHEM &
     & ,NCONST &
     & ,NFAM &
     & ,NMF &
     & ,NSS

      parameter (NACT   = 121)
      parameter (NCHEM  = 125)
      parameter (NCONST =   4)
      parameter (NFAM   =   2)
      parameter (NMF    = 122)
      parameter (NSS    =   0)

!.... Ancillary input tables

      integer &
     &  nSAD, nSADdust, nSADaer

      parameter (nSAD  =   5)
      parameter (nSADdust =   7)
      parameter (nSADaer =   6)


!.... Species individual identifications

      integer IBBC
      parameter (IBBC      =   0)
      integer IBOC
      parameter (IBOC      =   0)
      integer IDUST1
      parameter (IDUST1    =   0)
      integer IDUST2
      parameter (IDUST2    =   0)
      integer IDUST3
      parameter (IDUST3    =   0)
      integer IDUST4
      parameter (IDUST4    =   0)
      integer IDUST5
      parameter (IDUST5    =   0)
      integer IFBC
      parameter (IFBC      =   0)
      integer IFOC
      parameter (IFOC      =   0)
      integer IFSO2
      parameter (IFSO2     =   0)
      integer IFSO4A
      parameter (IFSO4A    =   0)
      integer IFSO4D1
      parameter (IFSO4D1   =   0)
      integer IFSO4D2
      parameter (IFSO4D2   =   0)
      integer IFSO4D3
      parameter (IFSO4D3   =   0)
      integer IFSO4D4
      parameter (IFSO4D4   =   0)
      integer IFSO4N1
      parameter (IFSO4N1   =   0)
      integer IFSO4N2
      parameter (IFSO4N2   =   0)
      integer IFSO4N3
      parameter (IFSO4N3   =   0)
      integer IFSO4N4
      parameter (IFSO4N4   =   0)
      integer IFSO4S1
      parameter (IFSO4S1   =   0)
      integer IFSO4S2
      parameter (IFSO4S2   =   0)
      integer IFSO4S3
      parameter (IFSO4S3   =   0)
      integer IFSO4S4
      parameter (IFSO4S4   =   0)
      integer INDMS
      parameter (INDMS     =   0)
      integer INOC
      parameter (INOC      =   0)
      integer INSO2
      parameter (INSO2     =   0)
      integer INSO4A
      parameter (INSO4A    =   0)
      integer INSO4D1
      parameter (INSO4D1   =   0)
      integer INSO4D2
      parameter (INSO4D2   =   0)
      integer INSO4D3
      parameter (INSO4D3   =   0)
      integer INSO4D4
      parameter (INSO4D4   =   0)
      integer INSO4N1
      parameter (INSO4N1   =   0)
      integer INSO4N2
      parameter (INSO4N2   =   0)
      integer INSO4N3
      parameter (INSO4N3   =   0)
      integer INSO4N4
      parameter (INSO4N4   =   0)
      integer INSO4S1
      parameter (INSO4S1   =   0)
      integer INSO4S2
      parameter (INSO4S2   =   0)
      integer INSO4S3
      parameter (INSO4S3   =   0)
      integer INSO4S4
      parameter (INSO4S4   =   0)
      integer ISSLT1
      parameter (ISSLT1    =   0)
      integer ISSLT2
      parameter (ISSLT2    =   0)
      integer ISSLT3
      parameter (ISSLT3    =   0)
      integer ISSLT4
      parameter (ISSLT4    =   0)
      integer IALD2
      parameter (IALD2     =   3)
      integer IALK4
      parameter (IALK4     =   4)
      integer IBRCL
      parameter (IBRCL     =   8)
      integer IBRONO2
      parameter (IBRONO2   =  10)
      integer IC2H6
      parameter (IC2H6     =  13)
      integer IC3H8
      parameter (IC3H8     =  14)
      integer ICH2O
      parameter (ICH2O     =  27)
      integer ICH4
      parameter (ICH4      =  32)
      integer ICL2
      parameter (ICL2      =  34)
      integer ICL2O2
      parameter (ICL2O2    =  35)
      integer ICLO
      parameter (ICLO      =  37)
      integer ICLONO2
      parameter (ICLONO2   =  38)
      integer ICO
      parameter (ICO       =  39)
      integer IH2
      parameter (IH2       =  47)
      integer IH2O2
      parameter (IH2O2     =  49)
      integer IH2O
      parameter (IH2O      =  50)
      integer IHCL
      parameter (IHCL      =  56)
      integer IHNO2
      parameter (IHNO2     =  59)
      integer IHNO3
      parameter (IHNO3     =  60)
      integer IHNO4
      parameter (IHNO4     =  61)
      integer IHO2
      parameter (IHO2      =  62)
      integer IHOBR
      parameter (IHOBR     =  63)
      integer IHOCL
      parameter (IHOCL     =  64)
      integer IC5H8
      parameter (IC5H8     =  69)
      integer IMACR
      parameter (IMACR     =  71)
      integer IMCO3
      parameter (IMCO3     =  75)
      integer IMEK
      parameter (IMEK      =  76)
      integer IMP
      parameter (IMP       =  80)
      integer IMVK
      parameter (IMVK      =  82)
      integer IAN
      parameter (IAN       =  83)
      integer IN2O
      parameter (IN2O      =  84)
      integer IN2O5
      parameter (IN2O5     =  85)
      integer INO
      parameter (INO       =  86)
      integer INO2
      parameter (INO2      =  87)
      integer INO3
      parameter (INO3      =  88)
      integer IO3
      parameter (IO3       =  91)
      integer IOCLO
      parameter (IOCLO     =  92)
      integer IOH
      parameter (IOH       =  94)
      integer IPMN
      parameter (IPMN      =  96)
      integer IC3H6
      parameter (IC3H6     = 101)
      integer IR4N2
      parameter (IR4N2     = 104)
      integer IRCHO
      parameter (IRCHO     = 109)
      integer ISO2
      parameter (ISO2      = 120)
      integer IACET
      parameter (IACET     = 122)
      integer INITROGEN
      parameter (INITROGEN = 123)
      integer IOXYGEN
      parameter (IOXYGEN   = 124)
      integer IHNO3COND
      parameter (IHNO3COND = 126)
      integer IMGAS
      parameter (IMGAS     = 125)

!.... Additional model-required parameters

      integer IN2
      parameter (IN2      = 123)

      integer IO2
      parameter (IO2      = 124)

      integer IDEHYD
      parameter (IDEHYD   =   0)

      integer ISYNOZ
      parameter (ISYNOZ   =   0)

      integer JBRONO21
      parameter (JBRONO21 =   0)

      integer KHNO31
      parameter (KHNO31   =   0)


!                                  --^--

