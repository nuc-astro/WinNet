!***********************************************************************
!
!     units_module
!
!***********************************************************************

      module units_module

      type units_in_cgs
        real :: e,pi,Na                !numbers
        real :: fm,km,MeV,s         !units
        real :: c,ec,ga,G,GMs_c2km,gv,kB,mb,Ms,Rs !cgs constants
        real :: h,GF,me,Q           !non-cgs constants
      end type units_in_cgs

      type(units_in_cgs), parameter :: units = units_in_cgs(
     &  2.718281828459045,          !e
     &  3.1415926535898e+00,        !pi
     &  6.0221367e23,               !Na Avogadro number
     &  1.e-13,                     !fm
     &  1.e+05,                     !km
     &  1.602192e-06,               !MeV
     &  1.,                         !s
     &  2.997924562e+10,            !c
     &  4.80324214e-10,             !ec
     &  1.23,                       !ga
     &  6.672387286039493e-08,      !G
     &  1.47664,                    !GMs_c2km
     &  1.,                         !gv
     &  1.38066244e-16,             !kB
     &  1.674e-24,                  !mb
     &  1.989e+33,                  !Ms
     &  6.9599e+10,                 !Rs
     &  4.1356943e-21,              !h [MeV*s]
     &  8.957e-44,                  !GF [MeV*cm^3]
     &  0.511,                      !me [MeV]
     &  1.2935)                     !Q [MeV]

      end module units_module

!***********************************************************************
