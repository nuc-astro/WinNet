!> @file format_class.f90
!!
!! The error file code for this file is ***W21***.
!! @brief Module \ref format_class with format statements for data files
!!

!> Define custom format statements used to read in major data files
!!
!! @author Christian Winteler
!! @date   22.01.09
!!
!! \b Editors:
!!   - 11.01.14: Oleg Korobkin
!!   - 28.07.22: MR - Changed format of reaclib to include chapter 9-11
!! .
module format_class
  implicit none

  character(60),dimension(100) :: my_format

contains

  subroutine load_format()

    my_format(1) = "(I2,3X,6A5,8X,A4,a1,a1,3X,1pE12.5)"
    my_format(2) = "(4E13.6)"
    my_format(3) = "(I1,4X,6I5,8X,A4,L1,L1,3X,1PE12.5,/,4E13.6,/, 3E13.6)"
    my_format(4) = "(A5,F12.3,2I4,F6.1,F10.3)"
    my_format(5) = "(8F9.2)"

  end subroutine load_format

end module format_class
