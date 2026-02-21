!> @file umfpack_interface.f90
!!
!! @brief Fortran interface to SuiteSparse UMFPACK via iso_c_binding
!!
!! Provides Fortran bindings to the UMFPACK sparse direct solver
!! as a replacement for Intel MKL PARDISO on non-Intel platforms.
!!
!! @author Maho (OpenClaw/Future Gadget Lab)
!! @date   2026-02-20
!!
module umfpack_interface
  use iso_c_binding
  implicit none

  ! UMFPACK constants
  integer(c_int), parameter :: UMFPACK_CONTROL = 20
  integer(c_int), parameter :: UMFPACK_INFO    = 90
  integer(c_int), parameter :: UMFPACK_A       = 0   ! Solve Ax=b
  integer(c_int), parameter :: UMFPACK_OK      = 0

  interface

    !> Set default control parameters
    subroutine umfpack_di_defaults(Control) bind(C, name='umfpack_di_defaults')
      import :: c_double
      real(c_double), intent(out) :: Control(*)
    end subroutine

    !> Symbolic factorization
    function umfpack_di_symbolic(n_row, n_col, Ap, Ai, Ax, Symbolic, &
                                 Control, Info) result(status) &
             bind(C, name='umfpack_di_symbolic')
      import :: c_int, c_double, c_ptr
      integer(c_int), value       :: n_row, n_col
      integer(c_int), intent(in)  :: Ap(*), Ai(*)
      real(c_double), intent(in)  :: Ax(*)
      type(c_ptr),    intent(out) :: Symbolic
      real(c_double), intent(inout) :: Control(*)
      real(c_double), intent(out) :: Info(*)
      integer(c_int)              :: status
    end function

    !> Numeric factorization
    function umfpack_di_numeric(Ap, Ai, Ax, Symbolic, Numeric, &
                                Control, Info) result(status) &
             bind(C, name='umfpack_di_numeric')
      import :: c_int, c_double, c_ptr
      integer(c_int), intent(in)  :: Ap(*), Ai(*)
      real(c_double), intent(in)  :: Ax(*)
      type(c_ptr),    value       :: Symbolic
      type(c_ptr),    intent(out) :: Numeric
      real(c_double), intent(inout) :: Control(*)
      real(c_double), intent(out) :: Info(*)
      integer(c_int)              :: status
    end function

    !> Solve the system
    function umfpack_di_solve(sys, Ap, Ai, Ax, X, B, Numeric, &
                              Control, Info) result(status) &
             bind(C, name='umfpack_di_solve')
      import :: c_int, c_double, c_ptr
      integer(c_int), value       :: sys
      integer(c_int), intent(in)  :: Ap(*), Ai(*)
      real(c_double), intent(in)  :: Ax(*), B(*)
      real(c_double), intent(out) :: X(*)
      type(c_ptr),    value       :: Numeric
      real(c_double), intent(inout) :: Control(*)
      real(c_double), intent(out) :: Info(*)
      integer(c_int)              :: status
    end function

    !> Free symbolic factorization
    subroutine umfpack_di_free_symbolic(Symbolic) &
               bind(C, name='umfpack_di_free_symbolic')
      import :: c_ptr
      type(c_ptr), intent(inout) :: Symbolic
    end subroutine

    !> Free numeric factorization
    subroutine umfpack_di_free_numeric(Numeric) &
               bind(C, name='umfpack_di_free_numeric')
      import :: c_ptr
      type(c_ptr), intent(inout) :: Numeric
    end subroutine

  end interface

end module umfpack_interface
