!> @file pardiso_shim.f90
!!
!! @brief Drop-in replacement for Intel MKL PARDISO using SuiteSparse UMFPACK
!!
!! Replaces pardisoinit/pardiso/mkl_dcscmv with UMFPACK-backed implementations.
!! For use on non-Intel platforms (e.g., ARM Mac with gfortran + SuiteSparse).
!!
!! When compiled without UMFPACK (-DUSE_DENSE_SOLVER), falls back to dense
!! LAPACK dgesv (original shim behavior, slow for large networks).
!!
!! @author Maho (OpenClaw/Future Gadget Lab)
!! @date   2026-02-20
!!

!> Initialize PARDISO-compatible state (no-op for UMFPACK backend)
subroutine pardisoinit(pt, mtype, iparm)
  implicit none
  integer*8, intent(out) :: pt(64)
  integer, intent(in)    :: mtype
  integer, intent(out)   :: iparm(64)
  pt = 0
  iparm = 0
  iparm(1) = 1   ! non-default values
  iparm(3) = 1   ! number of threads
  iparm(8) = 10  ! iterative refinement steps
end subroutine pardisoinit

!> Solve sparse linear system using UMFPACK (PARDISO-compatible interface)
!!
!! Replaces Intel MKL PARDISO call. Input is in CSC format:
!!   ia(1:n+1) = column pointers (1-based Fortran indexing)
!!   ja(1:nnz) = row indices (1-based Fortran indexing)
!!   a(1:nnz)  = values
!!
!! UMFPACK uses 0-based C indexing, so we convert on the fly.
!!
subroutine pardiso(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
                   perm, nrhs, iparm, msglvl, b, x, error)
#ifdef USE_DENSE_SOLVER
  ! Dense LAPACK fallback
  implicit none
  integer*8, intent(inout) :: pt(64)
  integer, intent(in)      :: maxfct, mnum, mtype, phase, n, nrhs, msglvl
  real(8), intent(in)      :: a(*)
  integer, intent(in)      :: ia(*), ja(*)
  integer, intent(inout)   :: perm(*)
  integer, intent(inout)   :: iparm(64)
  real(8), intent(inout)   :: b(n)
  real(8), intent(out)     :: x(n)
  integer, intent(out)     :: error

  real(8), allocatable :: A_dense(:,:)
  integer, allocatable :: ipiv(:)
  integer              :: i, j, col_start, col_end, row, info

  error = 0
  if (phase == -1) return

  allocate(A_dense(n, n), ipiv(n))
  A_dense = 0.0d0

  do j = 1, n
    col_start = ia(j)
    if (col_start .lt. 0) cycle  ! empty column (filtered by netsolve)
    col_end   = ia(j+1) - 1
    if (col_end .lt. col_start) cycle  ! also handle via next valid pointer
    do i = col_start, col_end
      row = ja(i)
      ! Store transposed: CSC input but PARDISO expects CSR,
      ! so we transpose to match PARDISO's interpretation
      A_dense(j, row) = a(i)
    end do
  end do

  x = b
  call dgesv(n, nrhs, A_dense, n, ipiv, x, n, info)
  if (info /= 0) error = info
  b = x

  deallocate(A_dense, ipiv)

#else
  ! UMFPACK sparse solver (default)
  use umfpack_interface
  implicit none
  integer*8, intent(inout) :: pt(64)
  integer, intent(in)      :: maxfct, mnum, mtype, phase, n, nrhs, msglvl
  real(8), intent(in)      :: a(*)
  integer, intent(in)      :: ia(*), ja(*)
  integer, intent(inout)   :: perm(*)
  integer, intent(inout)   :: iparm(64)
  real(8), intent(inout)   :: b(n)
  real(8), intent(out)     :: x(n)
  integer, intent(out)     :: error

  ! UMFPACK variables
  real(c_double)  :: Control(UMFPACK_CONTROL)
  real(c_double)  :: Info(UMFPACK_INFO)
  type(c_ptr)     :: Symbolic, Numeric
  integer(c_int)  :: status
  integer(c_int)  :: n_c

  ! 0-based index arrays for UMFPACK
  integer(c_int), allocatable :: Ap(:), Ai(:)
  real(c_double), allocatable :: Ax(:)
  real(c_double), allocatable :: X_umf(:), B_umf(:)
  integer :: i, j, cnt, nnz

  error = 0

  ! Phase -1: release memory (no-op)
  if (phase == -1) return

  n_c = int(n, c_int)

  ! Count non-zero entries
  nnz = ia(n+1) - 1

  ! Allocate and convert to 0-based indexing for UMFPACK
  allocate(Ap(n+1), Ai(nnz), Ax(nnz))
  allocate(X_umf(n), B_umf(n))

  ! Convert column pointers: Fortran 1-based → C 0-based
  ! Handle empty columns: netsolve() sets ia(j)=-1 for columns with no
  ! non-zero entries after filtering. UMFPACK requires monotonically
  ! non-decreasing column pointers, so fill empty columns with the
  ! next valid pointer value.
  Ap(n+1) = int(ia(n+1) - 1, c_int)
  do i = n, 1, -1
    if (ia(i) .lt. 0) then
      Ap(i) = Ap(i+1)  ! empty column: same pointer as next column
    else
      Ap(i) = int(ia(i) - 1, c_int)
    end if
  end do

  ! Convert row indices: Fortran 1-based → C 0-based
  do i = 1, nnz
    Ai(i) = int(ja(i) - 1, c_int)
  end do

  ! Copy values
  do i = 1, nnz
    Ax(i) = a(i)
  end do

  ! Copy RHS
  B_umf = b(1:n)

  ! Set defaults and enable iterative refinement (matching PARDISO iparm(8)=10)
  call umfpack_di_defaults(Control)
  Control(8) = 10.0d0  ! UMFPACK_IRSTEP = max 10 iterative refinement steps

  ! Symbolic factorization
  Symbolic = c_null_ptr
  status = umfpack_di_symbolic(n_c, n_c, Ap, Ai, Ax, Symbolic, Control, Info)
  if (status /= UMFPACK_OK) then
    error = int(status)
    deallocate(Ap, Ai, Ax, X_umf, B_umf)
    return
  end if

  ! Numeric factorization
  Numeric = c_null_ptr
  status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, Numeric, Control, Info)
  call umfpack_di_free_symbolic(Symbolic)
  if (status /= UMFPACK_OK) then
    error = int(status)
    deallocate(Ap, Ai, Ax, X_umf, B_umf)
    return
  end if

  ! Solve A'x = b (transpose) because netsolve() stores the Jacobian in CSC
  ! format but PARDISO interprets ia/ja as CSR (row pointers/column indices).
  ! This means PARDISO effectively solves the transposed system.
  ! UMFPACK correctly interprets our input as CSC, so we must explicitly
  ! request the transpose solve to match PARDISO's behavior.
  status = umfpack_di_solve(UMFPACK_At, Ap, Ai, Ax, X_umf, B_umf, &
                             Numeric, Control, Info)
  call umfpack_di_free_numeric(Numeric)
  if (status /= UMFPACK_OK) then
    error = int(status)
    deallocate(Ap, Ai, Ax, X_umf, B_umf)
    return
  end if

  ! Copy solution to output
  x(1:n) = X_umf(1:n)
  b(1:n) = X_umf(1:n)

  deallocate(Ap, Ai, Ax, X_umf, B_umf)
#endif

end subroutine pardiso


!> Sparse matrix-vector multiply: y = alpha * op(A) * x + beta * y
!! Replaces mkl_dcscmv for CSC format
!! transa='T' means A^T * x, 'N' means A * x
subroutine mkl_dcscmv(transa, m, k, alpha, matdescra, val, rowind, &
                       pntrb, pntre, x, beta, y)
  implicit none
  character, intent(in)    :: transa
  integer, intent(in)      :: m, k
  real(8), intent(in)      :: alpha, beta
  character(*), intent(in) :: matdescra
  real(8), intent(in)      :: val(*)
  integer, intent(in)      :: rowind(*), pntrb(*), pntre(*)
  real(8), intent(in)      :: x(*)
  real(8), intent(inout)   :: y(*)

  integer :: i, j, idx
  real(8) :: tmp

  ! Scale y by beta
  if (beta == 0.0d0) then
    do i = 1, m
      y(i) = 0.0d0
    end do
  else if (beta /= 1.0d0) then
    do i = 1, m
      y(i) = beta * y(i)
    end do
  end if

  ! CSC format: column j has values from pntrb(j) to pntre(j)-1
  if (transa == 'N' .or. transa == 'n') then
    do j = 1, k
      tmp = alpha * x(j)
      do idx = pntrb(j), pntre(j) - 1
        i = rowind(idx)
        y(i) = y(i) + tmp * val(idx)
      end do
    end do
  else
    do j = 1, k
      tmp = 0.0d0
      do idx = pntrb(j), pntre(j) - 1
        i = rowind(idx)
        tmp = tmp + val(idx) * x(i)
      end do
      y(j) = y(j) + alpha * tmp
    end do
  end if

end subroutine mkl_dcscmv
