!> @file effphase_class.f90
!!
!! The error file code for this file is ***W15***.
!! @brief Module \ref effphase_class for calculating effective phase space integral
!!

!>
!! Calculate the effective phase space integral needed for the log(ft)
!! weak rates.
!!
!! Ref  [1] <a href="http://adsabs.harvard.edu/abs/2001ADNDT..79....1L">
!!           Langanke and Martinez-Pinedo</a>, ADNDT 79, 1-46(2001)
!!
!! @author Christian Winteler
!! @date   28.07.2009
!! @uses   QUADPACK routines
#include "macros.h"
module effphase_class

  use parameter_class

  implicit none

  real(r_kind)      :: me_temp              !< me [MeV]/kT [MeV]
  real(r_kind)      :: chem_me              !< chem [MeV]/me [MeV]
  real(r_kind)      :: qval_me              !< Q [MeV]/me [MeV]

  !
  ! Public fields and methods of the module (nothing is private)
  !
  public:: &
     effphase

contains

  subroutine effphase (temp,chem,qec,phaseint)
    implicit none

    real(r_kind),intent(in)         :: temp     !< temperature in 10^9 K
    real(r_kind),intent(in)         :: chem     !< chemical potential in MeV. including electron rest mass
    real(r_kind),intent(in)         :: qec      !< Q-value for electron capture [MeV]
    real(r_kind),intent(out)        :: phaseint !< phase space integral

    real(r_kind)                    :: xeffc
!---- NAG variables
    real(r_kind)                    :: emin         ! finite limit of integration range
    integer                         :: inf          ! infinite limit (-1 -> (-inf,bound], +1 -> [bound, +inf)
    parameter (inf = 1)                           !                  2 -> (-inf,+inf))
    real(r_kind)                    :: epsabs       ! absolute accuracy required
    real(r_kind)                    :: epsrel       ! rel. accuracy required
    parameter (epsabs = 0.0d0, epsrel = 1.0d-10)
    real(r_kind)                    :: abserr       ! estimate of absolute error |I-Result|
    integer                         :: lw,liw       ! lw: #of subintervals rec.:(800-2000), liw: lw/4
    parameter (lw=1000,liw=lw/4)
    !real(r_kind), dimension(lw)     :: w
    !integer, dimension(liw)         :: iw
    integer                         :: ifail
    integer                         :: ierr, neval
    external xeffc

    me_temp = unit%mass_e / (temp*1.d9*unit%k_MeV)      ! mec^2 in units of kT
    chem_me = chem / unit%mass_e                        ! chemical potential in units of mec^2
    qval_me = qec / unit%mass_e                         ! Q value in units of mec^2

    ifail = -1

    emin = max(-qval_me,1.d0)                           !lower limit of integral ([1] Eq.16)


!---- Integration routines section: at the moment the use of 1) or 5) is recommended!

    !> 1) NAG routine
    !call d01amf(xeffc,emin,inf,epsabs,epsrel,phaseint,abserr,w,lw,iw,liw,ifail)

    !> 2) pure Gauss-Legendre (simple, untested)
    !!call semiint(xeffc,emin,phaseint)

    !> 3) QUADPACK routine (pulic domain)
    call qagi (xeffc, emin, inf, epsabs, epsrel, phaseint, abserr, neval, ierr)
    !if(ierr == 5) then
      !write(*,*) "Warning: Integrand is too small, switched to (large) finite integration range."
      !write(*,*) ierr, "Phaseint: ", phaseint, ", number of evaluations: ", neval, ", abserr: ", abserr
      !call qags (xeffc, emin, 10.d4, epsabs, epsrel, phaseint, abserr, neval, ierr)
      !write(*,*) ierr, "Phaseint: ", phaseint, ", number of evaluations: ", neval, ", abserr: ", abserr
      !write(*,*) qval_me,emin, me_temp, chem_me
    !endif

    return

  end subroutine effphase

end module effphase_class

  !> The integrand in effective psi, eq.[16] of Langanke & Pinedo (2001)
real(r_kind) function xeffc (w)
  use effphase_class
  implicit none

  real(r_kind),intent(in) :: w            !< total (kinetic+rest mass) energy of e-/e+
  real(r_kind)    :: fd_dist              ! Fermi-Dirac Distribution
  real(r_kind)    :: fd_expo              ! exponent in the fd_dist

  fd_expo = me_temp*(w - chem_me)         ! Exponent of [1] Eq(11) where Ee = w*me*c^2
  if (fd_expo < 36.d0) then
     fd_dist =1.d0/(1.d0+dexp(fd_expo))
  else if (fd_expo < 708.d0) then
     fd_dist = dexp(-fd_expo)
  else
     fd_dist = 0.d0
  end if

  xeffc = (w**2)*((qval_me + w)**2)*fd_dist   ! [1] Eq.16 with Snu = 0
  return

end function xeffc
