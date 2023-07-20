!> @file screengva2008.f90
!!
!! @brief Subroutines \ref screen and \ref screening for electron
!!        screening coefficients
!!

!> Module to calculate electron screening.
!!
!! This module calculates the screening function h (=ln(f_screen))
!!
!! @author Urs Frischknecht
!!
!! \b Edited:
!!   - 26.07.22, MR: Made it a module and made hv, htv, and hpv module variables.
!! .
#include "macros.h"
module screening_module


  real(r_kind), dimension(:), allocatable, public :: hv   !< Screening correction
  real(r_kind), dimension(:), allocatable, public :: htv  !< temp. derivative
  real(r_kind), dimension(:), allocatable, public :: hpv  !< density derivative

  logical, public :: iscreen !< Flag whether screening is enabled or not

  !
  ! Public and private fields and methods of the module
  !
  public:: &
     init_screening, screen, screening
  private:: &
     free_energy_kravchuk_yakovlev, screening_kravchuk_yakovlev
contains


  !> Initialize the screening module
  !!
  !! Allocates screening correction arrays.
  !!
  !! @author M. Reichert
  !! @date 26.07.22
  subroutine init_screening(nreac)
    use parameter_class, only: screening_mode
    use error_msg_class, only: raise_exception
    implicit none
    integer,intent(in) :: nreac !< Number of reactions
    integer            :: istat !< Status variable

    if (screening_mode .gt. 0) then
      iscreen = .true.
      ! Allocate the screening corrections
      allocate(hv(nreac),htv(nreac),hpv(nreac),stat=istat)

      if (istat .ne. 0) then
        call raise_exception('Could not allocate screening correction arrays.',&
                             'init_screening',400001)
      end if
    else
      iscreen = .false.
    end if

  end subroutine init_screening


  !> @brief This function calculates the screening coefficients hv.
  !!
  !! It serves as interface between the subroutine \ref screening and
  !! the \ref jacobian_class.
  !!
  !! @author Urs Frischknecht
  !!
  !! \b Edited:
  !!       - 12.01.14
  !!       - 26.07.22, MR: Made hv, htv, and hpv module variables.
  !!       - 28.07.22, MR: Introduced new chapters
  !! .
  subroutine screen(t9,rho,n,ye,mode)
    use global_class, only:rrate,isotope
    implicit none
    real(r_kind), intent(in) :: t9, rho, ye
    integer, intent(in)      :: mode,n
    real(r_kind) :: h1,ht1,hp1,h2,ht2,hp2,h3,ht3,hp3
    integer      :: i,z1,z2,z3,z4,z5,z6,a1,a2,a3,a4,a5,a6

    do i=1,n
       select case (rrate(i)%group)
       case(1:3,11)
          hv(i)=0.d0
          htv(i)=0.d0
          hpv(i)=0.d0
       case(4:7)
          z1=isotope(rrate(i)%parts(1))%p_nr
          z2=isotope(rrate(i)%parts(2))%p_nr
          a1=isotope(rrate(i)%parts(1))%mass
          a2=isotope(rrate(i)%parts(2))%mass
          call screening(t9,rho,z1,z2,a1,a2,ye,mode,h1,ht1,hp1)
          hv(i)=h1
          htv(i)=ht1
          hpv(i)=hp1
       case(8:9)
          z1=isotope(rrate(i)%parts(1))%p_nr
          z2=isotope(rrate(i)%parts(2))%p_nr
          z3=isotope(rrate(i)%parts(3))%p_nr
          z4=z1+z2
          a1=isotope(rrate(i)%parts(1))%mass
          a2=isotope(rrate(i)%parts(2))%mass
          a3=isotope(rrate(i)%parts(3))%mass
          a4=a1+a2
          call screening(t9,rho,z1,z2,a1,a2,ye,mode,h1,ht1,hp1)
          call screening(t9,rho,z3,z4,a3,a4,ye,mode,h2,ht2,hp2)
          hv(i)=h1+h2
          htv(i)=ht1+ht2
          hpv(i)=hp1+hp2
       case(10)
         z1=isotope(rrate(i)%parts(1))%p_nr
         z2=isotope(rrate(i)%parts(2))%p_nr
         z3=isotope(rrate(i)%parts(3))%p_nr
         z4=isotope(rrate(i)%parts(4))%p_nr
         z5=z1+z2
         z6=z3+z4
         a1=isotope(rrate(i)%parts(1))%mass
         a2=isotope(rrate(i)%parts(2))%mass
         a3=isotope(rrate(i)%parts(3))%mass
         a4=isotope(rrate(i)%parts(4))%mass
         a5=a1+a2
         a6=a3+a4
         call screening(t9,rho,z1,z2,a1,a2,ye,mode,h1,ht1,hp1)
         call screening(t9,rho,z3,z4,a3,a4,ye,mode,h2,ht2,hp2)
         call screening(t9,rho,z5,z6,a5,a6,ye,mode,h3,ht3,hp3)
         hv(i)=h1+h2+h3
         htv(i)=ht1+ht2+ht3
         hpv(i)=hp1+hp2+hp3
       end select
    end do

    return
  end subroutine screen



  !> Free energy according to parametrization of Kravchuk and Yakovlev
  !!
  !! This function calculates the free energy of a one-component plasma
  !! according to equation 19 of Kravchuk & Yakovlev (2014).
  !!
  !! @see [Kravchuk & Yakovlev 2014](https://ui.adsabs.harvard.edu/abs/2014PhRvC..89a5802K/abstract)
  !!
  !! @author M. Reichert
  !! @date 03.04.2023
  function free_energy_kravchuk_yakovlev(gamma) result(f_C)
    implicit none
    real(r_kind),intent(in) :: gamma  !< Ion coupling parameter
    real(r_kind)            :: f_C    !< Free energy
    ! Parameters from Kravchuk & Yakovlev (2014)
    real(r_kind),parameter  :: A_1= -0.907
    real(r_kind),parameter  :: A_2=  0.62954
    real(r_kind),parameter  :: A_3=  0.2771
    real(r_kind),parameter  :: B_1=  0.00456
    real(r_kind),parameter  :: B_2=  211.6
    real(r_kind),parameter  :: B_3= -0.0001
    real(r_kind),parameter  :: B_4=  0.00462

    f_C = A_1*( dsqrt(gamma*(A_2 + gamma)) - A_2 * dlog( dsqrt(gamma/A_2) + dsqrt(1d0+gamma/A_2))) &
         +2d0*A_3 * (dsqrt(gamma)-Datan(dsqrt(gamma))) &
         +B_1 * (gamma - B_2 * dlog(1d0+gamma/B_2))   &
         +B_3/2d0 * dlog(1d0 + gamma**2d0 / B_4)

  end function free_energy_kravchuk_yakovlev




  !> Interface for the screening.
  !!
  !! This subroutine is an interface for the different
  !! screening prescriptions.
  !!
  !! @author M. Reichert
  !! @date 03.04.2023
  subroutine screening(t9,rho,z1,z2,a1,a2,ye,mode,h,ht,hp)
    use error_msg_class, only: raise_exception, int_to_str
    use parameter_class, only: screening_mode
    implicit none
    integer,intent(in)       :: mode  !< =0 no derivatives are calculated,
                                      !! =1 ht=dln(f)/d(t9)=dh/dt9,
                                      !!    hp=dln(f)/d(rho)=dh/d(rho)
    real(r_kind),intent(in)  :: t9    !< temperature [GK]
    real(r_kind),intent(in)  :: rho   !< density [g/cm^3]
    integer,intent(in)       :: z1,z2 !< charge numbers of colliding nuclei
    integer,intent(in)       :: a1,a2 !< mass numbers of colliding nuclei
    real(r_kind),intent(in)  :: ye    !< electron fraction
    real(r_kind),intent(out) :: h     !< screening function
    real(r_kind),intent(out) :: ht    !< dln(f)/d(t9)=dh/dt9
    real(r_kind),intent(out) :: hp    !< dln(f)/d(rho)=dh/d(rho)
    !
    integer,parameter :: screen_mode =1

    if (screen_mode .eq. 1) then
        call screening_kravchuk_yakovlev(t9,rho,z1,z2,a1,a2,ye,mode,h,ht,hp)
    else
        call raise_exception('Unknown screening mode, got '//int_to_str(screening_mode)//".",&
                             "screening",400003)
    end if

  end subroutine screening



  !> Screening function according to Kravchuk & Yakovlev (2014)
  !!
  !! This function calculates the screening function according to
  !! equation 62 of Kravchuk & Yakovlev (2014) using the "combined" model.
  !!
  !! @note  This function is in principle only valid for the strong
  !!        screening regime does however not deviate too much for the weak
  !!        screening as well.
  !!
  !! @author M. Reichert
  !! @date 03.04.2023
  subroutine screening_kravchuk_yakovlev(t9,rho,z1,z2,a1,a2,ye,mode,h,ht,hp)
    use error_msg_class, only: raise_exception
    implicit none
    integer,intent(in)       :: mode  !< =0 no derivatives are calculated,
                                      !! =1 ht=dln(f)/d(t9)=dh/dt9,
                                      !!    hp=dln(f)/d(rho)=dh/d(rho)
    real(r_kind),intent(in)  :: t9    !< temperature [GK]
    real(r_kind),intent(in)  :: rho   !< density [g/cm^3]
    integer,intent(in)       :: z1,z2 !< charge numbers of colliding nuclei
    integer,intent(in)       :: a1,a2 !< mass numbers of colliding nuclei
    real(r_kind),intent(in)  :: ye    !< electron fraction
    real(r_kind),intent(out) :: h     !< screening function
    real(r_kind),intent(out) :: ht    !< dln(f)/d(t9)=dh/dt9
    real(r_kind),intent(out) :: hp    !< dln(f)/d(rho)=dh/d(rho)
    !
    real(r_kind) :: zz                !< helper variable for combined z's
    real(r_kind) :: tau               !< helper variable (Eq. 5 of K&Y 2014)
    real(r_kind) :: xi                !< xi = 3 gamma_12 / tau
    real(r_kind) :: gamma_12          !< Combined gamma
    real(r_kind) :: zfrac             !< Fraction of z1 and z2
    real(r_kind) :: gamma_1           !< Gamma of reacting nucleus 1
    real(r_kind) :: gamma_2           !< Gamma of reacting nucleus 2
    real(r_kind) :: gamma_C           !< Gamma of compound nucleus
    real(r_kind) :: b0                !< Salpeter parameter
    real(r_kind) :: b2                !< b2 parameter that gets multiplied to xi**2
    real(r_kind) :: b4                !< b4 parameter that gets multiplied to xi**4
    real(r_kind) :: z1r,z2r           !< z variables as doubles

    h    = 0.d0
    ht   = 0.d0
    hp   = 0.d0

    ! Don't screen if a neutron is involved
    if(z1.eq.0 .or. z2.eq.0) then
       return
    endif

    ! Initialize gammas
    gamma_12= 0.d0
    gamma_1 = 0.d0
    gamma_2 = 0.d0
    gamma_C = 0.d0

    ! Convert the integers to a double
    z1r = dble(z1)
    z2r = dble(z2)
    ! Calculate xi, gamma_12, and tau (Eq. 5 of K&Y 2014)
    zz   = 2.d0*z1r*z2r/(z1r**(1.d0/3.d0)+z2r**(1.d0/3.d0))
    tau  = 3.3722d0*(2.d0*a1*a2/(a1+a2)*(z1r*z2r)**2/t9)**(1.d0/3.d0)
    gamma_12= zz*0.22747d-3*(rho*ye)**(1.d0/3.d0)/t9
    xi=3.d0*gamma_12/tau

    ! Gamma of nucleus 1
    zz=z1r**(5.0/3.0)
    gamma_1= zz*0.22747e-3*(rho*ye)**(1.e0/3.e0)/t9
    ! Gamma of nucleus 2
    zz=z2r**(5.0/3.0)
    gamma_2= zz*0.22747e-3*(rho*ye)**(1.e0/3.e0)/t9
    ! Gamma of compound nucleus
    zz=(z1r+z2r)**(5.0/3.0)
    gamma_C= zz*0.22747e-3*(rho*ye)**(1.e0/3.e0)/t9

    ! Calculate b0 ( Eq. 18 of K&Y 2014)
    b0 = (free_energy_kravchuk_yakovlev(gamma_1) + &
          free_energy_kravchuk_yakovlev(gamma_2) - &
          free_energy_kravchuk_yakovlev(gamma_C))/gamma_12
    ! Calculate b2 ( Eq. 14 of K&Y 2014)
    zfrac = z1r/z2r
    b2 = -1d0/(16d0) * (1d0+zfrac**(1d0/3d0))**3d0 / (1d0 + zfrac)
    ! Calculate b4 ( Eq. 15 of K&Y 2014)
    b4 = zfrac/(64d0) * (1d0+zfrac**(1d0/3d0))**5d0 / ((1d0 + zfrac)**(11d0/3d0))

    ! Calculate the screening correction (Eq 62 of K&Y 2014)
    h = gamma_12 * ( b0 + 5d0/8d0*b2*xi**2d0 + 63d0/128d0*b4*xi**4d0 )

    ! Calculate derivatives if requested
    if (mode.eq.1) then
        call raise_exception('Derivatives not implemented yet.',&
                             'screening_kravchuk_yakovlev', 400004)
    end if

  end subroutine screening_kravchuk_yakovlev


end module screening_module
