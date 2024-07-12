!> @file winnse_module.f90
!!
!! The error file code for this file is ***W45***.
!!
!! @brief Module \ref winnse_module with subroutines needed for calculating
!!        nuclear statistical equilibrium
!!
!! @author Christian Winteler
!! @date   25.03.10
!!
!! \b Edited:
!!    - OK:  07.11.16
!!    - MR:  11.01.21  - implemented error_msg_class
!!    - MR:  15.01.21  - Made it a module
!!    - MR:  14.06.23  - Inplemented Powell's hybrid method

!> @brief Module to calculate NSE
!!
!! The winnse_module calculates NSE composition (with or without screening corrections).
!! For this it uses various inputs as, e.g., the partition functions,
!! binding energies, and the spin of the ground state of each nucleus.
!! The necessary data is stored in \ref global_class::isotope.
!!
!! The NSE composition favors different nuclei, dependent on the conditions.
!! A rough overview is given by:
!! @image html nse_com.png "Dominant nuclear species when assuming NSE for different conditions." width=600
!!
!! @see [Hix & Thielemann 1999](https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..321H/abstract),
!!      [Winteler 2013](https://edoc.unibas.ch/29895/),
!!      [Lippuner & Roberts 2017](https://ui.adsabs.harvard.edu/abs/2017ApJS..233...18L/abstract),
!!      [Smith et al. 2023](https://ui.adsabs.harvard.edu/abs/2022arXiv221009965S/abstract)
!!
#include "macros.h"
module winnse_module
   implicit none

   real(r_kind),dimension(:),allocatable,public   :: ynse   !< NSE abundances
   real(r_kind),dimension(:),allocatable,public   :: pf     !< Partition function for a given temperature
   real(r_kind),dimension(:),allocatable,private  :: gg     !< \f$ G(Z,A)=(2J_0+1) \f$ where \f$ J_0 \f$ is the spin of the ground state
   real(r_kind),dimension(:),allocatable,private  :: be     !< Binding energy, calculated using mass excesses
   real(r_kind),dimension(:),allocatable,private  :: scrn   !< Screening correction (details see \ref nse_screen)
   real(r_kind),dimension(:),allocatable,public   :: cnse   !< NSE coefficients (details see \ref cnsecalc)
   real(r_kind),public                            :: ye_ext !< Electron fraction to pass to external function

   integer,dimension(:),allocatable,public  :: aa   !< Mass number
   integer,dimension(:),allocatable,public  :: zz   !< Proton number
   integer,dimension(:),allocatable,public  :: nn   !< Neutron number

   integer,private :: zmax !< Highest proton number occuring in network

   !
   ! Public and private fields and methods of the module.
   !
   public:: &
        winnse_guess, winnse_descend, nse_init
   private:: &
        winnse_calc, nse_screen, wincnse_calc, winnse_calc_hybrid_powell, winnse_calc_NR

contains

!>
!! @brief Allocates and initialises various arrays needed for the nse calculation
!!
!! This subroutine also writes the content of \ref global_class::isotope
!! into shorter variables, such as \ref aa, \ref nn, and \ref gg. In
!! addition it calculates the binding energy
!! \f[ Be(i) = Z(i) \cdot (Z_{me} + e_{m}) + N(i) \cdot N_{me} - M_{exc}(i) \f]
!! with the mass excess of protons (Z), neutrons (N), and the isotope (i), \f$ Z_{me} \f$,
!! \f$ N_{me} \f$, and \f$ M_{exc} \f$. The mass of the electrons (e) is given by \f$ e_{m} \f$.
!! Furthermore it calculates \f$ 1 + 2 \cdot J \f$ with J being the spin
!! of the ground state.
!!
!! \b Edited: 12.01.14
subroutine nse_init()
   use global_class, only: isotope, net_size

   implicit none
   !
   integer                   :: i
   real(r_kind),parameter    :: nem = 8.071323d0 !neutron excess mass [MeV]
   real(r_kind),parameter    :: pem = 7.288969d0 !hydrogen excess mass (including the electron) [MeV]

!-----allocate arrays of dimension net_size
   allocate(ynse(net_size),pf(0:net_size),gg(net_size),be(net_size),cnse(net_size))
   allocate(aa(net_size),zz(net_size),nn(net_size))

!-----initialise arrays
   do i=1,net_size
      aa(i) = isotope(i)%mass
      nn(i) = isotope(i)%n_nr
      zz(i) = isotope(i)%p_nr
      gg(i)  = 1.d0 + 2.d0*isotope(i)%spin
      be(i) = zz(i)*pem + nn(i)*nem - isotope(i)%mass_exc
   end do
   be(1) = 0.d0
   be(2) = 0.d0

!-----find maximum value in zz
   zmax = maxval(zz)
   allocate(scrn(zmax))

   return

end subroutine nse_init


!>
!! @brief Calculate NSE composition with an initial guess.
!!
!! The subroutine tries to calculate nse with an initial guess and if it does
!! not converge, it descends from a high temperature. This routine is used
!! if the parameter \ref parameter_class::nse_calc_every is not equal to zero.
!!
!! @author Moritz Reichert
subroutine winnse_guess(t9, rho, ye, yn_guess, yp_guess, ysol)
   use parameter_class, only: nse_descend_t9start, nse_max_it
   use global_class,    only: net_size
   implicit none

   real(r_kind),intent(in)     :: t9          !< desired temperature in GK
   real(r_kind),intent(in)     :: rho         !< density in g/cm^3
   real(r_kind),intent(in)     :: ye          !< electron fraction
   real(r_kind),intent(in)     :: yn_guess    !< neutron abundance
   real(r_kind),intent(in)     :: yp_guess    !< proton abundance
   real(r_kind),dimension(net_size),intent(inout) :: ysol !< [out] nse abundances

   real(r_kind)                :: yn,yp
   integer                     :: kit,i
!    integer,parameter           :: kmax=25

   INFO_ENTRY("winnse_guess")


   ynse = ysol
   yn = yn_guess
   yp = yp_guess


   do
      call winnse_calc(t9, rho, ye, yn, yp, nse_max_it, kit)

      if ((kit .ge. 1) .and. (kit .le. nse_max_it))  then
         !-convergence successful
         ! Terminate
         exit
      else
         ! If no convergence is achieved, try to descend
         call winnse_descend(nse_descend_t9start, t9, rho, ye, yn, yp, ysol)
         ! If the process is still alive, terminate.
         exit
      end if

   end do

   do i=1,net_size
      if(ynse(i).lt.1.d-25) ynse(i) = 0.d0
   end do

   ysol = ynse

   INFO_EXIT("winnse_guess")
end subroutine winnse_guess


!>
!! @brief This routine descends from an initially high temperature, at which
!! the NSE abudances can be accurately predicted, to the desired
!! temperature
!!
!! At high temperatures, the solution is more stable and the solution of NSE
!! at a high temperature can be used as an initial guess of a lower temperature.
!! This increases the stability of the calculation and it can iterate to the
!! input temperature "t9fnl". The temperature step is adaptive in this process,
!! however, if it becomes too small (\f$ 10^{-20} \f$ GK) an error is raised.
!!
subroutine winnse_descend(t9strt, t9fnl, rho, ye, yni, ypi, ysol)
   use error_msg_class, only: raise_exception,num_to_str,int_to_str
   use parameter_class, only: nse_max_it, nse_delt_t9min
   use screening_module,only: iscreen
   use global_class,    only: net_size
   implicit none

   real(r_kind),intent(in)     :: t9strt !< high initial temperature in GK
   real(r_kind),intent(in)     :: t9fnl  !< desired temperature in GK
   real(r_kind),intent(in)     :: rho    !< density in g/cm^3
   real(r_kind),intent(in)     :: ye     !< electron fraction
   real(r_kind),intent(in)     :: yni    !< neutron abundance
   real(r_kind),intent(in)     :: ypi    !< proton abundance
   real(r_kind),dimension(net_size),intent(inout) :: ysol !< [out] nse abundances

   real(r_kind)                :: t9,delt9,t9hi
   real(r_kind)                :: yn,yn0,yp0,yp
   integer                     :: kit,i

   INFO_ENTRY("winnse_descend")

   t9hi = t9strt
   t9 = t9strt

   if(t9.eq.1.d2) then
      delt9 = sign(10.d0,(t9fnl-t9strt))
   else
      delt9 = 5.d-1*(t9fnl-t9hi)
   end if

   ynse = ysol
   yn0 = yni
   yp0 = ypi
   yn = yni
   yp = ypi


   do

      call winnse_calc(t9, rho, ye, yn, yp, nse_max_it, kit)

      if ((kit .ge. 1) .and. (kit .le. nse_max_it))  then
        !-convergence successful
         if(kit.lt.4) delt9 = 2.d0*delt9
         if(t9.eq.t9fnl) exit

!-----convergence not successful step back and retry with delt9 halved
      else
         if((t9.eq.t9fnl).and.(delt9.eq.0.d0)) then
            call raise_exception("Reached maximum number of iterations ("//&
                                 trim(adjustl(int_to_str(kit)))//") when calculating NSE.",&
                                 "winnse_descend",450003)
         end if
         t9 = t9-delt9
         delt9 = 5.d-1*delt9
      end if

!-----make sure delt9 is not bigger than t9-t9fnl
      delt9 = sign(min(2.d1,abs(delt9)),delt9)
      delt9 = sign(min(abs(t9-t9fnl),abs(delt9)),delt9)

      if (dabs(delt9).lt.nse_delt_t9min)then
         call raise_exception("Temperature difference was too small when calculating"//NEW_LINE("A")//&
                              "NSE composition (T9= "//trim(adjustl(num_to_str(dabs(t9))))//", delt9="//trim(adjustl(num_to_str(dabs(delt9))))//")."//NEW_LINE("A")//&
                              'Try to set "nse_descend_t9start", "nse_nr_tol", or "nse_max_it" to a higher value.',&
                              "winnse_descend",450004)
      end if
!-----descend in temperature
      t9 = t9+delt9
      if(t9.le.t9fnl) t9 = t9fnl !TODO: this can be wrong

   end do

   do i=1,net_size
      if(ynse(i).lt.1.d-25) ynse(i) = 0.d0
   end do

   ysol = ynse

   INFO_EXIT("winnse_descend")

end subroutine winnse_descend



!> Solves the nse-equations.
!!
!! This subroutine solves the NSE equations. It hereby
!! choses between two different solvers as defined by the
!! parameter \ref parameter_class::nse_solver that is defined as
!! in the following table:
!! | nse_solver | Solver                 |
!! |:----------:|:----------------------:|
!! | 0          | Newton-Raphson         |
!! | 1          | Powell's hybrid method |
!!.
!! @author M. Reichert
!! @date   14.06.23
subroutine winnse_calc(t9, rho, ye, yn, yp, imax, kit)
    use screening_module,only: iscreen
    use error_msg_class, only: raise_exception,int_to_str
    use parameter_class, only: nse_solver
    implicit none
    real(r_kind),intent(in)    :: t9   !< temperature in GK
    real(r_kind),intent(in)    :: rho  !< density in g/cm^3
    real(r_kind),intent(in)    :: ye   !< electron fraction
    real(r_kind),intent(inout) :: yn   !< neutron abundance
    real(r_kind),intent(inout) :: yp   !< proton abundance
    integer,intent(in)         :: imax !< max. number of iterations
    integer,intent(out)        :: kit  !< actual number of iterations

    INFO_ENTRY("winnse_calc")

    ! Calculate screening coefficients needed for wincnse_calc
    if (iscreen) call nse_screen(t9, rho, ye)

    ! Calculate nse coefficients cnse
    call wincnse_calc(t9, rho)

    ! Choose the NSE solver
    if (nse_solver .eq. 0) then
        call winnse_calc_NR(t9, rho, ye, yn, yp, imax, kit)
    elseif (nse_solver .eq. 1) then
        call winnse_calc_hybrid_powell(t9, rho, ye, yn, yp, imax, kit)
    else
        call raise_exception("Unknown NSE solver ("//int_to_str(nse_solver)//").", &
                             "winnse_calc",450006)
    end if

    INFO_EXIT("winnse_calc")
end subroutine winnse_calc





!> Solves the nse-equations using the Powell hybrid method.
!!
!! This subroutine solves the NSE equations. It's two dimensional as the unknowns
!! are the neutron abundances Yn and proton abundances Yp. The underlying
!! constraints are mass conservation \f$ \sum X_i = \sum Y_i \cdot A_i = 1 \f$
!! and charge conservation \f$ \sum Y_i \cdot Z_i = y_e \f$.
!! The subroutine makes use of the minpack routine fsolve.
!!
!! @author M. Reichert
!! @date   14.06.23
subroutine winnse_calc_hybrid_powell(t9, rho, ye, yn, yp, imax, kit)
    use parameter_class, only: nse_nr_tol
    implicit none
    ! Declare the pass
    real(r_kind),intent(in)    :: t9   !< temperature in GK
    real(r_kind),intent(in)    :: rho  !< density in g/cm^3
    real(r_kind),intent(in)    :: ye   !< electron fraction
    real(r_kind),intent(inout) :: yn   !< neutron abundance
    real(r_kind),intent(inout) :: yp   !< proton abundance
    integer,intent(in)         :: imax !< max. number of iterations
    integer,intent(out)        :: kit  !< actual number of iterations
    !
    integer                    :: info    !< Info variable for fsolve
    real(r_kind)               :: tol     !< Tolerance for convergence (see nse_nr_tol parameter)
    real(r_kind), dimension(2) :: x       !< Array for the solution
    real(r_kind), dimension(2) :: fvec    !< Array for the function values
    external mass_and_charge_conservation !< External function for the system of equations


    ! Set tolerance for convergence
    tol = nse_nr_tol

    ! Set initial guess
    x(1) = yn
    x(2) = yp

    ! Update external Ye
    ye_ext = ye

    ! Use Powell's hybrid method to solve the system of equations
    call fsolve ( mass_and_charge_conservation, 2, x, fvec, tol, info )

    ! Check if the solution is valid
    ! info = 1 means that the solution converged
    if (info .ne. 1) then
        kit = imax+10
    else
        ! Check if the function is really zero
        if (maxval(abs(fvec)) .gt. tol) then
            kit = imax+10
        else
            kit = 1
            yn = x(1)
            yp = x(2)
        end if
    end if

    return

 end subroutine winnse_calc_hybrid_powell



!>
!! @brief Solves the nse-equations using a 2-dimensional newton-raphson scheme
!!
!! This subroutine solves the NSE equations. It's two dimensional as the unknowns
!! are the neutron abundances Yn and proton abundances Yp. The underlying
!! constraints are mass conservation \f$ \sum X_i = \sum Y_i \cdot A_i = 1 \f$
!! and charge conservation \f$ \sum Y_i \cdot Z_i = y_e \f$.
!! Furthermore it uses the so called saha equation:
!! \f[
!! Y(N,Z) = G_{N,Z}(\rho N_A)^{N+Z-1}\frac{\left( N +Z \right)^{3/2}}{2^{N+Z}}\left(\frac{2\pi \hbar}
!!           {m_\mathrm{u} k_{\mathrm{B}}T}\right)^{\frac{3}{2}(N+Z-1)} \mathrm{exp}(B_{N,Z}/
!!           k_B T)Y^{N}_\mathrm{n} Y^Z_\mathrm{p},
!!\f]
!!
!! @warning Here, the abundance of neutrons and protons are used explicitely.
!! A simulation without neutrons and protons specified in the sunet
!! that tries to calculate NSE will result in an error!
!!
!! @see wincnse_calc
!!
!! \b Edited:
!!     - 12.01.14
!!     - 15.01.21
!!     - 23.01.21 MR, Fixed overflows in this routine
!! .
subroutine winnse_calc_NR(t9, rho, ye, yn, yp, imax, kit)
   use error_msg_class, only: raise_exception,num_to_str,int_to_str
   use parameter_class, only: nse_nr_tol
   use screening_module,only: iscreen
   use global_class,    only: net_size
   implicit none

   real(r_kind),intent(in)    :: t9   !< temperature in GK
   real(r_kind),intent(in)    :: rho  !< density in g/cm^3
   real(r_kind),intent(in)    :: ye   !< electron fraction
   real(r_kind),intent(inout) :: yn   !< neutron abundance
   real(r_kind),intent(inout) :: yp   !< proton abundance
   integer,intent(in)         :: imax !< max. number of iterations
   integer,intent(out)        :: kit  !< actual number of iterations
   !
   real(r_kind) :: yn0,ynl,ynr,delyn,testyn
   real(r_kind) :: yp0,ypl,ypr,delyp,testyp
   real(r_kind) :: f,dfdn,dfdp
   real(r_kind) :: g,dgdn,dgdp
   real(r_kind) :: det,detr
   real(r_kind) :: tol,testk
   real(r_kind) :: atst,ztst,ytst,abar,zbar
   real(r_kind) :: infty              !< Infinity to check for overflows
   integer      :: i,j

   infty = HUGE(infty)

!-----set tolerance for convergence
   tol = nse_nr_tol

!-----set maximum number of iteration steps
!      imax = 10

!-----save initial values
   yn0 = yn
   yp0 = yp


!-----calculate screening coefficients needed for wincnse_calc
   if (iscreen) call nse_screen(t9, rho, ye)

!-----calculate nse coefficients cnse
   call wincnse_calc(t9, rho)

!-----iteratively solve the coupled system of equations:
!     f...expression for charge conservation
!     g...expression for mass conservation

   kit = 0

   do i=1,imax
      ynl = dlog(yn)
      ypl = dlog(yp)
!-----calculate new abundances using the nse equation
      !! Take care of possible overflows
      ynse = cnse + nn*ynl + zz*ypl

      !< Ensure that ynse does not get to large
      if (maxval(ynse) .ge. dlog(infty)-5) then
         kit = imax+4
         exit
      end if

      ynse = dexp(ynse)

      ynr = 1.d0/yn
      ypr = 1.d0/yp

!-----calculate electron and mass conservation
      f = sum(zz*ynse) - ye
      g = sum(aa*ynse) -1.d0

!      print *, 'nn, zz, aa = ', nn, zz, aa
!-----calculate the entries in the jacobian, df/dYp,df/dYn,dg/dYp,dg/dYn
      dfdp = sum(ynse*zz*zz)*ypr
      dfdn = sum(ynse*zz*nn)*ynr
      dgdp = sum(ynse*aa*zz)*ypr
      dgdn = sum(ynse*aa*nn)*ynr

!-----calculate determinant of the jacobian

      ! The following large construct is just to ensure that there are
      ! no floating over/underflows
      ! The determinant is just given by
      ! det = dfdp*dgdn - dfdn*dgdp
      if ((dfdp .gt. 0) .and. (dgdn .gt. 0) .and.&
          (dfdn .gt. 0) .and. (dgdp .gt. 0)) then

         ! Check for overflows and correct them
         if ((dlog(dfdp)+dlog(dgdn) .ge. dlog(infty)) .or. &
             (dlog(dfdn)+dlog(dgdp) .ge. dlog(infty))) then
             ! Express it as dfdn*dfdp*(dgdn/dfdn - dgdp/dfdp)
             det = dlog(dfdn)+dlog(dfdp)+dlog(dgdn/dfdn - dgdp/dfdp)
             ! Check if there is still an overflow, if yes, this iteration has failed
             if (det .ge. dlog(infty)) then
                kit = imax+1
                exit
             else
                det = dexp(det)
             end if
         else
            det = dfdp*dgdn - dfdn*dgdp
         end if
      elseif (((dfdp .gt. infty/10.) .and. (dgdn .gt. infty/10.)) .or.&
              ((dfdn .gt. infty/10.) .and. (dgdp .gt. infty/10.))) then
              kit = imax+1
              exit
      else
         det = dfdp*dgdn - dfdn*dgdp
      end if



!-----solve the matrix equation for the changes in yp and yn
      if(det .ne. 0) then
         detr = 1.d0/det
         delyp = (dgdn*f - dfdn*g)*detr
         delyn = (dfdp*g - dgdp*f)*detr
      else
         call raise_exception('Did not converge when trying to calculate NSE. '//NEW_LINE("A")//&
                              'Found zero determinant.',&
                              "winnse_calc",450005)
      end if

!-----update the abundances, check if positive
      testyp = yp - delyp
      testyn = yn - delyn

      if ((testyp.gt.0.d0).and.(testyn.gt.0.d0)) then
         yp = testyp
         yn = testyn
      else
         kit = imax+4
         exit
      end if

!-----check for convergence
      testyp = delyp/yp
      testyn = delyn/yn
      testk = sqrt(testyp**2 + testyn**2 + f**2 + g**2)
      kit = i
      if (testk .lt. tol) exit

   end do

!-----if nse converges, update abundances
   if(kit.lt.imax) then
      ypl = dlog(yp)
      ynl = dlog(yn)
      ynse = dexp(cnse +ypl*zz +ynl*nn)

!-----calculate some composition properties, atst, ztst, ytst are the
!     total mass, electron fraction and abundances, abar and zbar are
!     the average baryon and proton number
      atst = sum(aa*ynse)
      ztst = sum(zz*ynse)
      ytst = sum(ynse)
      abar = atst/(ytst)
      zbar = ztst/(ytst)
   else
!-----set yn and yp to initial values and return
      kit = imax+1
      yn = yn0
      yp = yp0
   end if

   return

end subroutine winnse_calc_NR


!>
!! @brief Calls screening subroutine and adds subsequent proton capture screening
!! corrections in scrn(1:zmax).
!!
!! This subroutine iterates over all possible amount of protons of a nucleus
!! (up to \ref zmax). It uses the subroutine \ref screen to calculate the
!! screening corrections. After it has been called, the array \ref scrn contains
!! summed screening corrections for proton captures from z=1 to z=i for scrn(i)
!! on symmetric nuclei (\f$ A = 2 \cdot Z \f$).
!!
!! \b Edited:
!!    - 12.01.14
!! .
subroutine nse_screen(t9, rho, ye)
   use screening_module, only: screening
   implicit none
   real(r_kind),intent(in)     :: t9     !< temperature in GK
   real(r_kind),intent(in)     :: rho    !< density in g/cm^3
   real(r_kind),intent(in)     :: ye     !< electron fraction
   real(r_kind)                :: hp, ht, hh
   integer                     :: i,z1,z2,a1,a2

   scrn(1) = 0.d0

   z1 = 1
   a1 = 1
   do i=2,zmax
      z2 = i-1
      a2 = 2*z2
      if(z2.eq.1) a2=1
      call screening (t9,rho,z1,z2,a1,a2,ye,0,hh,ht,hp)
!----- add subsequent screening factors
      scrn(i) = scrn(i-1) + hh
   end do
   return
end subroutine nse_screen


!>
!! Calculates the nse coefficients C, as given in
!! [Hix,Thielemann '99](https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..321H/abstract), Eq.25
!!
!! This subroutine calculates:
!! \f[ \frac{G(^A Z)}{2^A} \left( \frac{ \rho N_A }{\theta} \right)^{A-1} A^{\frac{3}{2}}
!! \exp{ \left( \frac{B(^A Z)}{k_B T} \right)} \f]
!! ,where \f$ G \f$ are the partition functions, A the mass number, \f$ \rho \f$
!! the density, \f$ N_A \f$ avogadros number, B the binding energy, \f$ k_B \f$
!! the boltzman constant, and T the temperature. Theta is a helper variable and given
!! by:
!! \f[ \theta = \left( \frac{m_u k_B T}{2\pi \hbar^2} \right) \f]
!! The result is stored in the array \ref cnse.
!!
!! \b Edited:
!!   - 12.01.14
!! .
subroutine wincnse_calc(t9, rho)
   use parameter_class, only: unit
   use screening_module,only: iscreen
   use nucstuff_class,  only: inter_partf
   use global_class,    only: net_size
   implicit none

   real(r_kind),intent(in)     :: t9      !< temperature in GK
   real(r_kind),intent(in)     :: rho     !< density in g/cm^3
   !
   real(r_kind)                :: bkt, bktr
   real(r_kind)                :: ronath  !< rho*Na/theta
   real(r_kind)                :: ln2

!-----interpolate partition functions for current temperature
   call inter_partf (t9, pf)
   pf = dlog(pf)

!-----calculate boltzmann factor 1/kT in MeV^-1
   bkt = 1.d9*t9*unit%k_mev
   bktr= 1.d0/bkt

!-----calculate log[(rho*Na/theta)] see Hix '99 eq(26) for theta
   ronath = dlog(unit%n_a*rho)+1.5d0*dlog((2.d0*unit%pi*unit%hbar_mev**2)/  &
        (bkt*unit%hix))

   ln2 = dlog(2.d0)

!-----calculate the coefficients
   cnse(1) = 0.d0
   cnse(2) = 0.d0

   ! MR: Screening was thrown out previously. Why? Implemented it again...
   if (iscreen) then
      cnse(3:net_size) = dlog(gg(3:net_size)*aa(3:net_size)*dsqrt(dble(aa(3:net_size)))) &
           + pf(3:net_size) - ln2*aa(3:net_size) + (dble(aa(3:net_size))-1.d0)*ronath    &
           + be(3:net_size)*bktr + scrn(zz(3:net_size))
   else
      cnse(3:net_size) = dlog(gg(3:net_size)*aa(3:net_size)*dsqrt(dble(aa(3:net_size)))) &
           + pf(3:net_size) - ln2*aa(3:net_size) + (dble(aa(3:net_size))-1.d0)*ronath    &
           + be(3:net_size)*bktr
   end if



end subroutine wincnse_calc


end module winnse_module



!> Determines the mass and charge conservation for the NSE calculation
!!
!! This subroutine calculates the mass and charge conservation for the NSE
!! calculation. It is used in the Powell hybrid method to solve the system of
!! equations.
!! The first entry of fvec will contain
!! \f[ \sum_i Z_i Y_i - y_e \f]
!! and the second entry will contain
!! \f[ \sum_i A_i Y_i - 1 \f]
!! where \f$ Z_i \f$ is the proton number, \f$ A_i \f$ the mass number, \f$ Y_i \f$
!! the abundance of the isotope i, and \f$ y_e \f$ the electron fraction.
!!
!! @author M. Reichert
!! @date 14.06.23
subroutine mass_and_charge_conservation(n,x,fvec)
    use winnse_module, only: cnse, nn, zz, aa, ynse, ye_ext
    implicit none
    integer                     :: n    !< Number of equations and variables (2)
    real (r_kind), dimension(n) :: fvec !< Vector that contains function evaluations
    real (r_kind), dimension(n) :: x    !< Vector of variables

    real (r_kind) :: ynl !< Logarithm of the neutron abundance
    real (r_kind) :: ypl !< Logarithm of the proton abundance


    ! Define some helpers
    ynl = dlog(x(1))
    ypl = dlog(x(2))
    ! Calculate new abundances using the nse equation
    ! Take care of possible overflows
    ynse = cnse + nn*ynl + zz*ypl

    !< Ensure that ynse does not get to large
    if (maxval(ynse) .ge. dlog(huge(ynl))-5) then
        fvec(1) = huge(ynl)
        fvec(2) = huge(ynl)
    else
        ! Calculate the actual abundances
        ynse = dexp(ynse)

        ! Calculate electron and mass conservation
        fvec(1) = sum(zz*ynse) -ye_ext
        fvec(2) = sum(aa*ynse) -1.d0
    end if


end subroutine mass_and_charge_conservation