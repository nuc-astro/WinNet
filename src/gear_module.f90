!> @file gear_module.f90
!!
!! The error file code for this file is ***W22***.
!! @brief Contains the module \ref gear_module
!!

!>
!! @brief \ref gear_module contains adaptive high-order Gear solver
!!
!! The solver is described by, e.g.,
!! [Byrne & Hindmarsh (1975), ACM TOMS 1(1), p.71](https://www.sciencedirect.com/science/article/pii/0021999187900015).
!! The basic idea is to estimate the most efficient way between a larger timestep,
!! but a higher order (and therefore computationally more expensive) or a
!! smaller timestep and a lower order solution.
!! The Gear solver therefore calculates the timestep based on an error estimation
!! originating from the solution using different orders.
!!
!! @author Dirk Martin
!! @date   01.08.17
!!
!! @see
!! [Byrne & Hindmarsh 1975](https://www.sciencedirect.com/science/article/pii/0021999187900015),
!! [Longland et al. 2014](https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..67L/abstract),
!! [Martin, D. 2017](https://tuprints.ulb.tu-darmstadt.de/6301/)
!!
!! \b Edited:
!!   - OK 16.08.17: privatised local fields and methods
!!   - MR 18.04.19: fixed gear for termination criterion 0
!!   - MR 15.01.21: fixed gear for all termination criterions in timestep_module
!!   - MR 22.01.21: Extented comments
!!   - MR 17.07.22: Fixed a bug related to switching the orders. This bug was only visible
!!                  with the compiler flag check_bounds which lead to slightly different results
!!   - MR 17.02.23: Introduced "revert_timestep" to revert the timestep if the solution did not converge
!! .
#include "macros.h"

module gear_module

  use global_class, only: net_size
  use error_msg_class, only: raise_exception
  use parameter_class, only: gear_eps, gear_escale, gear_cFactor, &
                             gear_nr_maxcount, gear_nr_eps
  implicit none

  integer,parameter  :: histsize = 13  !< length of history in tau, e, ...
  integer,parameter  :: qmax     = 5   !< maximum order of q

  real(r_kind),dimension(6,6) :: Amat = &
     reshape( (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, &
                 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, &
                 0.0d0, 0.0d0, 1.0d0, 3.0d0, 6.0d0,10.0d0, &
                 0.0d0, 0.0d0, 0.0d0, 1.0d0, 4.0d0,10.0d0, &
                 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 5.0d0, &
                 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/6, 6/) ) !< Helper Matrix
     !< This Matrix is defined by:
     !< \f[
     !< A_{ij}=\begin{cases}0  & \text{if } i<j\\\begin{pmatrix}i\\j\end{pmatrix}
     !< = \frac{i!}{j!(i-j)!} & \text{if } i\ge j\end{cases} \qquad \text{with } i,j \in [0,1,...,q]
     !< \f]
     !< The Taylor-series of \f$\vec{y}_n\f$ truncated at order q, can be calculated as
     !< the product of the Nordsieck vector (\ref z).

  integer :: ifac
  integer, dimension(qmax+1) :: helperArr = (/ (ifac, ifac=1,qmax+1) /)
  integer, dimension(qmax+1) :: factorial

  !!! Newton-Raphson options
  integer                :: gear_nr_count     = 0       !< counter for NR steps

  !!! iteration variables
  integer      :: q, n, nq
  real(r_kind) :: ti, h

  !!! initial iteration values
  integer      :: q_ini   = 1      !< start with order 1
  integer      :: n_ini   = 1      !< iteration step
  integer      :: nq_ini  = 1      !< iteration step at current order

  !!! initialize vectors
  real(r_kind),dimension(histsize)        :: tau     !< hlist inverted
  real(r_kind),dimension(:,:),allocatable :: z       !< Nordsieck vector
      !< It is defined as
      !< \f[
      !< \vec{z}_n = \left[\vec{y}_n,h\dot{\vec{y}}_n,\frac{h^2 \ddot{\vec{y}}_n}{2!}
      !<             ,...,\frac{h^q \vec{y}_n^{(q)}}{q!} \right]
      !< \f]

  real(r_kind),dimension(:,:),allocatable :: znew !< Predictor step
  real(r_kind),dimension(histsize)        :: xi
  real(r_kind),dimension(:),allocatable   :: en
  real(r_kind),dimension(:,:),allocatable :: e   !list of en
  real(r_kind),dimension(histsize)        :: l, el
  real(r_kind),dimension(:),allocatable   :: errq,errqm1,errqp1

  integer :: alloc_stat
  real(r_kind) :: hq, hqp1, hqm1

  !
  ! Public and private fields and methods of the module
  !
  public:: &
     gear_escale,gear_cFactor,gear_nr_maxcount,gear_nr_eps, &
     gear_nr_count
  public:: &
     init_gear_solver,  init_dYdt,        get_time,   get_timestep, &
     get_solution,    get_predictor_Y,    get_l1,     set_timestep, &
     get_solution_at, get_predictor_dYdt, shift_tau,  set_xi,set_l, &
     set_new_result,  nordsieck_product,  prepare_next_step,        &
     revert_timestep


  private:: &
     histsize,qmax,Amat,ifac,helperArr,factorial,gear_eps,q,n,nq,ti,h, &
     q_ini,n_ini,nq_ini,tau,z,znew,xi,en,e,l,el,errq,errqm1,errqp1,    &
     alloc_stat,hq,hqp1,hqm1
  private:: &
     geterrors,errorqm1,errorqp1,qfun,shiftorder,get_cFactor

contains

!>
!! Initialize iteration variables
!!
!! Here, all necessary allocations are done. And
!! the initial abundance is written into the abundance
!! history array (\ref z(1,:)).
!!
!! @author D. Martin
!! @date   01.08.17
subroutine init_gear_solver(Y,t_init,h_init)
  implicit none

  real(r_kind), dimension(net_size), intent(in) :: Y
  real(r_kind), intent(in) :: t_init
  real(r_kind), intent(in) :: h_init

  ! initialize factorial function
  do ifac=1,qmax+1
     factorial(ifac) = product(helperArr(1:ifac))
  enddo

  ! allocate arrays involving net_size
  if (.not. allocated(z)) then
     allocate(z(histsize,net_size),stat=alloc_stat)
     if ( alloc_stat /= 0) call raise_exception('allocating z failed','init_gear_solver',220001)
     allocate(znew(histsize,net_size),stat=alloc_stat)
     if ( alloc_stat /= 0) call raise_exception('allocating znew failed','init_gear_solver',220001)
     allocate(en(net_size),stat=alloc_stat)
     if ( alloc_stat /= 0) call raise_exception('allocating en failed','init_gear_solver',220001)
     allocate(e(histsize,net_size),stat=alloc_stat)
     if ( alloc_stat /= 0) call raise_exception('allocating e failed','init_gear_solver',220001)
     allocate(errq(net_size),stat=alloc_stat)
     if ( alloc_stat /= 0) call raise_exception('allocating errq failed','init_gear_solver',220001)
     allocate(errqm1(net_size),stat=alloc_stat)
     if ( alloc_stat /= 0) call raise_exception('allocating errqm1 failed','init_gear_solver',220001)
     allocate(errqp1(net_size),stat=alloc_stat)
     if ( alloc_stat /= 0) call raise_exception('allocating errqp1 failed','init_gear_solver',220001)
  end if

  q  = q_ini   ! start with order 1
  ti = t_init  ! initial time
  h  = h_init  ! initial step size
  n  = n_ini   ! iteration step
  nq = nq_ini  ! iteration step at current order

  z(1:histsize,1:net_size) = 0.0
  z(1,1:net_size) = Y(1:net_size)


end subroutine init_gear_solver


!>
!! Initialize dYdt component of the Nordsieck vector.
!!
!! Stores h * dYdt in the Nordsieck vector \ref z(2,:).
!!
!! @author D. Martin
!! @date   01.08.17
subroutine init_dYdt(dYdt)
  implicit none

  real(r_kind), dimension(net_size), intent(in) :: dYdt

  z(2,1:net_size) = h*dYdt

end subroutine init_dYdt


!>
!! The current time
!!
!! @author D. Martin
!! @date   01.08.17
function get_time()
   implicit none

   real(r_kind) :: get_time

   get_time = ti
end function get_time


!>
!! The current timestep
!!
!! @author D. Martin
!! @date   01.08.17
function get_timestep()
   implicit none

   real(r_kind) :: get_timestep

   get_timestep = h
end function get_timestep


!>
!! Set the current timestep
!!
!! @warning When changing the timestep
!! of Gear, the Nordsieck vector (\ref z) has
!! to be rescaled because it depends on the timestep.
!! The timestep of gear should therefore only changed
!! by this routine!
!!
!! @see timestep_module::restrict_timestep,
!!      network_init_module::switch_evolution_mode.
!!
!! @author D. Martin
!! @date   01.08.17
subroutine set_timestep(timestep)
   implicit none

   real(r_kind), intent(in) :: timestep

   integer :: j


   ! rescale Nordsieck vector
   do j=2,q+1
      z(j,1:net_size) = z(j,1:net_size)*((timestep/h)**(j-1.d0))
   enddo

   h = timestep

   ! Check the order afterwards again
   nq = q
end subroutine set_timestep


!> Revert a timestep
!!
!! This routine reverts the timestep of Gear
!! to the previous timestep.
!! It is used in case the result did not converge
!!
!! @see timestep_module::advance_gear
!!
!! @author M. Reichert
subroutine revert_timestep(timestep)
    implicit none
    real(r_kind), intent(in) :: timestep

    ti = ti - h

    call set_timestep(timestep)
    tau = eoshift(tau, 1)

end subroutine revert_timestep


!>
!! The current l_1 (in Fortran: l(2))
!!
!! @see set_l
!!
!! @author D. Martin
!! @date   01.08.17
function get_l1()
   implicit none

   real(r_kind) :: get_l1

   get_l1 = l(2)
end function get_l1


!>
!! Gives a copy of the solution at the current time
!!
!! The solution abundance is stored at \ref z(1,:).
!!
!! @author D. Martin
!! @date   01.08.17
function get_solution()
   implicit none

   real(r_kind), dimension(net_size) :: get_solution

   get_solution = z(1,1:net_size)
end function get_solution


!>
!! Determines the predicted y_next
!!
!! The predicted abundance is stored in \ref z_new(1,:)
!!
!!
!! @author D. Martin
!! @date   01.08.17
function get_predictor_Y()
   implicit none

   real(r_kind), dimension(net_size) :: get_predictor_Y

   get_predictor_Y = znew(1,1:net_size)

end function get_predictor_Y


!>
!! Determines the predicted dydt_next
!!
!! This returns the predicted dydt from the predictor step
!! from \ref z_new.
!!
!! @author D. Martin
!! @date   01.08.17
function get_predictor_dYdt()
   implicit none

   real(r_kind), dimension(net_size) :: get_predictor_dYdt

   get_predictor_dYdt = znew(2,1:net_size)
end function get_predictor_dYdt


!>
!! Determines the solution at a given point in time for [ti-h, ti]
!!
!! @author D. Martin
!! @date   01.08.17
subroutine get_solution_at(time_inter, Y_inter)
   implicit none

   real(r_kind), intent(in) :: time_inter
   real(r_kind), dimension(net_size), intent(out) :: Y_inter

   integer :: j

   Y_inter(1:net_size) = 0.0d0
   if ((time_inter.gt.ti-tau(1)) .and. (time_inter.le.ti)) then
      do j=1,q+1
         Y_inter = Y_inter + (((time_inter - ti)/tau(1))**(j-1.0d0))*z(j,1:net_size)
      enddo
   endif

   return

end subroutine get_solution_at


!>
!! Shifts tau vector and controls the current time
!!
!! @author D. Martin
!! @date   01.08.17
subroutine shift_tau
  implicit none

  tau = eoshift(tau, -1)
  tau(1) = h

  ti = ti + h

end subroutine shift_tau


!>
!! Nordsieck vector-matrix product
!!
!! This calculates:
!! \f[
!! \vec{z}_{n+1}^{(0)} = \mathbf{A} \cdot \vec{z}_n
!! \f]
!! and stores the result in \ref z_new.
!!
!! @author D. Martin
!! @date   01.08.17
subroutine nordsieck_product
  implicit none

  integer :: qp1

  integer :: i,j
  znew = 0.0d0

  qp1 = q+1

  do i=1,qp1
    do j=1,qp1
      znew(i,1:net_size) = znew(i,1:net_size)+Amat(j,i)*z(j,1:net_size)
    enddo
  enddo
end subroutine nordsieck_product


!>
!! Function to calculate xi (needed for getting l in the corrector step)
!!
!! xi is defined as:
!! \f[
!! \xi _i = \frac{t_{n+1} - t_{n+1-i}}{h}
!! \f]
!!
!! @author D. Martin
!! @date   01.08.17
subroutine set_xi
  implicit none

  integer      :: j

  xi(1:histsize) = 0.0d0
  do j=1,q
    xi(j) = xi(j)+sum(tau(1:j))
    xi(j) = xi(j)/h
  enddo

end subroutine set_xi


!>
!! Function to calculate \f$ \ell \f$ (needed for the corrector step)
!! \f[
!! \begin{align}
!! \ell _0 (q) &=1, \nonumber \\
!! \quad \ell_1(q)&=\sum \limits _{i=1}^{q} \left( \xi _i ^{-1} \right), \nonumber\\
!! \quad \ell_j (q) &= \ell_j(q-1)+\ell_{j-1}(q-1)/\xi _q, \nonumber\\
!! \quad \ell _q (q) &= \left( \prod \limits_{i=1}^{q}\xi _i \right) ^{-1}. \nonumber
!! \end{align}
!! \f]
!!
!! @author D. Martin
!! @date   01.08.17
subroutine set_l
  implicit none

  integer      :: j,iback

  l    = 0.0d0
  l(1) = 1.0d0
  l(2) = 1.0d0/xi(1) ! CAREFUL: this is l_1 in the formulas
  do j=2,q
    l(j+1) = l(j)/xi(j)
    do iback=j,2,-1
      l(iback) = l(iback) + l(iback-1)/xi(j)
    enddo !iback
  enddo !j

end subroutine set_l


!>
!! Updates the solution (corrector step)
!!
!! This routine calculates
!! \f[ \vec{z}_{n+1}=\vec{z}_{n+1}^{(0)}+\vec{e}_{n+1}\vec{\ell} \f]
!!
!! @see set_l, timestep_module::advance_gear
!!
!! @author D. Martin
!! @date   01.08.17
subroutine set_new_result(ydiff)
  implicit none

  real(r_kind),dimension(net_size),intent(in) :: ydiff

  integer :: j

  en = ydiff

  ! shift en list
  e = eoshift(e, -1)
  e(1,1:net_size) = en

  ! z -> znew
  z = 0.0d0
  do j=1,q+1
    z(j,1:net_size) = znew(j,1:net_size) + en*l(j)
  enddo

end subroutine set_new_result


!>
!! Stepsize and order control to pepare the next step
!!
!! @author D. Martin
!! @date   01.08.17
subroutine prepare_next_step
  implicit none

  integer :: j

  if(q>1 .and. nq==q) then
    call geterrors


    call shiftorder

    ! rescale Nordsieck vector, if necessary
    do j=2,q+1
      z(j,1:net_size) = z(j,1:net_size)*((h/tau(1))**(j-1.d0))
    enddo

    nq = 0

  else if(q==1) then
    !!! first increase of the order
    q = q+1
    z(q+1,1:net_size) = 0.0d0
    nq = 0
  endif

  nq = nq+1
  n = n+1
end subroutine prepare_next_step


!>
!! Calculate errors (helper functions below)
!!
!! Calculates the errors (E, E(q+1), and E(q-1)) by calling the
!! respective subroutines \ref errorq, \ref errorqm1, \ref errorqp1.
!! Also estimates the resulting timestep from each order.
!!
!! @author D. Martin
!! @date   01.08.17
!!
!! \b Edited:
!!   - M. Reichert, 15.07.22: Set time step of qp1 to zero if already in maximum order.
!! .
subroutine geterrors
  implicit none

  call errorq
  call errorqm1
  call errorqp1 ! only calculate for q<qmax

  errq(1:net_size) = errq(1:net_size)/max(dabs(z(1,1:net_size)),gear_escale)
  errqm1(1:net_size) = errqm1(1:net_size)/max(dabs(z(1,1:net_size)),gear_escale)
  errqp1(1:net_size) = errqp1(1:net_size)/max(dabs(z(1,1:net_size)),gear_escale)

  if(maxval(errq) < 1.d-100) then
    hq = 10.d0*h
  else
    hq = gear_cFactor*h*((gear_eps/maxval(errq))**(1.0d0/(q+1.0d0)))
  endif
  if(maxval(errqm1) < 1.d-100) then
    hqm1 = 10.d0*h
  else
    hqm1 = gear_cFactor*h*((gear_eps/maxval(errqm1))**(1.0d0/q))
  endif

  ! Take care that one does not increase the order for the maximum order
  if (q .ne. qmax) then
    if(maxval(errqp1) < 1.d-100) then
      hqp1 = 10.d0*h
    else
      hqp1 = gear_cFactor*h*((gear_eps/maxval(errqp1))**(1.0d0/(q+2.0d0))) ! see above
    endif
  else
    hqp1 = 0
  end if
!  print *, hq
!  print *, hqm1
!  print *, hqp1

end subroutine geterrors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Error functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! Calculate the errors for the current order
!!
!! The error is defined as:
!! \f[
!! \vec{E}_{n+1}(q) = -\frac{1}{\ell _1}\left[ 1+ \prod _{i=2}^{q}\left(
!! \frac{t_{n+1}-t_{n+1-i}}{t_{n}-t_{n+1-i}} \right) \right] ^{-1} \vec{e}_{n+1}
!! \f]
!!
!! @author D. Martin
!! @date   01.08.17
subroutine errorq
  implicit none

  integer :: i
  real(r_kind) :: tdiff1, tdiff2, tprod

  !errq = 999999.0d0
  !if(q<2) return

  tprod = 1.0d0
  do i=2,q
    tdiff1 = 0.0d0
    tdiff2 = 0.0d0
    tdiff1 = tdiff1 + sum(tau(1:i))  !(t[n]-t[n-i])
    tdiff2 = tdiff1 - tau(1)    !(t[n-1]-t[n-i])
    tprod = tprod*tdiff1/tdiff2 !tprod = tprod*(t[n]-t[n-i])/(t[n-1]-t[n-i])
  enddo

  errq(1:net_size) = dabs(-(1.0d0/l(2))*(en(1:net_size)/(1.0d0+tprod)))

  return
end subroutine errorq


!>
!! Calculate the errors for the order (q-1)
!!
!! The error is defined as:
!! \f[
!! \vec{E}_{n+1}(q-1) = - \frac{\prod \limits _{i=1}^{q-1} \xi _i }
!!                       {\ell _1 (q-1)}  \frac{h^q \vec{y}_{n+1}^{(q)}}{q!}
!! \f]
!!
!! @author D. Martin
!! @date   01.08.17
subroutine errorqm1
  implicit none

  integer :: qm1
  real(r_kind) :: xiprod

  qm1 = q-1

  errqm1 = 999999.0d0
  if(qm1<1) return

  xiprod = 1.0d0
  xiprod = xiprod*product(xi(1:qm1))           !xiprod = xiprod*xi[i]
  xiprod = xiprod/(l(2)*qm1)

  errqm1(1:net_size) = dabs(-xiprod*z(q+1,1:net_size)) !abs(-xiprod*zq[k])

  return
end subroutine errorqm1


!>
!! Calculate the errors for the order (q+1)
!!
!! The error is defined as:
!! \f[
!! \vec{E}_{n+1}(q+1) = \frac{-\xi_{q+1}(\vec{e}_{n+1}-Q_{n+1}\vec{e}_n) }
!!                      {(q+2)\ell _1 (q+1)\left[1+ \prod \limits_{i=2}^{q} \frac{t_{n+1}-t_{n+1-i}}
!!                      {t_n-t_{n+1-i}} \right] }
!! \f]
!!
!! @author D. Martin
!! @date   01.08.17
subroutine errorqp1
  implicit none

  integer :: qp1

  integer      :: i
  real(r_kind) :: xiqp1,tdiff1,tdiff2,helpprod,qval

  qp1 = q+1

  xiqp1 = 0.0d0
  xiqp1 = xiqp1 + sum(tau(1:qp1))
  xiqp1 = xiqp1/tau(1)

  helpprod = 1.0d0
  do i=2,q
    tdiff1 = 0.0d0
    tdiff2 = 0.0d0
    tdiff1 = tdiff1 + sum(tau(1:i))   !(t[n]-t[n-i])
    tdiff2 = tdiff1 - tau(1)            !(t[n-1]-t[n-i])
    helpprod = helpprod*tdiff1/tdiff2
  enddo
  helpprod = helpprod + 1.0d0

  call qfun(helpprod,qval)

  errqp1(1:net_size) = dabs(-xiqp1*(en(1:net_size)-qval*e(2,1:net_size))/((qp1+1)*l(2)*qp1*helpprod))

  return
end subroutine errorqp1


!> Calculate Helper function Q
!!
!! Q is defined as:
!! \f[
!! Q_{n+1} =\frac{C_{n+1}}{C_n}\left( \frac{h_{n+1}}{h_n}\right)^{q+1} \\
!! \f]
!!
!! @author D. Martin
!! @date   01.08.17
subroutine qfun(helpprod,qval)
  implicit none
  real(r_kind),intent(in)                     :: helpprod
  real(r_kind),intent(out)                    :: qval

  integer :: qp1

  integer      :: i
  real(r_kind) :: helpprod2, tdiff1, tdiff2, cfac1, cfac2, helpfac

  qp1 = q+1

  !!! 2nd helpprod for n-1
  helpprod2 = 1.0d0
  do i=2,q
    tdiff1 = 0.0d0
    tdiff2 = 0.0d0
    tdiff1 = tdiff1 + sum(tau(2:i+1))   !(t[n-1]-t[n-1-i] = t[n-1]-t[n-(i+1)])
    tdiff2 = tdiff1 - tau(2)            !(t[n-2]-t[n-1-i] = t[n-2]-t[n-(i+1))
    helpprod2 = helpprod2*tdiff1/tdiff2
  enddo
  helpprod2 = helpprod2 + 1.0d0

  call cfun(qp1,xi(1:histsize),helpprod,cfac1)
  call cfun(qp1,xi(1:histsize),helpprod2,cfac2)
  helpfac = cfac1/cfac2
  qval = helpfac*((tau(1)/tau(2))**qp1)

  return
end subroutine qfun

!> Calculate helper function C
!!
!! The function is defined as:
!! \f[
!! C_{n+1} = \frac{\prod \limits _{i=1}^{q} \xi _{i}}{(q+1)!}\left[1+ \prod
!!           \limits _{i=2}^{q} \frac{t_{n+1} - t_{n+1-i}}{t_n -t_{n+1-i}} \right].
!! \f]
!!
!! @author D. Martin
!! @date   01.08.17
subroutine cfun(qp1,xi_in,helpprod,cval)
  implicit none

  integer,intent(in)                          :: qp1
  real(r_kind),dimension(histsize),intent(in) :: xi_in
  real(r_kind),intent(in)                     :: helpprod
  real(r_kind),intent(out)                    :: cval

  !integer      :: j
  real(r_kind) :: facqp1,xiprod

  ! simple factorial calculation: (q+1)!
  facqp1 = factorial(qp1)
  !facqp1 = 1.0d0
  !do j=2,qp1
  !   facqp1 = facqp1*real(j)
  !enddo

  xiprod = 1.0d0
  xiprod = xiprod*product(xi_in(1:q))
  cval = xiprod*helpprod/facqp1

  return
end subroutine cfun
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!>
!! Adjusts the order, if possible and necessary.
!! This is done every qth step, or earlier in case the time
!! step was modified.
!!
!!
!! @author D. Martin
!! @date   01.08.17
subroutine shiftorder
  implicit none

  integer :: j, iback


  if(hqm1 .ge. hq .and. hqm1 .ge. hqp1) then

    el    = 0.0d0
    el(3) = 1.0d0
    do j=3,q
      el(j+1) = 1.0d0
      do iback=j,3,-1
        el(iback) = el(iback)*xi(j-2) + el(iback-1)
      enddo
    enddo

    ! subtract from Nordsieck vector for j=2,3,...,q-1
    do j=3,q
      z(j,1:net_size) = z(j,1:net_size) - z(q+1,1:net_size)*el(j)
    enddo
    q = q-1
    h = hqm1
    ! print *, "The order has decreased!",q,h
  else if(hqp1>hq .and. q<qmax) then
    q = q+1
    ! set (q+2)-th entry [=(q+1)-th order] of z to 0.0
    z(q+1,1:net_size) = 0.0d0
    h = hqp1
    ! print *, "The order has been increased!",q,h
  else
    h = hq
    ! print *, "The order stayed the same!",q,h
  endif

  !gear_cFactor = get_cFactor(q)

end subroutine shiftorder

!>
!! Function to set the conservative timestep factor depending on the current
!! order q.
!!
!! @note This function is not used!
!!
!! @author D. Martin
!! @date   01.08.17
function get_cFactor(q_val)
  implicit none

  integer, intent(in) :: q_val
     real(r_kind) :: get_cFactor

  get_cFactor = 0.075d0*(q_val-1) + 0.1d0
  return

  if(q_val==5) then
     get_cFactor = 0.30d0
  elseif(q_val==4) then
     get_cFactor = 0.25d0
  elseif(q_val==3) then
     get_cFactor = 0.20d0
  elseif(q_val==2) then
     get_cFactor = 0.15d0
  else
     get_cFactor = 0.10d0
  endif

end function get_cFactor

end module gear_module
