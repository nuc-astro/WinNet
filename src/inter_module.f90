!> @file inter_mod.f90
!!
!! The error file code for this file is ***W27***.
!! This file contains the module \ref inter_module

!>
!! @brief Module \ref inter_module with interpolation routines
!!
!! This routines are used to interpolate the hydrodynamic
!! quantaties (e.g., temperature, density,..)
!!
#include "macros.h"
module inter_module
  implicit none

  !
  ! Public and private fields and methods of the module (everything is public)
  !
  public:: &
      lininter_2D, cubinter_2D, get_indice_2D, &
      calc_derivative_2D, interp1d, inverse_interp1d
  private:: &
      lininter, cubinter, makimainter, akimainter, pchipinter

contains

!>
!! Linear interpolation in lin-log space.
!!
!! This routine assumes that the x and y arrays are given
!! in the form of {x_i, log(y_i)}.
!! An optional argument flin allows to control the return value:
!!  - if flin is absent: return exp(log(y)) = y;
!!  - if flin is present: return log(y).
!!
!! \b Edited:
!!       - 11.01.14
!! .
subroutine lininter(n,ref_array,ref_value,array,res,flin)
   implicit none

   integer,intent(in)                    :: n         !< array sizes
   real(r_kind),dimension(n),intent(in)  :: ref_array !< array of x-values {x_i}
   real(r_kind),dimension(n),intent(in)  :: array     !< array of log(y)-values {log(y_i)}
   real(r_kind),intent(in)               :: ref_value !< value of log(x)
   real(r_kind),intent(out)              :: res       !< interpolated value of y
   logical,intent(in),optional           :: flin      !< if present, return log(y) instead of y
   !
   real(r_kind)                          :: grad
   integer                               :: i,ii
   logical                               :: linear

   ! Default value
   if (.not. present(flin)) then
     linear = .false.
   else
     linear = flin
   end if

   ii = 1
   do i = 1,n
      if (ref_value .le. ref_array(i)) exit
      ii = i+1
   end do

   if (ii.eq.1) then
      res = array(1)
   else if (ii.gt.n) then
      res = array(n)
   else
      grad = (ref_value-ref_array(ii-1))/(ref_array(ii)-ref_array(ii-1))
      res = array(ii-1) + grad*(array(ii)-array(ii-1))
   end if

   if(linear) then
     return
   else
     res= dexp(res)
     return
   end if

end subroutine lininter



!> Interface for 1D interpolation routines.
!!
!! This routine is a wrapper for the 1D interpolation routines
!! lininter, cubinter, akimainter, makimainter, and pchipinter.
!! The interpolation type is controlled by the optional argument
!! itype. If itype is absent, \ref interp_mode is used.
!! The optional argument flin allows to control the return value:
!! - if flin is absent: return exp(log(y)) = y;
!! - if flin is present: return log(y).
!!
!! @author M. Reichert
!! @date 01.05.2023
subroutine interp1d(n,xp,xb,yp,res,flin,itype)
    use error_msg_class, only: raise_exception, int_to_str
    use parameter_class, only: interp_mode
    integer,intent(in)                    :: n    !< array sizes
    real(r_kind),dimension(n),intent(in)  :: xp   !< array of x-values
    real(r_kind),intent(in)               :: xb   !< value of x
    real(r_kind),dimension(n),intent(in)  :: yp   !< array of log(y)-values
    real(r_kind),intent(out)              :: res  !< interpolated value of y
    logical,intent(in),optional           :: flin !< if present, return log(y) instead of y
    integer,intent(in),optional           :: itype!< Type of interpolation (1: linear, 2: cubic)
    ! Internal variables
    integer :: interp_type
    logical :: linear

    ! Set default values
    if (.not. present(itype)) then
        interp_type = interp_mode
    else
        interp_type = itype
    end if

    if (.not. present(flin)) then
        linear = .false.
    else
        linear = flin
    end if

    ! Decide on desired interpolation type
    if (interp_type .eq. itype_LINEAR) then
        call lininter(n,xp,xb,yp,res,linear)
    else if (interp_type .eq. itype_CUBIC) then
        call cubinter(n,xp,xb,yp,res,linear)
    else if (interp_type .eq. itype_AKIMA) then
        call akimainter(n,xp,xb,yp,res,linear)
    else if (interp_type .eq. itype_MAKIMA) then
        call makimainter(n,xp,xb,yp,res,linear)
    else if (interp_type .eq. itype_PCHIP) then
        call pchipinter(n,xp,xb,yp,res,linear)
    else
        call raise_exception("Interpolation type not known. Got "//&
                             int_to_str(interp_type)//&
                             ". Try to change interp_mode to a valid value.",&
                             "interp1d",270004)
    end if

end subroutine interp1d



!>
!! Cubic interpolation in log-log space. Here, it is assumed that
!!
!! the x and y arrays are given in the form of {x_i, log(y_i)},
!! but before interpolating the x coordinate, it is transformed
!! to a log space: {log(x_i), log(y_i)}.
!! An optional argument flin allows to control the return value:
!!  - if flin is absent: return exp(log(y)) = y;
!!  - if flin is present: return log(y).
!!
!! \b Edited:
!!     - 14.01.14
!! .
subroutine cubinter (n,xp,xb,yp,res,flin)
   implicit none

   integer,intent(in)                    :: n    !< array sizes
   real(r_kind),dimension(n),intent(in)  :: xp   !< array of x-values
   real(r_kind),intent(in)               :: xb   !< value of x
   real(r_kind),dimension(n),intent(in)  :: yp   !< array of log(y)-values
   real(r_kind),intent(out)              :: res  !< interpolated value of y
   logical,intent(in),optional           :: flin !< if present, return log(y) instead of y
   !
   integer      :: ii,k
   real(r_kind) :: x,x1,x2,x3,x4,x12,x13,x14,x23,x24,x34
   logical      :: linear

   if(present(flin)) then
      linear= flin
   else
      linear= .false.
   end if

   x= xb
   do ii= 1,n
     if (xb .le. xp(ii)) exit
   end do
   if (x.le.xp(1)) then
     res= yp(1)
   else if (x.ge.xp(n)) then
     res= yp(n)
   else if(ii.le.2) then ! linear interpolation in lin-log space
     x1=  x - xp(1)
     x2=  x - xp(2)
     x12= x2-x1
     res= (yp(1)*x2 - yp(2)*x1)/x12
   else if(ii.eq.n) then ! linear interpolation in log-log space
     x1=  log(x/xp(n-1))
     x2=  log(x/xp(n))
     x12= x2-x1
     res= (yp(n-1)*x2 - yp(n)*x1)/x12
   else                 ! cubic interpolation in log-log space
     k= min(max(ii-2,2),n-3)
     x= log(x)
     x1= x - log(xp(k))
     x2= x - log(xp(k+1))
     x3= x - log(xp(k+2))
     x4= x - log(xp(k+3))

     x12= x2-x1
     x13= x3-x1
     x14= x4-x1
     x23= x3-x2
     x24= x4-x2
     x34= x4-x3

     res= yp(k)  *x2*x3*x4/(x12*x13*x14) &
        - yp(k+1)*x1*x3*x4/(x12*x23*x24) &
        + yp(k+2)*x1*x2*x4/(x13*x23*x34) &
        - yp(k+3)*x1*x2*x3/(x14*x24*x34)
   endif

   if(linear) then
     return
   else
     res= dexp(res)
     return
   end if
end subroutine cubinter


!> PCHIP interpolation in lin-log space.
!!
!! The Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)
!! interpolation is a method for interpolating data points where
!! every section between points is monoton. The interpolation
!! avoids overshoots and undershoots and is close to the linear
!! interpolation, but having the first derivative continuous.
!!
!! @see F. N. Fritsch and R. E. Carlson, Monotone Piecewise Cubic
!!      Interpolation, SIAM Journal on Numerical Analysis, 17(2),
!!      1980, pp. 238-246.
!!      [Steffen 1990](https://ui.adsabs.harvard.edu/abs/1990A%26A...239..443S/abstract)
!! @author: M. Reichert
!! @date: 01.05.2023
subroutine pchipinter(n,xp,xb,yp,res,flin)
    implicit none
    integer,intent(in)                    :: n    !< array sizes
    real(r_kind),dimension(n),intent(in)  :: xp   !< array of x-values
    real(r_kind),intent(in)               :: xb   !< value of x
    real(r_kind),dimension(n),intent(in)  :: yp   !< array of log(y)-values
    real(r_kind),intent(out)              :: res  !< interpolated value of y
    logical,intent(in),optional           :: flin !< if present, return log(y) instead of y
    ! Internal variables
    integer :: ii, i
    real(r_kind) :: a, b, c, d
    real(r_kind) :: hi, him1, hip1
    real(r_kind) :: si, sim1, sip1
    real(r_kind) :: pi, pip1
    real(r_kind) :: yi_prime, yip1_prime
    real(r_kind) :: x1, x2
    real(r_kind) :: delt
    logical      :: linear

    if(present(flin)) then
       linear= flin
    else
       linear= .false.
    end if

    ! Find the index for the interpolation interval
    do ii= 1,n
        if (xb .le. xp(ii)) exit
    end do
    ii = ii-1

    ! Extrapolation, value is smaller than grid
    if (xb.le.xp(1)) then
        res= yp(1)
    ! Extrapolation, value is larger than grid
    else if (xb.ge.xp(n)) then
        res= yp(n)
    ! Value is at left border
    else if(ii.lt.2) then
        ! Interpolate linear
        x1=  xb - xp(ii)
        x2=  xb - xp(ii+1)
        res= (yp(ii)*x2 - yp(ii+1)*x1)/(x2-x1)
    ! Value is at left border
    else if(ii.gt.n-2) then
        ! Interpolate linear
        x1=  xb - xp(ii)
        x2=  xb - xp(ii+1)
        res= (yp(ii)*x2 - yp(ii+1)*x1)/(x2-x1)
    else
        ! Normal case, use pchip interpolation
        ! Values in interval
        hi = xp(ii+1)-xp(ii)
        si = (yp(ii+1)-yp(ii))/hi
        ! Values to the left
        him1 = xp(ii)-xp(ii-1)
        sim1 = (yp(ii)-yp(ii-1))/him1
        ! Values to the right
        hip1 = xp(ii+2)-xp(ii+1)
        sip1 = (yp(ii+2)-yp(ii+1))/hip1

        pi   = (hi*sim1+him1*si)/(hi+him1)
        pip1 = (hip1*si+hi*sip1)/(hip1+hi)

        ! Restrict slopes to ensure monotonicity
        if (si*sim1 .le. 0) then
            yi_prime = 0
        elseif ((dabs(pi)>2*dabs(sim1)) .or. (dabs(pi)>2*dabs(si))) then
            yi_prime = 2*sign(1.d0,si)*min(dabs(sim1),dabs(si))
        else
            yi_prime=pi
        end if

        if (sip1*si .le. 0) then
            yip1_prime = 0
        elseif ((dabs(pip1)>2*dabs(si)) .or. (dabs(pip1)>2*dabs(sip1))) then
            yip1_prime = 2*sign(1.d0,si)*min(dabs(sip1),dabs(si))
        else
            yip1_prime=pip1
        end if

        a = (yi_prime+yip1_prime-2*si)/hi**2
        b = (3*si-2*yi_prime-yip1_prime)/hi
        c = yi_prime
        d = yp(ii)

        delt = xb-xp(ii)
        ! Calculate interpolated value
        res = a*delt**3 + b*delt**2 + c*delt + d

    end if

    if(linear) then
      return
    else
      res= dexp(res)
      return
    end if

end subroutine pchipinter



!> Interpolation using prescription of Akima
!!
!! The interpolation is a spline interpolation that tries
!! to avoid overshoots.
!!
!! Here, it is assumed that
!! the x and y arrays are given in the form of {x_i, log(y_i)},
!! but before interpolating the x coordinate, it is transformed
!! to a log space: {log(x_i), log(y_i)}.
!! An optional argument flin allows to control the return value:
!!  - if flin is absent: return exp(log(y)) = y;
!!  - if flin is present: return log(y).
!! .
!!
!! @see [Akima, Journal of the ACM, 17: 589–602, 1970](https://web.archive.org/web/20201218210437/http://www.leg.ufpr.br/lib/exe/fetch.php/wiki:internas:biblioteca:akima.pdf),
!!      [Wikipedia: Akima Spline](https://en.wikipedia.org/wiki/Akima_spline)
!! @author M. Reichert
!! @date 30.04.2023
subroutine akimainter(n,xp,xb,yp,res,flin)
    implicit none
    integer,intent(in)                    :: n    !< array sizes
    real(r_kind),dimension(n),intent(in)  :: xp   !< array of x-values
    real(r_kind),intent(in)               :: xb   !< value of x
    real(r_kind),dimension(n),intent(in)  :: yp   !< array of log(y)-values
    real(r_kind),intent(out)              :: res  !< interpolated value of y
    logical,intent(in),optional           :: flin !< if present, return log(y) instead of y
    ! Internal variables
    real(r_kind) :: xval(6), yval(6), mval(5), sval(2)
    integer :: ii, i
    real(r_kind) :: a, b, c, d
    real(r_kind) :: x1,x2
    real(r_kind) :: denom
    logical      :: linear

    if(present(flin)) then
       linear= flin
    else
       linear= .false.
    end if

    xval = 0.0d0
    yval = 0.0d0
    mval = 0.0d0
    sval = 0.0d0


    ! Find the index for the interpolation interval
    do ii= 1,n
        if (xb .le. xp(ii)) exit
    end do
    ii = ii-1

    ! Extrapolation, value is smaller than grid
    if (xb.le.xp(1)) then
        res= yp(1)
    ! Extrapolation, value is larger than grid
    else if (xb.ge.xp(n)) then
        res= yp(n)
    ! Value is at left border
    else if(ii.lt.3) then
        ! Interpolate linear
        x1=  xb - xp(ii)
        x2=  xb - xp(ii+1)
        res= (yp(ii)*x2 - yp(ii+1)*x1)/(x2-x1)
    ! Value is at left border
    else if(ii.gt.n-3) then
        ! Interpolate linear
        x1=  xb - xp(ii)
        x2=  xb - xp(ii+1)
        res= (yp(ii)*x2 - yp(ii+1)*x1)/(x2-x1)
    else
        ! Normal case, use Akima interpolation

        ! Save the required values of the grid for the interpolation
        xval(1:6) = xp(ii-2:ii+3)
        yval(1:6) = yp(ii-2:ii+3)
        ! Calculate the derivative
        mval(1:5) = (yval(2:6) - yval(1:5)) / (xval(2:6) - xval(1:5))

        ! Get si and si+1
        i = 3
        sval(1) = ((dabs(mval(i+1)-mval(i))*mval(i-1) + dabs(mval(i-1)-mval(i-2))*mval(i)))
        denom   = (dabs(mval(i+1)-mval(i)) + dabs(mval(i-1)-mval(i-2)))
        ! Avoid division by something very small
        if (denom .gt. 1e-8) then
            sval(1) = sval(1)/denom
        else
            sval(1) = (mval(i-1)+mval(i))/2d0
        end if
        sval(2) = ((dabs(mval(i+2)-mval(i+1))*mval(i) + dabs(mval(i)-mval(i-1))*mval(i+1)))
        denom   = (dabs(mval(i+2)-mval(i+1)) + dabs(mval(i)-mval(i-1)))
        if (denom .gt. 1e-8) then
            sval(2) = sval(2)/denom
        else
            sval(2) = (mval(i)+mval(i+1))/2d0
        end if
        ! Get the coefficients
        a = yval(i)
        b = sval(1)
        c = (3.0d0 * mval(i) - 2.0d0 * sval(1) - sval(2)) / (xval(i+1) - xval(i))
        d = (sval(1) + sval(2) - 2.0d0 * mval(i)) / (xval(i+1) - xval(i))**2
        ! Calculate the results
        res = a + b * (xb - xval(i)) + c * (xb - xval(i))**2 + d * (xb - xval(i))**3
    end if

    if(linear) then
      return
    else
      res= dexp(res)
      return
    end if

end subroutine akimainter



!> Interpolation using modified prescription of Akima
!!
!! The interpolation is a spline interpolation that tries
!! to avoid overshoots. It is doing this slightly more
!! radically than the Akima interpolation, but not as
!! strong as the pchip interpolation.
!!
!! Here, it is assumed that
!! the x and y arrays are given in the form of {x_i, log(y_i)},
!! but before interpolating the x coordinate, it is transformed
!! to a log space: {log(x_i), log(y_i)}.
!! An optional argument flin allows to control the return value:
!!  - if flin is absent: return exp(log(y)) = y;
!!  - if flin is present: return log(y).
!! .
!!
!! @see https://blogs.mathworks.com/cleve/2019/04/29/makima-piecewise-cubic-interpolation/
!! @author M. Reichert
!! @date 30.04.2023
subroutine makimainter(n,xp,xb,yp,res,flin)
    implicit none
    integer,intent(in)                    :: n    !< array sizes
    real(r_kind),dimension(n),intent(in)  :: xp   !< array of x-values
    real(r_kind),intent(in)               :: xb   !< value of x
    real(r_kind),dimension(n),intent(in)  :: yp   !< array of log(y)-values
    real(r_kind),intent(out)              :: res  !< interpolated value of y
    logical,intent(in),optional           :: flin !< if present, return log(y) instead of y
    ! Internal variables
    real(r_kind) :: xval(6), yval(6), mval(5), sval(2)
    integer :: ii, i
    real(r_kind) :: a, b, c, d
    real(r_kind) :: x1,x2
    real(r_kind) :: denom
    logical      :: linear

    if(present(flin)) then
       linear= flin
    else
       linear= .false.
    end if

    xval = 0.0d0
    yval = 0.0d0
    mval = 0.0d0
    sval = 0.0d0


    ! Find the index for the interpolation interval
    do ii= 1,n
        if (xb .le. xp(ii)) exit
    end do
    ii = ii-1

    ! Extrapolation, value is smaller than grid
    if (xb.le.xp(1)) then
        res= yp(1)
    ! Extrapolation, value is larger than grid
    else if (xb.ge.xp(n)) then
        res= yp(n)
    ! Value is at left border
    else if(ii.lt.3) then
        ! Interpolate linear
        x1=  xb - xp(ii)
        x2=  xb - xp(ii+1)
        res= (yp(ii)*x2 - yp(ii+1)*x1)/(x2-x1)
    ! Value is at left border
    else if(ii.gt.n-3) then
        ! Interpolate linear
        x1=  xb - xp(ii)
        x2=  xb - xp(ii+1)
        res= (yp(ii)*x2 - yp(ii+1)*x1)/(x2-x1)
    else
        ! Normal case, use Akima interpolation

        ! Save the required values of the grid for the interpolation
        xval(1:6) = xp(ii-2:ii+3)
        yval(1:6) = yp(ii-2:ii+3)
        ! Calculate the derivative
        mval(1:5) = (yval(2:6) - yval(1:5)) / (xval(2:6) - xval(1:5))

        ! Get si and si+1
        i = 3
        sval(1) = ((dabs(mval(i+1)-mval(i)) + dabs(mval(i+1)+mval(i))/2d0)*mval(i-1) + &
                   (dabs(mval(i-1)-mval(i-2)) + dabs(mval(i-1)+mval(i-2))/2d0)*mval(i))

        denom   = ((dabs(mval(i+1)-mval(i)) + dabs(mval(i+1)+mval(i))/2d0) + &
                   (dabs(mval(i-1)-mval(i-2)) + dabs(mval(i-1)+mval(i-2))/2d0))
        ! Avoid division by something very small
        if (denom .gt. 1e-8) then
            sval(1) = sval(1)/denom
        else
            sval(1) = (mval(i-1)+mval(i))/2d0
        end if
        sval(2) = (((dabs(mval(i+2)-mval(i+1)) + dabs(mval(i+2)+mval(i+1))/2d0)*mval(i) + &
                    (dabs(mval(i)-mval(i-1))   + dabs(mval(i)+mval(i-1))/2d0) *mval(i+1)))
        denom   = ((dabs(mval(i+2)-mval(i+1)) + dabs(mval(i+2)+mval(i+1))/2d0)  + &
                   (dabs(mval(i)-mval(i-1))   + dabs(mval(i)+mval(i-1))/2d0))
        if (denom .gt. 1e-8) then
            sval(2) = sval(2)/denom
        else
            sval(2) = (mval(i)+mval(i+1))/2d0
        end if
        ! Get the coefficients
        a = yval(i)
        b = sval(1)
        c = (3.0d0 * mval(i) - 2.0d0 * sval(1) - sval(2)) / (xval(i+1) - xval(i))
        d = (sval(1) + sval(2) - 2.0d0 * mval(i)) / (xval(i+1) - xval(i))**2
        ! Calculate the results
        res = a + b * (xb - xval(i)) + c * (xb - xval(i))**2 + d * (xb - xval(i))**3
    end if

    if(linear) then
      return
    else
      res= dexp(res)
      return
    end if

end subroutine makimainter


!> Get indices of the grid for given x- and y- values
!!
!! @author: M. Reichert
!! @date: 06.07.22
subroutine get_indice_2D(xval,yval,x,y,x_dim,y_dim,indices)
  implicit none
  integer,intent(in)                       :: x_dim, y_dim !< Dimension of x- and y-
  real(r_kind),intent(in)                  :: xval,yval    !< Function values
  real(r_kind),dimension(x_dim),intent(in) :: x            !< X-grid values
  real(r_kind),dimension(y_dim),intent(in) :: y            !< Y-grid values
  integer,dimension(2,2),intent(out)       :: indices      !< Index in the grid
  ! Internal variables
  integer :: i !< Loop variables

  ! Index in x-direction
  do i=2,x_dim-1
    if (xval .lt. x(i)) exit
  end do
  indices(1,1)=i-1
  indices(1,2)=i

  ! Index in y-direction
  do i=2,y_dim-1
    if (yval .lt. y(i)) exit
  end do
  indices(2,1)=i-1
  indices(2,2)=i

end subroutine get_indice_2D


!> Calculate derivatives at all points in a given grid.
!!
!! These derivatives are necessary for the cubic interpolation.
!! The subroutine returns the derivative df/dx, df/dy, and ddf/dxdy.
!! They are second order accurate, i.e.,
!! \f[ df/dx = \frac{h_1^2 f(x+h_2) + (h_2^2-h_1^2)f(x)-h_2^2 f(x-h_1) }{h_1 h_2 (h_1+h_2)} \f]
!!
!! @author: M. Reichert
!! @date 06.07.22
subroutine calc_derivative_2D(f,x,y,dfx,dfy,dfxy,x_dim,y_dim)
  implicit none
  integer,intent(in)                                :: x_dim, y_dim !< Dimension of x- and y-
  real(r_kind),dimension(x_dim, y_dim),intent(in)   :: f            !< Function values
  real(r_kind),dimension(x_dim),intent(in)          :: x            !< x values
  real(r_kind),dimension(y_dim),intent(in)          :: y            !< y values
  real(r_kind),dimension(x_dim,y_dim),intent(inout) :: dfx          !< Derivative with respect to x
  real(r_kind),dimension(x_dim,y_dim),intent(inout) :: dfy          !< Derivative with respect to y
  real(r_kind),dimension(x_dim,y_dim),intent(inout) :: dfxy         !< Derivative with respect to x and y
  ! Internal variables
  integer         :: i, j   !< Loop variables
  real(r_kind)    :: h1,h2  !< steps

  ! Calculate derivatives for every point
  do i=1,x_dim
    do j=1,y_dim
      ! X-direction
      if ((i .gt. 1)  .and. (i .lt. x_dim) ) then
        h1 = abs(x(i)  -x(i-1))
        h2 = abs(x(i+1)-x(i))
        dfx(i,j) = h1**2*f(i+1,j)+(h2**2-h1**2)*f(i,j)-h2**2*f(i-1,j)
        dfx(i,j) = dfx(i,j) / (h1*h2*(h1+h2))
      elseif (i .eq. 1) then
        h2 = abs(x(i+1)-x(i))
        dfx(i,j) = (f(i+1,j)-f(i,j)) / (h2)
      elseif (i .eq. x_dim) then
        h1 = abs(x(i)-x(i-1))
        dfx(i,j) = (f(i,j)-f(i-1,j)) / (h1)
      end if

      ! Y-direction
      if ((j .gt.1) .and. (j .lt. y_dim)) then
        h1 = abs(y(j)  -y(j-1))
        h2 = abs(y(j+1)-y(j))
        dfy(i,j) = h1**2*f(i,j+1)+(h2**2-h1**2)*f(i,j)-h2**2*f(i,j-1)
        dfy(i,j) = dfy(i,j) / (h1*h2*(h1+h2))
      elseif (j .eq. 1) then
        h2 = abs(y(j+1)-y(j))
        dfy(i,j) = (f(i,j+1)-f(i,j)) / (h2)
      elseif (j .eq. y_dim) then
        h1 = abs(y(j)-y(j-1))
        dfy(i,j) = (f(i,j)-f(i,j-1)) / (h1)
      end if
    end do
  end do

  ! Mixed partial derivatives
  do i=1, x_dim
    do j=1, y_dim
      if ((j .gt. 1)  .and. (j .lt. y_dim) ) then
        h1 = abs(y(j)  -y(j-1))
        h2 = abs(y(j+1)-y(j))
        dfxy(i,j) = h1**2*dfx(i,j+1)+(h2**2-h1**2)*dfx(i,j)-h2**2*dfx(i,j-1)
        dfxy(i,j) = dfxy(i,j) / (h1*h2*(h1+h2))
      elseif (j .eq. 1) then
        h2 = abs(y(j+1)-y(j))
        dfxy(i,j) = (dfx(i,j+1)-dfx(i,j)) / (h2)
      elseif (j .eq. y_dim) then
        h1 = abs(y(j)-y(j-1))
        dfxy(i,j) = (dfx(i,j)-dfx(i,j-1)) / (h1)
      end if
    end do
  end do

end subroutine calc_derivative_2D


!> Interpolate linearly on 2 Dimensional grid
!!
!! Bilinear interpolation
!!
!! @see weak_inter
!!
!! @author: M. Reichert
!! @date: 06.07.22
!! .
function lininter_2D(xin,yin,x,y,f,x_dim,y_dim,indice,extrapolation) result (wr)
  use error_msg_class, only: raise_exception
  implicit none
  integer,intent(in)                             :: x_dim, y_dim   !< Dimension of x- and y-
  real(r_kind),intent(in)                        :: xin,yin        !< Function input values
  real(r_kind),dimension(x_dim),intent(in)       :: x              !< X values
  real(r_kind),dimension(y_dim),intent(in)       :: y              !< Y values
  real(r_kind),dimension(x_dim,y_dim),intent(in) :: f              !< Function and derivatives
  integer,optional,intent(in)                    :: extrapolation  !< 1: None; 2: constant
  integer,dimension(2,2),optional,intent(in)     :: indice         !< Index where to interpolate
  real(r_kind)                                   :: wr             !< Interpolated value (output)
  ! Internal Variables
  integer,dimension(2,2)      :: ind            !< Index where to interpolate
  real(r_kind)                :: xval,yval      !< Internal function input values
  real(r_kind)                :: x1, x2, y1, y2 !< Corner values of the grid
  integer                     :: extr           !< Internal extr. flag
  real(r_kind),dimension(2,2) :: weights        !< Weights for interpolation

  ! Store variables in internal variable in order to be able to change them
  xval=xin; yval=yin

  ! Get the index
  if (.not. present(indice)) then
    call get_indice_2D(xval,yval,x,y,x_dim,y_dim,ind)
  else
    ind(:,:) = indice(:,:)
  end if

  ! Check that an extrapolation value is present
  if (.not. present(extrapolation)) extr=1
  if (present(extrapolation))       extr=extrapolation

  ! Check if extrapolation is necessary
  if ((xval .lt. x(1)) .or. (yval .lt. y(1))) then
    select case(extr)
      case (1)
        call raise_exception("Value for interpolation was out of bounds.",&
                             "inter_module",270003)
      case (2)
        xval= max(x(1),xval)
        yval= max(y(1),yval)
      case default
        call raise_exception("Unknown interpolation type.",&
                             "inter_module",270004)
    end select
  elseif ((xval .gt. x(x_dim)) .or. (yval .gt. y(y_dim))) then
    select case(extr)
    case (1)
      call raise_exception("Value for interpolation was out of bounds.",&
                           "inter_module",270003)
    case (2)
      xval= min(x(x_dim),xval)
      yval= min(y(y_dim),yval)
    case default
      call raise_exception("Unknown interpolation type.",&
                           "inter_module",270004)
    end select
  end if

  ! Define x and y border points
  x1 = x(ind(1,1))
  x2 = x(ind(1,2))
  y1 = y(ind(2,1))
  y2 = y(ind(2,2))

  ! Weights for linear interpolation
  weights(1,1) = ((x2-xval)*(y2-yval))/((x2-x1)*(y2-y1))
  weights(1,2) = ((x2-xval)*(yval-y1))/((x2-x1)*(y2-y1))
  weights(2,1) = ((xval-x1)*(y2-yval))/((x2-x1)*(y2-y1))
  weights(2,2) = ((xval-x1)*(yval-y1))/((x2-x1)*(y2-y1))

  ! Get the interpolated value
  wr = weights(1,1)*f(ind(1,1),ind(2,1)) + &
       weights(1,2)*f(ind(1,1),ind(2,2)) + &
       weights(2,1)*f(ind(1,2),ind(2,1)) + &
       weights(2,2)*f(ind(1,2),ind(2,2))

   return

end function lininter_2D




!> Cubic interpolate on 2-Dimensional grid
!!
!! This interpolation is a bicubic interpolation. For values outside the grid,
!! either a constant extrapolation is performed or an error is raised.
!!
!! @see weak_inter_bilinear, weak_inter, readweak_logft
!!
!! @author: Moritz Reichert
!! @date:  14.03.22
!!
!! @see [Bicubic interpolation](https://en.wikipedia.org/wiki/Bicubic_interpolation)
!! .
function cubinter_2D(xin,yin,x,y,f,dfx,dfy,dfxy,x_dim,y_dim,indice,extrapolation) result (wr)
  use error_msg_class, only: raise_exception
  implicit none
  integer,intent(in)                              :: x_dim, y_dim   !< Dimension of x- and y-
  real(r_kind),intent(in)                         :: xin,yin        !< Function input values
  real(r_kind),dimension(x_dim),intent(in)        :: x              !< X values
  real(r_kind),dimension(y_dim),intent(in)        :: y              !< Y values
  real(r_kind),dimension(x_dim,y_dim),intent(in)  :: f,dfx,dfy,dfxy !< Function and derivatives
  integer,optional,intent(in)                     :: extrapolation  !< 1: None; 2: constant
  integer,dimension(2,2),optional,intent(in)      :: indice         !< Index where to interpolate
  real(r_kind)                                    :: wr             !< Interpolated value (output)
  ! Internal Variables
  real(r_kind)                  :: xval,yval      !< Internal function input values
  real(r_kind),dimension(16,16) :: A_inv          !< Inverse matrix to solve eq.
  real(r_kind),dimension(16)    :: vec,alphas     !< Vectors for lin. eq.
  integer,dimension(2,2)        :: ind            !< Index where to interpolate
  real(r_kind),dimension(4,4)   :: alphas_r       !< Reshaped alpha coefficients
  real(r_kind)                  :: x1, x2, y1, y2 !< Corner values of the grid
  real(r_kind)                  :: xbar,ybar      !< Helper variables
  integer                       :: i,j            !< Loop variable
  integer                       :: extr           !< Internal value for extr. flag

  ! Store variables in internal variable in order to be able to change them
  xval=xin; yval=yin

  ! Get the index
  if (.not. present(indice)) then
    call get_indice_2D(xval,yval,x,y,x_dim,y_dim,ind)
  else
    ind(:,:) = indice(:,:)
  end if

  ! Check that an extrapolation value is present
  if (.not. present(extrapolation)) extr=2
  if (present(extrapolation))       extr=extrapolation

  ! Check if extrapolation is necessary
  if ((xval .lt. x(1)) .or. (yval .lt. y(1))) then
    select case(extr)
      case (1)
        call raise_exception("Value for interpolation was out of bounds.",&
                             "inter_module",270003)
      case (2)
        xval= max(x(1),xval)
        yval= max(y(1),yval)
      case default
        call raise_exception("Unknown interpolation type.",&
                             "inter_module",270004)
    end select
  elseif ((xval .gt. x(x_dim)) .or. (yval .gt. y(y_dim))) then
    select case(extr)
    case (1)
      call raise_exception("Value for interpolation was out of bounds.",&
                           "inter_module",270003)
    case (2)
      xval= min(x(x_dim),xval)
      yval= min(y(y_dim),yval)
    case default
      call raise_exception("Unknown interpolation type.",&
                           "inter_module",270004)
    end select
  end if

  ! Define x and y border points
  x1 = x(ind(1,1))
  x2 = x(ind(1,2))
  y1 = y(ind(2,1))
  y2 = y(ind(2,2))

  ! Define xbar and ybar
  xbar = (xval-x1)/(x2-x1)
  ybar = (yval-y1)/(y2-y1)

  ! Define helper inverse matrix to solve equation for determining
  ! alpha coefficients later
  A_inv = reshape( (/ 1, 0,-3, 2, 0, 0, 0, 0,-3, 0, 9,-6, 2, 0,-6, 4, &
                      0, 0, 3,-2, 0, 0, 0, 0, 0, 0,-9, 6, 0, 0, 6,-4, &
                      0, 0, 0, 0, 0, 0, 0, 0, 3, 0,-9, 6,-2, 0, 6,-4, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-6, 0, 0,-6, 4, &
                      0, 1,-2, 1, 0, 0, 0, 0, 0,-3, 6,-3, 0, 2,-4, 2, &
                      0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 3,-3, 0, 0,-2, 2, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 3,-6, 3, 0,-2, 4,-2, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 2,-2, &
                      0, 0, 0, 0, 1, 0,-3, 2,-2, 0, 6,-4, 1, 0,-3, 2, &
                      0, 0, 0, 0, 0, 0, 3,-2, 0, 0,-6, 4, 0, 0, 3,-2, &
                      0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 3,-2, 1, 0,-3, 2, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 2, 0, 0, 3,-2, &
                      0, 0, 0, 0, 0, 1,-2, 1, 0,-2, 4,-2, 0, 1,-2, 1, &
                      0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 2,-2, 0, 0,-1, 1, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1, 0, 1,-2, 1, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1  &
                      /), (/16, 16/) ) !< Helper Matrix

  ! Values for cubic interpolation
  vec  =  (/f(ind(1,1),ind(2,1)),f(ind(1,2),ind(2,1)),f(ind(1,1),ind(2,2)),f(ind(1,2),ind(2,2)),&
          (x2-x1)*dfx(ind(1,1),ind(2,1)),(x2-x1)*dfx(ind(1,2),ind(2,1)),(x2-x1)*dfx(ind(1,1),ind(2,2)),(x2-x1)*dfx(ind(1,2),ind(2,2)),&
          (y2-y1)*dfy(ind(1,1),ind(2,1)),(y2-y1)*dfy(ind(1,2),ind(2,1)),(y2-y1)*dfy(ind(1,1),ind(2,2)),(y2-y1)*dfy(ind(1,2),ind(2,2)),&
          (x2-x1)*(y2-y1)*dfxy(ind(1,1),ind(2,1)),(x2-x1)*(y2-y1)*dfxy(ind(1,2),ind(2,1)),&
          (x2-x1)*(y2-y1)*dfxy(ind(1,1),ind(2,2)),(x2-x1)*(y2-y1)*dfxy(ind(1,2),ind(2,2))/)

  ! Solve lin. equation to get alpha coefficients
  alphas = matmul(A_inv,vec)
  ! Reshape to 4 x 4
  alphas_r = reshape(alphas, (/4,4/))

  ! Calculate the interpolation result
  wr = 0
  do i=1,4
    do j=1,4
      wr = wr + alphas_r(i,j)*xbar**(i-1)*ybar**(j-1)
    end do
  end do

end function cubinter_2D





!>
!! The inverse of the 1D interpolation function.
!!
!! Finds the inverse of cubinter and returns the value of x.
!! This is done by using a Newton-Raphson method.
!!
!! @note The result depends on xb, because x_i does not necessarily have
!!       to be monotonic. The result is the closest value to xb
!!
!! @author M. Reichert
!!
!! @see timestep_module::restrict_timestep
!!
!! \b Edited:
!!        - 24.11.15
!!        - 30.01.23
!! .
subroutine inverse_interp1d (n,xp,xb,yp,res,flin,itype)
   use parameter_class, only: interp_mode
   implicit none

   integer,intent(in)                    :: n    !< array sizes
   real(r_kind),dimension(n),intent(in)  :: xp   !< array of log(x)-values
   real(r_kind),intent(in)               :: xb   !< value of x
   real(r_kind),dimension(n),intent(in)  :: yp   !< array of y-values
   real(r_kind),intent(inout)            :: res  !< interpolated value of y
   logical,intent(in),optional           :: flin !< flag for linear interpolation
   integer,intent(in),optional           :: itype !< interpolation type
   !
   integer                   :: ii
   integer, parameter        :: maxiter=50
   real(r_kind)              :: x1,x2
   real(r_kind)              :: m,y1,y2,xnew,b
   real(r_kind),parameter    :: diff=1d-4
   real(r_kind),parameter    :: tol =1d-10
   logical                   :: converged
   logical                   :: linear
   integer                   :: interp_type

     ! Set default values
    if (.not. present(itype)) then
       interp_type = interp_mode
    else
       interp_type = itype
    end if

    if (.not. present(flin)) then
       linear = .false.
    else
       linear = flin
    end if

    x1 = res
    x2=x1+diff
    converged = .false.

    inner_loop: do ii=1,maxiter
       ! Calculate slope
       call interp1d(n,yp,x1,xp,y1,linear,interp_type)
       call interp1d(n,yp,x2,xp,y2,linear,interp_type)
       y1 = y1-xb
       y2 = y2-xb
       if (x1 .ne. x2) then
            m  = (y1-y2)/(x1-x2)
       else
            converged = .False.
            exit inner_loop
       endif

        ! calculate new x value
        b  = y1-m*x1

        ! Exit if the slope is zero
        if (m .eq. 0) then
            converged = .False.
            ! Give a very large negative number
            exit
        end if

        xnew = -b/m

        ! Exit if converged
        if (abs(x1-xnew) .lt. tol) then
            converged = .True.
            exit inner_loop
        end if

        ! Reshuffle variables if it is not converged
        x1 = xnew
        x2 = x1+diff
    end do inner_loop

    ! Take NR value if it converged
    if (converged) then
        res = xnew
    else
        res = res
    end if

   return
end subroutine inverse_interp1d

end module inter_module
