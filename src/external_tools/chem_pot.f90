!> @file chem_pot.f90
!!
!!
!! The error file code for this file is ***W13***.
!! @brief Subroutine \ref chempot to compute chemical potentials of electrons and positrons
!!
!! This file originates
!! from [Cococubed](https://cococubed.com/code_pages/chemical_potential.shtml).
!! The chemical potentials are only used if \ref parameter_class::use_timmes_mue is
!! enabled or heating is turned on.

!>
!! Given a temperature temp [K], density den [g/cm**3], and a composition
!! characterized by y_e, this routine returns the electron and positron
!! chemical potentials. derivatives with respect to temperature and
!! density are computed, but not returned.
!!
!! References: <a href="http://adsabs.harvard.edu/abs/1999ApJS..125..277T">Timmes & Arnett ApJS 1999</a>
!!
!! \b Edited:
!! - 02.02.21 - MR - gave intent(in/out) and let it crash when it does not work
!! .
#include "../macros.h"
subroutine chempot(temp,den,ye,etaele,etapos)
  use parameter_class, only: unit,chem_pot_file
  use error_msg_class, only: num_to_str, raise_exception
  use file_handling_class
  implicit none
  save

!..declare the pass
  real(r_kind),intent(in)  :: temp   !< Temperature [K]
  real(r_kind),intent(in)  :: den    !< Density [g/ccm]
  real(r_kind),intent(in)  :: ye     !< Electron fraction
  real(r_kind),intent(out) :: etaele !< Electron chemical potential
  real(r_kind),intent(out) :: etapos !< Positron chemical potential

!..declare local variables
!..for the tables
  integer          i,j,imax,jmax
  parameter        (imax = 271, jmax = 101)

!..for the density and temperature table
  real(r_kind) d(imax),t(jmax)

!..for chemical potential table
  real(r_kind) ef(imax,jmax),efd(imax,jmax),                   &
       eft(imax,jmax),efdt(imax,jmax)

!..for storing the differences
  real(r_kind) dt_sav(jmax),                                   &
       dti_sav(jmax),                                  &
       dd_sav(imax),                                   &
       ddi_sav(imax)

!..for the interpolations
  integer          iat,jat
  real(r_kind) tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi,          &
       tsav,dsav,deni,tempi,kt,ktinv,beta,dbetadt
  real(r_kind) dth,dti,dd,ddi,xt,xd,mxt,mxd,                   &
       si0t,si1t,si0mt,si1mt,                          &
       si0d,si1d,si0md,si1md,                          &
       dsi0t,dsi1t,dsi0mt,dsi1mt,                      &
       dsi0d,dsi1d,dsi0md,dsi1md,                      &
       x,z,din,fi(16),                                 &
       xpsi0,xdpsi0,xpsi1,xdpsi1,h3,                   &
       w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md,            &
       detadd,detadt,detapdd,detapdt

!..physical constants and parameters
  real(r_kind) mecc,positron_start

!..for initialization
  integer          ifirst
  data             ifirst/0/
  integer :: chem_pot_unit

!..cubic hermite polynomial statement functions
!..psi0 & derivatives
  xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.d0
  xdpsi0(z) = z * (6.0d0*z - 6.0d0)

!..psi1 & derivatives
  xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
  xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


!..bicubic hermite polynomial statement function
  h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) =                    &
       fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t                     &
       + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt                    &
       + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t                     &
       + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt                    &
       + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t                     &
       + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt                    &
       + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t                     &
       + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt

!.. mecc = me*c^2 in g*cm^2/s^2
!.. 1.602176462d-6 ..conversion factor between MeV and g*cm^2/s^2
  mecc    = unit%mass_e*1.602176462d-6
  positron_start = 0.02d0

!..read the table and construct the deltas only once
  if (ifirst .eq. 0) then
     ifirst = 1
     chem_pot_unit= open_infile(chem_pot_file)

     tlo   = 3.0d0
     thi   = 13.0d0
     tstp  = (thi - tlo)/float(jmax-1)
     tstpi = 1.0d0/tstp
     dlo   = -12.0d0
     dhi   = 15.0d0
     dstp  = (dhi - dlo)/float(imax-1)
     dstpi = 1.0d0/dstp

     do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
           dsav = dlo + (i-1)*dstp
           d(i) = 10.0d0**(dsav)
           read(chem_pot_unit,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
     enddo

!..construct the temperature and density deltas and their inverses
     do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dti         = 1.0d0/dth
        dt_sav(j)   = dth
        dti_sav(j)  = dti
     end do
     do i=1,imax-1
        dd          = d(i+1) - d(i)
        ddi         = 1.0d0/dd
        dd_sav(i)   = dd
        ddi_sav(i)  = ddi
     enddo

!..close the file
     close(chem_pot_unit)
  end if

!..normal execution starts here
!..enter the table with ye*den

  din = ye*den


!..bomb proof the input
  if (temp .gt. t(jmax)) then
     call raise_exception("Temperature ("//num_to_str(temp)//&
          ") was off the grid (hotter than "//num_to_str(t(jmax))//"K) "//NEW_LINE("A")//&
          "when trying to interpolate electron and positron chemical potentials."//NEW_LINE("A")//&
          "Check your conditions or consider to switch off theoretical weak rates (iwformat = 0).", &
          "chempot",&
          130003)
     ! write(6,'(1x,5(a,1pe11.3))') 'temp=',temp,' t(jmax)=',t(jmax)
     ! write(6,*) 'temp too hot, off grid, returning'
     ! return
  end if
  if (temp .lt. t(1)) then
     call raise_exception("Temperature ("//num_to_str(temp)//&
          ") was off the grid (colder than "//num_to_str(t(1))//"K) "//NEW_LINE("A")//&
          "when trying to interpolate electron and positron chemical potentials."//NEW_LINE("A")//&
          "Check your conditions or consider to switch off theoretical weak rates (iwformat = 0).", &
          "chempot",&
          130004)
     ! write(6,'(1x,5(a,1pe11.3))') 'temp=',temp,' t(1)=',t(1)
     ! write(6,*) 'temp too cold, off grid, returning'
     ! return
  end if
  if (din  .gt. d(imax)) then
     call raise_exception("Density * Ye ("//num_to_str(din)//&
          ") was off the grid (larger than "//num_to_str(d(imax))//") "//NEW_LINE("A")//&
          "when trying to interpolate electron and positron chemical potentials."//NEW_LINE("A")//&
          "Density: "//num_to_str(den)//" g/ccm, Ye: "//num_to_str(ye)//NEW_LINE("A")//&
          "Check your conditions or consider to switch off theoretical weak rates (iwformat = 0).", &
          "chempot",&
          130005)
     ! write(6,'(1x,5(a,1pe11.3))') 'den*ye=',din,' d(imax)=',d(imax)
     ! write(6,*) 'ye*den too big, off grid, returning'
     ! return
  end if
  if (din  .lt. d(1)) then
     call raise_exception("Density * Ye ("//num_to_str(din)//&
          ") was off the grid (smaller than "//num_to_str(d(1))//") "//NEW_LINE("A")//&
          "when trying to interpolate electron and positron chemical potentials."//NEW_LINE("A")//&
          "Density: "//num_to_str(den)//" g/ccm, Ye: "//num_to_str(ye)//NEW_LINE("A")//&
          "Check your conditions or consider to switch off theoretical weak rates (iwformat = 0).", &
          "chempot",&
          130006)
     ! write(6,'(1x,5(a,1pe11.3))') 'ye*den=',din,' d(1)=',d(1)
     ! write(6,*) 'ye*den too small, off grid, returning'
     ! return
  end if

!..initialize
  deni    = 1.0d0/den
  tempi   = 1.0d0/temp
  kt      = unit%kerg * temp
  ktinv   = 1.0d0/kt
  beta    = kt/mecc
  dbetadt = unit%kerg/mecc

!..hash locate this temperature and density
  jat = int((log10(temp) - tlo)*tstpi) + 1
  jat = max(1,min(jat,jmax-1))
  iat = int((log10(din) - dlo)*dstpi) + 1
  iat = max(1,min(iat,imax-1))

!..look in the electron chemical potential table only once
  fi(1)  = ef(iat,jat)
  fi(2)  = ef(iat+1,jat)
  fi(3)  = ef(iat,jat+1)
  fi(4)  = ef(iat+1,jat+1)
  fi(5)  = eft(iat,jat)
  fi(6)  = eft(iat+1,jat)
  fi(7)  = eft(iat,jat+1)
  fi(8)  = eft(iat+1,jat+1)
  fi(9)  = efd(iat,jat)
  fi(10) = efd(iat+1,jat)
  fi(11) = efd(iat,jat+1)
  fi(12) = efd(iat+1,jat+1)
  fi(13) = efdt(iat,jat)
  fi(14) = efdt(iat+1,jat)
  fi(15) = efdt(iat,jat+1)
  fi(16) = efdt(iat+1,jat+1)

!..get the interpolation weight functions

  xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
  xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
  mxt = 1.0d0 - xt
  mxd = 1.0d0 - xd

  si0t   =  xpsi0(xt)
  si1t   =  xpsi1(xt)*dt_sav(jat)

  si0mt  =  xpsi0(mxt)
  si1mt  =  -xpsi1(mxt)*dt_sav(jat)

  si0d   =  xpsi0(xd)
  si1d   =  xpsi1(xd)*dd_sav(iat)

  si0md  =  xpsi0(mxd)
  si1md  =  -xpsi1(mxd)*dd_sav(iat)

!..derivatives of weight functions
  dsi0t  = xdpsi0(xt)*dti_sav(jat)
  dsi1t  = xdpsi1(xt)

  dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
  dsi1mt = xdpsi1(mxt)

  dsi0d  = xdpsi0(xd)*ddi_sav(iat)
  dsi1d  = xdpsi1(xd)

  dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
  dsi1md = xdpsi1(mxd)

!..electron chemical potential
  etaele  = h3(iat,jat,                                            &
       si0t,   si1t,   si0mt,   si1mt,                     &
       si0d,   si1d,   si0md,   si1md)

!..derivative with respect to density
  x       = h3(iat,jat,                                            &
       si0t,   si1t,   si0mt,   si1mt,                     &
       dsi0d,  dsi1d,  dsi0md,  dsi1md)
  detadd  = ye * x

!..derivative with respect to temperature
  detadt  = h3(iat,jat,                                            &
       dsi0t,  dsi1t,  dsi0mt,  dsi1mt,                     &
       si0d,   si1d,   si0md,   si1md)

!..positron chemical potential
  etapos  = 0.0d0
  detapdd = 0.0d0
  detapdt = 0.0d0

  if (beta .gt. positron_start) then
     etapos = -etaele - 2.0d0/beta
     detapdd = -detadd
     detapdt = -detadt + 2.0d0/beta**2 * dbetadt
  end if

  return

end subroutine chempot
