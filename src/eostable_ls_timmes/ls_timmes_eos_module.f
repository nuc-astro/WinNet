#include "../macros.h"
      module ls_timmes_eos_module

      use units_module

      implicit none

!.....LS & Timmes switch parameters.....................................
      real(r_kind),parameter :: Tref=5.00e9 ![K] temperature LS-Timmes switch
                                    ! 5.80e9 K approx. 0.5 MeV
      real(r_kind),parameter :: dref=1.00e7 ![g/cm3] density LS-Timmes switch

!.....EoS input variable parameters.....................................
      integer,parameter :: ink=0
      integer,parameter :: ine=1
      integer,parameter :: ins=2

!.....nr & bisection parameters.........................................
      real(r_kind),parameter :: T1=1.00e2  !lower temperature bracket [K]
      real(r_kind),parameter :: T2=1.00e12 !upper temperature bracket [K]
                                   ! 1.17e12 K approx. 100 MeV

      type timmes_eos_state
        integer  :: ef      !eos flag
        real(r_kind)     :: T       !Temperature                [K]
        real(r_kind)     :: F       !Helmholtz free energy      [erg/g]
        real(r_kind)     :: S       !Entropy                    [kB/baryon]
        real(r_kind)     :: P       !Pressure                   [erg/cm3]
        real(r_kind)     :: d2FdTdd !d^2F/dTdrho                [erg*cm3/(g2*K)]
        real(r_kind)     :: d2FdT2  !d^2F/dT^2                  [erg/(g*K2)]
        real(r_kind)     :: cs      !Sound speed                [cm/s]
        real(r_kind)     :: e       !Internal energy            [erg/g]
        real(r_kind)     :: abar    !Average mass number        [-]
        real(r_kind)     :: abarxx  !Average mass number        [-]
        real(r_kind)     :: zbar    !Average charge number      [-]
        real(r_kind)     :: xa      !Mass fraction of alpha particles   [-]
        real(r_kind)     :: cv      !Specific heat capacity at constant volume
        real(r_kind)     :: xp      !Mass fraction of exterior protons  [-]
        real(r_kind)     :: xn      !Mass fraction of exterior neutrons [-]
        real(r_kind)     :: muhat   !Neutron chemical potential minus
                            ! proton chemical potential         [MeV]
        real(r_kind)     :: mun     !Neutron chemical potential         [MeV]
        real(r_kind)     :: mue     !Electron chemical potential        [MeV]
      end type

      contains

      subroutine timmes_eos(var,vin,d,Ye,state,status)

      implicit none

      integer,intent(in)                      :: var
      real(r_kind),intent(in)                         :: vin,d,Ye
      type(timmes_eos_state),intent(inout)    :: state
      integer,intent(out)                     :: status

!-----------------------------------------------------------------------
!
!     LS & Timmes EoS routine for temperature, int. energy or entropy
!
!     Input:  var       ... selected input variable (see definitions
!                                                  up in module)
!          vin=T,e,s    ... Temperature [K], internal energy [erg/g],
!                           entropy [kB/baryon]
!             d         ... Density            [g/cm3]
!             Ye        ... Electron abundance [-]
!
!     Output: state  ... see timmes_eos_state type declaration
!             status ... 0 for zero problems and 1 for trouble
!
!-----------------------------------------------------------------------

      real(r_kind) :: T

!.....call relevant routine.............................................
      if ( var .eq. ink ) then
        call timmes_eos_interface(vin,d,Ye,state,status)
      elseif ( var .eq. ine ) then
        !find corresponding temperature
        call timmes_eos_e_bisec(vin,d,Ye,state%abar,T,status)
        if ( status .eq. 0 ) call timmes_eos_interface(T,d,Ye,state,
     1                                                 status)
      elseif ( var .eq. ins ) then
        !find corresponding temperature
        call timmes_eos_s_bisec(vin,d,Ye,state%abar,T,status)
        if ( status .eq. 0 ) call timmes_eos_interface(T,d,Ye,state,
     1                                                 status)
      else
        if (VERBOSE_LEVEL .gt. 1) then
            write(6,*) 'timmes_eos_module: Check input variable '
     1                ,'setting...'
        endif
        status = 1
      endif

      return

      end subroutine timmes_eos

      subroutine timmes_eos_interface(T,d,Ye,state,status)

      use timmes_eos_module

      implicit none

      real(r_kind),intent(in)                         :: T,d,Ye
      type(timmes_eos_state),intent(inout) :: state
      integer,intent(out)                     :: status

!-----------------------------------------------------------------------
!
!     Timmes EoS interface
!
!     NOTE: state%abar is used as an input!!!
!
!     Input:  T          ... Temperature         [K]
!             d          ... Density             [g/cm3]
!             Ye         ... Electron abundance  [-]
!             state%abar ... average mass number [-]
!
!     Output: state  ... see timmes_eos_state type declaration
!             status ... 0 for zero problems and 1 for trouble
!
!-----------------------------------------------------------------------

      include 'vector_eos.dek'

!.....set input for timmes eos..........................................
      temp_row(1) = T
      den_row(1)  = d
      abar_row(1) = state%abar
      zbar_row(1) = state%abar*Ye

!.....call timmes eos...................................................
      !here the pipeline is only 1 element long
      jlo_eos = 1
      jhi_eos = 1

      call eosfxt

!.....set output variables..............................................
      state%T = T                               ![K]
      state%e = etot_row(1)                     ![erg/g]
      state%S = stot_row(1)/(units%kB/units%mb) ![kB/baryon]
      state%F = etot_row(1) - T*stot_row(1)     ![erg/g]
      state%p = ptot_row(1)                     ![erg/cm3]
      state%d2FdTdd = dpt_row(1)/(d*d)          ![erg*cm3/(g2*K)]
      state%d2FdT2 = -dst_row(1)                ![erg/(g*K2)]
      state%cs = cs_row(1)                      ![cm/s]
      state%zbar = state%abar*Ye                ![-]
      state%xa = 0.                             ![-]
      state%xp = 0.                             ![-]
      state%xn = 0.                             ![-]
      state%cv = cv_row(1)/units%Na/units%MeV   ![Mev/(baryon*K)]
      state%muhat = 0.                          ![MeV]
      state%mun = 0.                            ![MeV]
      state%mue = etaele_row(1)                 ![MeV]
!      state%mue = 0.                           ![MeV]

!debug:
!        write(6,*) 'pion = ',pion_row(1)
!        write(6,*) 'prad = ',prad_row(1)
!        write(6,*) 'pele = ',pele_row(1) + ppos_row(1)
!        write(6,*) 'pcou = ',pcou_row(1)
!        write(6,*) 'ptot = ',ptot_row(1)

!.....set status according EoS..........................................
      if ( .not.eosfail ) then
        status = 0
      else
        status = 1
      endif

      end subroutine timmes_eos_interface

      subroutine timmes_eos_e_bisec(e,d,Ye,abar,T,status)

      implicit none

      real(r_kind),intent(in)        :: e,d,Ye,abar
      real(r_kind),intent(out)       :: T
      integer,intent(out)    :: status

!-----------------------------------------------------------------------
!
!     Find Temperature for given internal energy of LS-Timmes EoS
!       with bisection algorithm ... see Numerical Recipes
!
!     Input:  e          ... Internal energy     [erg/g]
!             d          ... Density             [g/cm3]
!             Ye         ... Electron abundance  [-]
!
!     Output: T      ... Temperature [K]
!             status ... 0 for zero problems and 1 for trouble
!
!-----------------------------------------------------------------------

      integer,parameter      :: jmax=100 !max number of bisec iterations
      real(r_kind),parameter            :: precision=5.e-4 !eps*(T1+T2)/2
      integer                   :: j
      type(timmes_eos_state) :: state
      real(r_kind)                      :: Tmid,dT,f,fmid

!.....start bisection...................................................
      state%abar = abar
      call timmes_eos_interface(T1,d,Ye,state,status)
      if ( status .eq. 0 ) then
        f = e - state%e
      else
        if (VERBOSE_LEVEL .gt. 1) then
            write(6,*) 'ls_timmes_eos_module: Check bracket T1 =',T1
     1                 ,'in eos_find_T.'
        end if
        status = 1
        return
      endif
      call timmes_eos_interface(T2,d,Ye,state,status)
      if ( status .eq. 0 ) then
        fmid = e - state%e
      else
        if (VERBOSE_LEVEL .gt. 1) then
            write(6,*) 'ls_timmes_eos_module: Check bracket T2 =',T2
     1                ,'in eos_find_T.'
        end if
        status = 1
        return
      endif
      if( f*fmid .ge. 0. ) then
        if (VERBOSE_LEVEL .gt. 1) then
            write(6,*) 'ls_timmes_eos_module: Root not bracketed... '
     1                ,'try other values of T1,T2 ...'
        end if
        status = 1
        return
      endif
      if ( f .lt. 0. ) then
        T = T1
        dT = T2 - T1
      else
        T = T2
        dT = T1 - T2
      endif
      do j=1,jmax
        dT = 0.5*dT
        Tmid = T + dT
        call timmes_eos_interface(Tmid,d,Ye,state,status)
        if ( status .eq. 0 ) then
          fmid = e - state%e
        else
          if (VERBOSE_LEVEL .gt. 1) then
            write(6,*) 'ls_timmes_eos_module: Error while bisecting...'
          end if
          status = 1
          return
        endif
        if ( fmid .le. 0. ) T = Tmid
        if ( (abs(dT) .lt. precision) .or. (fmid .eq. 0.) ) then
          status = 0
          return
        endif
      enddo

!.....too many iterations...............................................
      status = 1
      if (VERBOSE_LEVEL .gt. 1) then
        write(6,*) 'ls_timmes_eos_module: Too many iterations...'
      end if

      end subroutine timmes_eos_e_bisec

      subroutine timmes_eos_e_nr(e,d,Ye,abar,T,status)

      implicit none

      real(r_kind),intent(in)     :: e,d,Ye,abar
      real(r_kind),intent(inout)  :: T
      integer,intent(out) :: status

!-----------------------------------------------------------------------
!
!     Find Temperature for given internal energy of LS-Timmes EoS
!       with Newton-Raphson algorithm ... see Numerical Recipes
!
!     Input:  e          ... Internal energy     [erg/g]
!             d          ... Density             [g/cm3]
!             Ye         ... Electron abundance  [-]
!
!     Output: T      ... Temperature [K]
!             status ... 0 for zero problems and 1 for trouble
!
!-----------------------------------------------------------------------

      integer,parameter         :: jmax=100 !max number of nr iterations
      real(r_kind),parameter            :: precision=5.e-4 !eps*(T1+T2)/2
      integer                   :: j
      type(timmes_eos_state) :: state
      real(r_kind)                      :: f,df


!.....check the eventually given initial guess for T....................
      state%abar = abar
      call timmes_eos_interface(T,d,Ye,state,status)
      if ( status .eq. 1 ) then
        T = T2
        call timmes_eos_interface(T,d,Ye,state,status)
      endif
      f  = state%e - e
      df = -T*state%d2FdT2

!.....nr iterations.....................................................
      do j=1,jmax
        T = T - f/df
        if ( T .le. 0. ) T = T2
        call timmes_eos_interface(T,d,Ye,state,status)
        if ( status .ne. 0 ) then
          if (VERBOSE_LEVEL .gt. 1) then
              write(6,*) 'ls_timmes_eos_module: NR iteration failed '
     1                  ,'with T =',T
          end if
          status = 1
          return
        endif
        f  = state%e - e
        df = -T*state%d2FdT2
        if ( (f .eq. 0.) .or. (abs(f/df) .le. precision) ) then
          status = 0
          return
        endif
      enddo

!.....too many iterations...............................................
      status = 1
      if (VERBOSE_LEVEL .gt. 1) then
          write(6,*) 'ls_timmes_eos_module: Too many iterations...'
      end if
      end subroutine timmes_eos_e_nr

      subroutine timmes_eos_s_bisec(s,d,Ye,abar,T,status)

      implicit none

      real(r_kind),intent(in)        :: s,d,Ye,abar
      real(r_kind),intent(out)       :: T
      integer,intent(out)    :: status

!-----------------------------------------------------------------------
!
!     Find Temperature for given entropy of LS-Timmes EoS
!       with bisection algorithm ... see Numerical Recipes
!
!     Input:  s          ... Entropy     [kB/baryon]
!             d          ... Density             [g/cm3]
!             Ye         ... Electron abundance  [-]
!
!     Output: T      ... Temperature [K]
!             status ... 0 for zero problems and 1 for trouble
!
!-----------------------------------------------------------------------

      integer,parameter      :: jmax=100 !max number of bisec iterations
      real(r_kind),parameter            :: precision=5.e-4 !eps*(T1+T2)/2
      integer                   :: j
      type(timmes_eos_state) :: state
      real(r_kind)                      :: Tmid,dT,f,fmid

!.....start bisection...................................................
      state%abar = abar
      call timmes_eos_interface(T1,d,Ye,state,status)
      if ( status .eq. 0 ) then
        f = s - state%s
      else
        if (VERBOSE_LEVEL .gt. 1) then
            write(6,*) 'ls_timmes_eos_module: Check bracket T1 =',T1
     1                ,'in eos_find_T.'
        end if
        status = 1
        return
      endif
      call timmes_eos_interface(T2,d,Ye,state,status)
      if ( status .eq. 0 ) then
        fmid = s - state%s
      else
        if (VERBOSE_LEVEL .gt. 1) then
            write(6,*) 'ls_timmes_eos_module: Check bracket T2 =',T2
     1                ,'in eos_find_T.'
        end if
        status = 1
        return
      endif
      if( f*fmid .ge. 0. ) then
        if (VERBOSE_LEVEL .gt. 1) then
            write(6,*) 'ls_timmes_eos_module: Root not bracketed... '
     1                ,'try other values of T1,T2 ...'
        end if
        status = 1
        return
      endif
      if ( f .lt. 0. ) then
        T = T1
        dT = T2 - T1
      else
        T = T2
        dT = T1 - T2
      endif
      do j=1,jmax
        dT = 0.5*dT
        Tmid = T + dT
        call timmes_eos_interface(Tmid,d,Ye,state,status)
        if ( status .eq. 0 ) then
          fmid = s - state%s
        else
          if (VERBOSE_LEVEL .gt. 1) then
            write(6,*) 'ls_timmes_eos_module: Error while bisecting...'
          end if
          status = 1
          return
        endif
        if ( fmid .le. 0. ) T = Tmid
        if ( (abs(dT) .lt. precision) .or. (fmid .eq. 0.) ) then
          status = 0
          return
        endif
      enddo

!.....too many iterations...............................................
      status = 1
      if (VERBOSE_LEVEL .gt. 1) then
          write(6,*) 'ls_timmes_eos_module: Too many iterations...'
      end if
      end subroutine timmes_eos_s_bisec

      subroutine timmes_eos_s_nr(s,d,Ye,abar,T,status)

      implicit none

      real(r_kind),intent(in)     :: s,d,Ye,abar
      real(r_kind),intent(inout)  :: T
      integer,intent(out) :: status

!-----------------------------------------------------------------------
!
!     Find Temperature for given entropy of LS-Timmes EoS
!       with Newton-Raphson algorithm ... see Numerical Recipes
!
!     Input:  s          ... Entropy             [kB/baryon]
!             d          ... Density             [g/cm3]
!             Ye         ... Electron abundance  [-]
!
!     Output: T      ... Temperature [K]
!             status ... 0 for zero problems and 1 for trouble
!
!-----------------------------------------------------------------------

      integer,parameter         :: jmax=100 !max number of nr iterations
      real(r_kind),parameter            :: precision=5.e-4 !eps*(T1+T2)/2
      integer                   :: j
      type(timmes_eos_state) :: state
      real(r_kind)                      :: f,df


!.....check the eventually given initial guess for T....................
      state%abar = abar
      call timmes_eos_interface(T,d,Ye,state,status)
      if ( status .eq. 1 ) then
        T = T2
        call timmes_eos_interface(T,d,Ye,state,status)
      endif
      f  = (state%s - s)*(units%kB/units%mb)
      df = -state%d2FdT2

!.....nr iterations.....................................................
      do j=1,jmax
        T = T - f/df
        if ( T .le. 0. ) T = T2
        call timmes_eos_interface(T,d,Ye,state,status)
        if ( status .ne. 0 ) then
          if (VERBOSE_LEVEL .gt. 1) then
              write(6,*) 'ls_timmes_eos_module: NR iteration failed '
     1                  ,'with T =',T
          end if
          status = 1
          return
        endif
        f  = (state%s - s)*(units%kB/units%mb)
        df = -state%d2FdT2
        if ( (f .eq. 0.) .or. (abs(f/df) .le. precision) ) then
          status = 0
          return
        endif
      enddo

!.....too many iterations...............................................
      status = 1
      if (VERBOSE_LEVEL .gt. 1) then
          write(6,*) 'ls_timmes_eos_module: Too many iterations...'
      end if
      end subroutine timmes_eos_s_nr

      end module ls_timmes_eos_module
