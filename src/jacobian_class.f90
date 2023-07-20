!> @file jacobian_class.f90
!!
!! The error file code for this file is ***W28***.
!! @brief Module jacobian_class for calculating and solving the linearized equations
!!

!> Contains subroutines to assemble the system Jacobian and changes in the abundances
!!
!! \b Edited:
!!     - 23.01.21, MR - fixed overflows and turned off fission in NSE
#include "macros.h"
module jacobian_class

  use parameter_class, only: nuflag,use_timmes_mue, freeze_rate_temp
  use pardiso_class,   only: vals,jind
  use benam_class
  use nucstuff_class
  use effphase_class
  use nuflux_class
  use error_msg_class, only: raise_exception,num_to_str,int_to_str
  use global_class,    only: nreac
  implicit none

  ! helper variables
  real(r_kind),parameter,private :: rate_min_cutoff  = 1d-50 !< Minimum cutoff for the reaction rates
  real(r_kind),parameter,private :: rate_max_cutoff  = 1d99  !< Maximum cutoff for the reaction rates
  logical,private                :: freeze_rate_indicator = .False. !< Indicator for debug statement
  real(r_kind),private           :: old_time=-99, old_temp=-99, old_dens=-99, old_ye=-99, old_rad=-99
  !
  ! Public and private fields and methods of the module (all routines are publis)
  !
  public:: &
     jacobi_init, abchange
  private:: &
     calculate_reaction_rate

contains



!> Subroutine to calculate the reaction rates
!!
!! The calculation depends on the source flag of the reaction rate.
!! While reaclib rates are done by a simple fit formula, tabulated rates
!! have to be interpolated, theoretical weak rates have to interpolate
!! on a temperature-density grid, etc.
!!
!! @see reaclib_rate_module::calculate_reacl_rate, tabulated_rate_module::calculate_tab_rate,
!!      tw_rate_module::calculate_twr_rate
!!
!! \b Edited:
!!   - 26.07.22, MR: Created this subroutine from code snippets
!!                   that were flying around in Jacobi_init and abchange.
!! .
subroutine calculate_reaction_rate(time, temp, rho, Ye, rkm, rrate_array, idx, rat_calc)
  use global_class,           only: reactionrate_type, nreac
  use tw_rate_module,         only: mue, calculate_twr_rate, weak
  use tabulated_rate_module,  only: calculate_tab_rate, tabulated, tabulated_index
  use nuflux_class,           only: calculate_nu_rate
  use nucstuff_class,         only: pf, inter_partf, calc_t9_pow
  use reaclib_rate_module,    only: calculate_reacl_rate
  use parameter_class,        only: unit, iwformat, use_timmes_mue, nuflag, heating_mode
  use detailed_balance,       only: calculate_detailed_balance
  use screening_module,       only: screen, hv, iscreen
  use error_msg_class,        only: raise_exception,int_to_str
  implicit none
  ! Declare the pass
  type(reactionrate_type),dimension(nreac),intent(inout) :: rrate_array!< Rate array
  real(r_kind),intent(in)                                :: time       !< Time [s]
  real(r_kind),intent(in)                                :: temp       !< Temperature
  integer,intent(in)                                     :: idx        !< Index of rate in array
  real(r_kind),intent(in)                                :: rho        !< Density
  real(r_kind),intent(in)                                :: Ye         !< Electron fraction
  real(r_kind),intent(in)                                :: rkm        !< Radius in km
  real(r_kind),intent(out)                               :: rat_calc   !< Calculated reaction rate
  ! Internal variables
  character*4                         :: src          !< Source flag of reaction
  logical                             :: is_same_step !< Flag that determines if it is the same step
                                                      !  or if the conditions changed
  real(r_kind)                        :: eta_e,eta_p
  real(r_kind)                        :: elnd         !< Electron density
  integer                             :: reac_origin  !< Origin of the reaclib rate (tabulated, reaclib, ...)


  ! Check if indices and other things have to be updated
  if ((old_time .eq. time) .and. (old_temp .eq. temp) .and. &
      (old_dens .eq. rho)  .and. (old_ye .eq. Ye) .and.&
      (old_rad .eq. rkm)) then
    is_same_step = .True.
  else
    old_time = time
    old_temp = temp
    old_dens = rho
    old_ye   = Ye
    old_rad  = rkm
    is_same_step = .False.
  end if

  ! NOTE, MR: Probably one could implement a flag for an hydrostatic run here
  !           to avoid calculating the rates again for this case.

  ! Update general things as indices and the chemical potential
  if (is_same_step .eqv. .False.)  then
     ! t9_pow calculates powers of t9 and rho needed to evaluate reaction rates
     call calc_t9_pow(temp)

     ! inter_partf interpolates the partition functions for a given temperature
     call inter_partf (temp, pf)

     ! screen (by U.Frischknecht) calculates the screening coefficients hv
     if(iscreen) call screen(temp,rho,nreac,Ye,0) ! Calculates coefficient hv

     ! tabulated_index determines the boundaries for the temp interpolation
     if (tabulated) call tabulated_index (temp)

     ! Get the chemical potential for theoretical weak rates
     if (weak .and. (iwformat .eq. 2) .and. use_timmes_mue) then
        call chempot(1.d9*temp,rho,ye,eta_e,eta_p)
        mue = eta_e * temp * 1.d9 * unit%k_MeV
        mue = mue + unit%mass_e
     end if

     ! Interpolate neutrino quantities at given time
     if (nuflag.ge.1) then
        call nuflux(time, rkm)   ! returns fluxnu(4)
        call nutemp(time)        ! returns tempnu(4)
        call nucs()
     end if
  end if

  ! Get the source flag and origin of the rate, depending on the flag the calculation
  ! of the rate differs
  src = rrate_array(idx)%source
  reac_origin = rrate_array(idx)%reac_src

  ! Calculate of the reaction rates
  ! theoretical weak rates.
  ! Here is the place to implement new type of reactions
  if (reac_origin .eq. rrs_twr) then
    call calculate_twr_rate(rrate_array(idx), temp, rho, Ye, rat_calc,.True.)
  ! tabulated rates
  else if (reac_origin .eq. rrs_tabl) then
    call calculate_tab_rate(rrate_array(idx), temp, rat_calc)
  ! Neutrino reactions
  else if (reac_origin .eq. rrs_nu) then
    call calculate_nu_rate(rrate_array(idx),rat_calc)
  ! Detailed balance rate
  else if (reac_origin .eq. rrs_detb) then
    call calculate_detailed_balance(temp,rrate_array,rrate_array(idx),rat_calc)
  ! Fission rate (treat like reaclib)
  else if (reac_origin .eq. rrs_fiss) then
    call calculate_reacl_rate(rrate_array(idx),rat_calc)
  ! Beta delayed neutron emission format (treat like reaclib)
  else if (reac_origin .eq. rrs_wext) then
    call calculate_reacl_rate(rrate_array(idx),rat_calc)
  ! Additional alpha decays (treat like reaclib)
  else if (reac_origin .eq. rrs_aext) then
    call calculate_reacl_rate(rrate_array(idx),rat_calc)
  ! Reaclib rate
  else if (reac_origin .eq. rrs_reacl) then
    call calculate_reacl_rate(rrate_array(idx),rat_calc)
  ! Not known -> throw exception
  else
    call raise_exception("Reaction origin not known, got '"//int_to_str(reac_origin)//&
                         "'."//NEW_LINE("A")//"Reaction source is: "//src,&
                         "calculate_reaction_rate", 280011)
  end if

  ! Cut the rate here and at the end already to avoid overflows
  if(rat_calc.gt.rate_max_cutoff) rat_calc=rate_max_cutoff

  ! Save the result somewhere
  rrate_array(idx)%cached_p = rat_calc

  ! Apply screening corrections
  if(iscreen .eqv. .true.) rat_calc=rat_calc*dexp(hv(idx))

  ! Take care of double counting factor
  rat_calc = rat_calc * rrate_array(idx)%one_over_n_fac

  ! For electron captures the amount of electrons, i.e., the electron density
  ! has to be known.
  elnd = rho * Ye
  if (trim(adjustl(src)) .eq. 'ec') rat_calc = rat_calc * elnd

  ! Apply partition functions for inverse rates and multiply density
  select case (rrate_array(idx)%group)
  case(1:3,11)
    if (rrate_array(idx)%is_reverse) then
      rat_calc = rat_calc * (pf(rrate_array(idx)%parts(2))*pf(rrate_array(idx)%parts(3))* &
                             pf(rrate_array(idx)%parts(4))*pf(rrate_array(idx)%parts(5)))/ &
                             pf(rrate_array(idx)%parts(1))
    end if
  case(4:7)
    rat_calc = rat_calc * rho
    if (rrate_array(idx)%is_reverse) then
      rat_calc = rat_calc * (pf(rrate_array(idx)%parts(3))*pf(rrate_array(idx)%parts(4))  * &
                             pf(rrate_array(idx)%parts(5))*pf(rrate_array(idx)%parts(6))) / &
                            (pf(rrate_array(idx)%parts(1))*pf(rrate_array(idx)%parts(2)))
    end if
  case(8:9)
    rat_calc = rat_calc * rho**2
    if (rrate_array(idx)%is_reverse) then
      rat_calc = rat_calc * (pf(rrate_array(idx)%parts(4))*pf(rrate_array(idx)%parts(5))    &
                          / (pf(rrate_array(idx)%parts(1))*pf(rrate_array(idx)%parts(2))    &
                            *pf(rrate_array(idx)%parts(3))))
    end if
  case(10)
    rat_calc = rat_calc * rho**3
    if (rrate_array(idx)%is_reverse) then
      rat_calc = rat_calc * (pf(rrate_array(idx)%parts(5))*pf(rrate_array(idx)%parts(6)))   &
                          / (pf(rrate_array(idx)%parts(1))*pf(rrate_array(idx)%parts(2))    &
                            *pf(rrate_array(idx)%parts(3))*pf(rrate_array(idx)%parts(4)))
    end if
  case default
    ! Throw exception
    call raise_exception("Got unknown reaclib chapter ("//   &
                         int_to_str(rrate_array(idx)%group)//"). ", &
                         "calculate_reaction_rate",280009)
  end select

  ! Finally check overflows or too low rates
  ! to prevent bad rate fits producing NaN errors
  if(rat_calc.gt.rate_max_cutoff) rat_calc=rate_max_cutoff

  ! Cut the rate if it is too low
  if(rat_calc.lt.rate_min_cutoff) rat_calc=0.d0

  ! Cache the rate
  rrate_array(idx)%cached = rat_calc

end subroutine calculate_reaction_rate



!>
!! This subroutine calculates the change in abundances (dYdt) due to
!! reaction rates at given temperature (temp) and density (rho).
!!
!! The underlying differential equation is given by:
!!
!! \f[
!! \mathrm{dYdt}_i = \underbrace{\sum _j N_j ^i \lambda _j Y_j}_\text{Decays and photodisintegrations} +
!! \underbrace{\sum _{j,k} \frac{N^i _{j,k}}{1+\delta _{jk}}\rho N_{\mathrm{A}} \langle \sigma v
!! \rangle _{j,k}Y_jY_k}_\text{two-body reactions}+\underbrace{\sum _{j,k,l}\frac{N^i_{j,k,l}}
!! {1+\Delta_{jkl}}\rho ^2N_{\mathrm{A}}^2\langle \sigma v\rangle _{j,k,l} Y_j Y_k Y_l}_\text{three-body reactions}
!! \f]
!!
!! where the amount of synthesized and destroyed nuclei N_j is stored in
!! global_class::rrate\%ch_amount, the reaction rates come from various sources,
!! but mostly from reaclib (see \ref network_init_module::network_init, nucstuff_class::calc_t9_pow).
!! The individual terms in the equation above are treated with the help of the reaclib chapters,
!! the first term is given by chapter 1-3, the second one by 4-7, the third one 8-9
!! (see also [Reaclib chapters](https://reaclib.jinaweb.org/help.php?topic=reaclib_format&intCurrentNum=0)).
!! <table>
!! <caption id="multi_row">Reaclib chapters</caption>
!! <tr><th> Chapter   <th> Equation
!! <tr><td> 1         <td> \f$ e_1             \rightarrow e_2                 \f$
!! <tr><td> 2         <td> \f$ e_1             \rightarrow e_2 + e_3           \f$
!! <tr><td> 3         <td> \f$ e_1             \rightarrow e_2 + e_3 + e_4     \f$
!! <tr><td> 4         <td> \f$ e_1 + e_2       \rightarrow e_3                 \f$
!! <tr><td> 5         <td> \f$ e_1 + e_2       \rightarrow e_3 + e_4           \f$
!! <tr><td> 6         <td> \f$ e_1 + e_2       \rightarrow e_3 + e_4 + e_5     \f$
!! <tr><td> 7         <td> \f$ e_1 + e_2       \rightarrow e_3 + e_4 + e_5 +e_6\f$
!! <tr><td> 8         <td> \f$ e_1 + e_2 + e_3 \rightarrow e_4 (+ e_5 + e_6)   \f$
!! </table>
!!
!! @note When single_zone_vars::evolution_mode is equal to EM_NSE,
!! this subroutine will ignore strong reactions and only
!! calculate the contribution of weak reactions.
!!
!! @note This routine can be used to calculate the jacobian of the original Reaclib format
!!       without chapter 9-11 as well as the modern format with chapter 9-11.
!!
!! @see jacobi_init,
!!      [Hix & Thielemann 1999](https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..321H/abstract),
!!      [Winteler 2013](https://edoc.unibas.ch/29895/),
!!      [Lippuner & Roberts 2017](https://ui.adsabs.harvard.edu/abs/2017ApJS..233...18L/abstract)
!!
!! \b Edited:
!!     - 12.01.14
!!     - 26.07.22, MR: Cleaned up the subroutine
!!     - 06.02.23, MR: Implemented parameter freeze_rate_temp.
!! .
subroutine abchange (time, itemp, rho, Ye, rkm, Y, dYdt, evolution_mode)
   use global_class,        only: nreac, rrate, net_names
   use fission_rate_module, only: fissrate, nfiss, fiss_neglect
   implicit none

   real(r_kind),intent(in)                 :: time     !< time [s]
   real(r_kind),intent(in)                 :: itemp    !< temperature in units of 10**9 K
   real(r_kind),intent(in)                 :: rho      !< baryon density [g/cm^3]
   real(r_kind),intent(in)                 :: Ye       !< electron fraction
   real(r_kind),intent(in)                 :: rkm      !< current radius im km
   real(r_kind),dimension(:),intent(in)    :: Y        !< isotope abundances
   real(r_kind),dimension(:),intent(inout) :: dYdt     !< change in abundances due to reaction rates
   integer,intent(in)                      :: evolution_mode !< state of the network
   !
   real(r_kind)                        :: rat
   integer                             :: i, j
   real(r_kind)                        :: infty
   real(r_kind)                        :: temp

   infty = HUGE(infty)

   dYdt = 0.d0

   ! Freeze the rates at the minimum valid temperature, 1d-2GK
   if (itemp.le.freeze_rate_temp) then
      temp = freeze_rate_temp
      if (.not. freeze_rate_indicator) then
         print *,"Freezing rates at "//num_to_str(freeze_rate_temp)//" GK"
         freeze_rate_indicator = .True.
      endif
   else
      temp = itemp
   end if

   ! Loop through all reactions
   outer: do i=1,nreac

      ! Only consider weak reactions in NSE
      if ((evolution_mode.eq.EM_NSE).and.(rrate(i)%is_weak.eqv..false.)) cycle outer

      ! Calculate the reaction rate
      call calculate_reaction_rate(time, temp, rho, Ye, rkm, rrate, i, rat)

      ! Check if one can ignore the rate
      if (skip_rate(rat, rrate(i),Y)) cycle outer

      ! Ignore negligible rates
      if (rat .lt. rate_min_cutoff) cycle outer

      ! Don't do fission
      if ((rrate(i)%reac_type .eq. rrt_nf) .or. &
          (rrate(i)%reac_type .eq. rrt_bf) .or. &
          (rrate(i)%reac_type .eq. rrt_sf)) then
          cycle outer
      end if

      ! Create DGL
      select case (rrate(i)%group)
      case(1:3,11)
         do j = 1, 5
            if (rrate(i)%parts(j) .eq. 0) exit
            dYdt(rrate(i)%parts(j)) = dYdt(rrate(i)%parts(j)) + rat * &
                 rrate(i)%ch_amount(j) * Y(rrate(i)%parts(1))
         end do
      case(4:7)
         do j = 1, 6
            if (rrate(i)%parts(j) .eq. 0) exit
            dYdt(rrate(i)%parts(j)) = dYdt(rrate(i)%parts(j)) + rat * &
                 rrate(i)%ch_amount(j) * Y(rrate(i)%parts(1)) *       &
                 Y(rrate(i)%parts(2))
         end do
      case(8:9)
         do j = 1, 6
            if (rrate(i)%parts(j) .eq. 0) exit
            dYdt(rrate(i)%parts(j)) = dYdt(rrate(i)%parts(j)) + rat * &
                 rrate(i)%ch_amount(j) * Y(rrate(i)%parts(1)) *       &
                 Y(rrate(i)%parts(2)) * Y(rrate(i)%parts(3))
         end do
      case(10)
         do j = 1, 6
            if (rrate(i)%parts(j) .eq. 0) exit
            dYdt(rrate(i)%parts(j)) = dYdt(rrate(i)%parts(j)) + rat * &
                 rrate(i)%ch_amount(j) * Y(rrate(i)%parts(1)) *       &
                 Y(rrate(i)%parts(2)) * Y(rrate(i)%parts(3))  *       &
                 Y(rrate(i)%parts(4))
         end do
      case default
        ! Throw exception
        call raise_exception("Got unknown reaclib chapter ("//&
                             int_to_str(rrate(i)%group)//"). "&
                             ,"abchange",280009)
      end select

      ! Check for infinity rates
      if ((rat .ge. infty) .or. (rat .ne. rat)) then
         call raise_exception("The rate: "//NEW_LINE("A")//&
                              trim(adjustl(reaction_string(rrate(i))))//&
                              NEW_LINE("A")//"was infinity or NaN.",&
                              "abchange",280003)
      end if
   end do outer

   ! MR: Changed that fission reactions are not taken into account in NSE
   if ((fissflag.ne.0) .and. (evolution_mode.ne.EM_NSE)) then
      do i=1,nfiss                    ! loop over fission reactions
         rat = 0.d0
         do j=1,9
            rat = rat +fissrate(i)%param(j)*t9_pow(j)
         end do

         ! Check overflows
         if (rat .lt. dlog(rate_max_cutoff)) then
            rat = dexp(rat)
         else
            rat = rate_max_cutoff
         end if

         if(rat.lt.rate_min_cutoff) rat=0.d0

         if ((fissrate(i)%reac_type .eq. rrt_bf) .or. (fissrate(i)%reac_type .eq. rrt_sf)) then
               ! Check rate for infinity or nan
               if ((rat .ge. infty) .or. (rat .ne. rat)) then
                  call raise_exception("Fission rate (spont. and b-delayed) of "//&
                                        trim(adjustl(net_names(fissrate(i)%fissparts(1))))//&
                                        " was infinity or NaN.","abchange",280004)
               end if
               if (Y(fissrate(i)%fissparts(1)) .eq. 0) then
                  do j = 1,fissrate(i)%dimens
                      ! Calculate the equation
                      dYdt(fissrate(i)%fissparts(j)) =     &
                          dYdt(fissrate(i)%fissparts(j)) + &
                          rat * fissrate(i)%ch_amount(j) * Y(fissrate(i)%fissparts(1))
                  end do
               end if
         else if (fissrate(i)%reac_type .eq. rrt_nf) then                          ! n-induced fission
            rat = rat * rho
            ! Avoid multiplying when it will be zero anyways, this also prevents things to get done
            ! when a rate is possible infinity
            ! Check rate for infinity or nan
            if ((rat .ge. infty) .or. (rat .ne. rat)) then
                call raise_exception("Neutron induced fission rate of "//&
                                    trim(adjustl(net_names(fissrate(i)%fissparts(2))))//&
                                    " was infinity or NaN.","abchange",280005)
            end if
            if (Y(fissrate(i)%fissparts(1))*Y(fissrate(i)%fissparts(2)) .eq. 0) then
                do j = 1,fissrate(i)%dimens
                    ! Calculate the equation
                    dYdt(fissrate(i)%fissparts(j)) =    &
                        dYdt(fissrate(i)%fissparts(j)) + &
                        rat * fissrate(i)%ch_amount(j) * Y(fissrate(i)%fissparts(1)) * Y(fissrate(i)%fissparts(2))
                end do
            end if
         end if

      end do
   end if


   return

end subroutine abchange



!> Function to decide whether a rate can contribute or not
!!
!! This function returns true if the rate should not
!! be considered in the differential equation.
!! This can happen by a small value of the rate or
!! by zero abundance of the reactants.
!!
!! @author M. Reichert
!! @date 14.02.23
function skip_rate(rat, rrate_in,Y) result(res)
    use nucstuff_class, only: get_nr_reactants
    implicit none
    real(r_kind),intent(in)                       :: rat
    type(reactionrate_type),intent(in)            :: rrate_in
    real(r_kind), dimension(net_size), intent(in) :: Y
    logical :: res
    integer :: i
    integer :: zero_count

    res = .false.
      ! Ignore negligible rates
    if (rat .lt. rate_min_cutoff) res = .True.

    ! Don't do fission
    if ((rrate_in%reac_type .eq. rrt_nf) .or. &
        (rrate_in%reac_type .eq. rrt_bf) .or. &
        (rrate_in%reac_type .eq. rrt_sf)) then
        res = .True.
    end if

    ! Check if the reaction can be skipped
    zero_count = 0
    do i=1,get_nr_reactants(rrate_in%group)
        if (Y(rrate_in%parts(i)) .le. 0) zero_count = zero_count+1
        if (zero_count .gt. 1) exit
    end do
    if (zero_count .gt. 1) res = .True.


end function skip_rate


!>
!! This subroutine calculates the entries in the jacobian and the
!! right-hand side of the numerical integration scheme (rhs).
!!
!! The jacobian is calculated by
!! \f[ J = \frac{\mathrm{dYdt}}{dY}.  \f]
!! The values are directly stored into a sparse format (see \ref pardiso_class::vals).
!!
!! @note When single_zone_vars::evolution_mode is equal to EM_NSE,
!! this subroutine will ignore strong reactions and only
!! calculate the contribution of weak reactions.
!!
!! @note This routine can be used to calculate the jacobian of the original Reaclib format
!!       without chapter 9-11 as well as the modern format with chapter 9-11.
!!
!! @see abchange,
!!      [Hix & Thielemann 1999](https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..321H/abstract),
!!      [Winteler 2013](https://edoc.unibas.ch/29895/),
!!      [Lippuner & Roberts 2017](https://ui.adsabs.harvard.edu/abs/2017ApJS..233...18L/abstract)
!!
!! \b Edited:
!!     - 26.07.22, MR: Cleaned up the subroutine
!!     - 06.02.23, MR: Implemented parameter freeze_rate_temp.
!! .
subroutine jacobi_init (time, itemp, rho, rkm, Y, Y_p, dYdt, rhs, h, evolution_mode)
   use global_class,        only: net_size, isotope, net_names, ineu
   use global_class,        only: reactionrate_type, nreac, rrate
   use fission_rate_module, only: fissrate, nfiss, fiss_neglect
   use pardiso_class,       only: dia, pt_e, pt_b, rows
   use gear_module,         only: get_l1, get_predictor_Y, get_predictor_dYdt
   use nuclear_heating,     only: calculate_qdot, reset_qdot
   use screening_module,    only: iscreen, hv
   use parameter_class,     only: heating_mode

   real(r_kind),intent(in)                 :: time     !< time [s]
   real(r_kind),intent(in)                 :: itemp    !< initial temperature in units of 10**9 K
   real(r_kind),intent(in)                 :: rho      !< baryon density [g/cm^3]
   real(r_kind),intent(in)                 :: rkm      !< current radius im km
   real(r_kind),dimension(:),intent(in)    :: Y        !< current isotope abundances
   real(r_kind),dimension(:),intent(in)    :: Y_p      !< isotope abundances for the previous step
   real(r_kind),dimension(:),intent(inout) :: dYdt     !< change in abundances due to reaction rates
   real(r_kind),dimension(:),intent(inout) :: rhs      !< right-hand side of the num. integration scheme
   real(r_kind),intent(in)                 :: h        !< timestep size
   integer,intent(in)                      :: evolution_mode !< Evolution mode of the network
   ! Internal variables
   integer                                  :: i, j    !< Loop variables
   real(r_kind)                             :: d       !< Diagonals
   real(r_kind)                             :: temp    !< Temperature [GK]
   real(r_kind)                             :: rat     !< Reaction rate
   real(r_kind)                             :: Ye      !< electron fraction
   type(reactionrate_type)                  :: rr_tmp  !< Reaction rate instance
   character*1,dimension(6)                 :: descra  !< Descriptor for sparse matrix multiplication
   real(r_kind)                             :: infty   !< Infinity to check for errors
   integer                                  :: fl_c    !< fission loop count
   real(r_kind)                             :: l_1


   INFO_ENTRY("jacobi_init")

   ! Initialize variable for infinity checks
   infty = HUGE(infty)

   ! Descriptor for sparse matrix multiplication
   descra = (/"G"," "," ","F"," "," "/)
   ! The diagonal value
   d = 1.d0/h

!-----compute electron fraction and molar density
   Ye= sum(isotope(1:net_size)%p_nr*Y(1:net_size))

   ! Reset the qdot for the nuclear heating
   call reset_qdot(time)

   ! Check Ye
   if (Ye.ne.Ye) then
      call raise_exception('Ye is NaN','jacobi_init',280006)
   elseif (Ye.lt.0) then
      call raise_exception('Ye is negative ('//trim(adjustl(num_to_str(ye)))//").",&
                           'jacobi_init',280007)
   end if

   ! Check Rho * ye
   if(dlog10(rho*ye).ne.dlog10(rho*ye)) then
      ! Raise an exception
      call raise_exception('log10(rho*ye) is NaN.'//NEW_LINE("A")//&
                            "Density is : "//trim(adjustl(num_to_str(rho)))//" [g/ccm]."//NEW_LINE("A")//&
                            "Ye is      : "//trim(adjustl(num_to_str(ye))),&
                            'jacobi_init',280009)
   end if

   ! Freeze the rates at the minimum valid temperature, 1d-2GK
   if (itemp.le.freeze_rate_temp) then
      temp = freeze_rate_temp
      if (.not. freeze_rate_indicator) then
         print *,"Freezing rates at "//num_to_str(freeze_rate_temp)//" GK"
         freeze_rate_indicator = .True.
      endif
   else
      temp = itemp
   end if

   ! Initialize dydt and vals
   vals = 0.d0
   dYdt = 0.d0

   ! switch between different solvers
   if(solver==0) then ! implicit Euler
      do i=1,net_size
         vals(dia(i)) = 1.d0/h
      end do
   elseif(solver==1) then ! Gear's method
      l_1 = get_l1()
      do i=1,net_size
         vals(dia(i)) = l_1/h
      end do
   endif

   ! Loop through reactions
   outer: do i=1,nreac
      rr_tmp = rrate(i)

      ! Consider only weak reactions in NSE
      if ((evolution_mode.eq.EM_NSE).and.(rr_tmp%is_weak.eqv..false.)) cycle outer

      ! Calculate the reaction rate
      call calculate_reaction_rate(time, temp, rho, Ye, rkm, rrate, i, rat)

      ! Check if one can ignore the rate
      if (skip_rate(rat, rrate(i),Y)) cycle outer

      ! Create the jacobian and the DGLs
      select case (rr_tmp%group)
      ! Chapter 1-3: one educt
      case(1:3,11)
         do j = 1, 5
            if (rr_tmp%parts(j) .eq. 0) exit
            dYdt(rr_tmp%parts(j)) = dYdt(rr_tmp%parts(j)) + rat *       &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(1))

            vals(rr_tmp%cscf_ind(1,j)) = vals(rr_tmp%cscf_ind(1,j)) - rat *      &
                 rr_tmp%ch_amount(j)
         end do
         ! Chapter 1-4: two educts
      case(4:7)
         do j = 1, 6
            if (rr_tmp%parts(j) .eq. 0) exit
            dYdt(rr_tmp%parts(j)) = dYdt(rr_tmp%parts(j)) + rat *       &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(1)) *              &
                 Y(rr_tmp%parts(2))

            vals(rr_tmp%cscf_ind(1,j)) = vals(rr_tmp%cscf_ind(1,j)) - rat *      &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(2))
            vals(rr_tmp%cscf_ind(2,j)) = vals(rr_tmp%cscf_ind(2,j)) - rat *      &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(1))
         end do
      ! Chapter 8-9: three educts
    case(8:9)
         do j = 1, 6
            if (rr_tmp%parts(j) .eq. 0) exit
            dYdt(rr_tmp%parts(j)) = dYdt(rr_tmp%parts(j)) + rat *       &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(1)) *              &
                 Y(rr_tmp%parts(2)) * Y(rr_tmp%parts(3))

            vals(rr_tmp%cscf_ind(1,j)) = vals(rr_tmp%cscf_ind(1,j)) - rat *      &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(2)) *        &
                 Y(rr_tmp%parts(3))
            vals(rr_tmp%cscf_ind(2,j)) = vals(rr_tmp%cscf_ind(2,j)) - rat *      &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(1)) *        &
                 Y(rr_tmp%parts(3))
            vals(rr_tmp%cscf_ind(3,j)) = vals(rr_tmp%cscf_ind(3,j)) - rat *      &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(1)) *        &
                 Y(rr_tmp%parts(2))
         end do
      ! Chapter 10: four educts
    case(10)
         do j = 1, 6
            if (rr_tmp%parts(j) .eq. 0) exit
            dYdt(rr_tmp%parts(j)) = dYdt(rr_tmp%parts(j)) + rat *       &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(1)) *              &
                 Y(rr_tmp%parts(2)) * Y(rr_tmp%parts(3)) * Y(rr_tmp%parts(4))

            vals(rr_tmp%cscf_ind(1,j)) = vals(rr_tmp%cscf_ind(1,j)) - rat *      &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(2)) *        &
                 Y(rr_tmp%parts(3)) * Y(rr_tmp%parts(4))
            vals(rr_tmp%cscf_ind(2,j)) = vals(rr_tmp%cscf_ind(2,j)) - rat *      &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(1)) *        &
                 Y(rr_tmp%parts(3)) * Y(rr_tmp%parts(4))
            vals(rr_tmp%cscf_ind(3,j)) = vals(rr_tmp%cscf_ind(3,j)) - rat *      &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(1)) *        &
                 Y(rr_tmp%parts(2)) * Y(rr_tmp%parts(4))
            vals(rr_tmp%cscf_ind(4,j)) = vals(rr_tmp%cscf_ind(4,j)) - rat *      &
                 rr_tmp%ch_amount(j) * Y(rr_tmp%parts(1)) *        &
                 Y(rr_tmp%parts(2)) * Y(rr_tmp%parts(3))
         end do
      case default
        ! Throw exception
        call raise_exception("Got unknown reaclib chapter ("//&
                             int_to_str(rr_tmp%group)//"). "  &
                             ,"jacobi_init",280009)
      end select

      ! Cache the rate with the densities included
      rrate(i)%cached = rat

      ! Check for infinity rates
      ! In jacobi init it may not occur as there are several statements setting
      ! the rate to 1e99 if it is greater than this
      if ((rat .gt. infty) .or. (rat .ne. rat)) then
        if (iscreen) then
         call raise_exception("The rate: "//NEW_LINE("A")//&
                              trim(adjustl(reaction_string(rrate(i))))//&
                              NEW_LINE("A")//"was infinity or NaN."//new_line("A")//&
                              "Logarithm of screening factor was "//num_to_str(hv(i)),&
                              "jacobi_init",280003)
        else
         call raise_exception("The rate: "//NEW_LINE("A")//&
                              trim(adjustl(reaction_string(rrate(i))))//&
                              NEW_LINE("A")//"was infinity or NaN.",&
                              "jacobi_init",280003)
        end if
      end if

      ! Call the function to calculate the qdot for heating
      call calculate_qdot(rrate(i),Y,h)

   end do outer

   ! MR: Changed that fission reactions are not taken into account in NSE
   if ((fissflag.ne.0) .and. (evolution_mode.ne.EM_NSE)) then
      fissloop: do i=1,nfiss ! loop over fission reactions to implement correct equations for the fragments
         rat = 0.d0
         do j=1,9
            rat = rat +fissrate(i)%param(j)*t9_pow(j)
         end do

         ! Check overflows
         if (rat .lt. dlog(rate_max_cutoff)) then
            rat = dexp(rat)
         else
            ! Say something for higher verbose levels
            if (VERBOSE_LEVEL .ge. 2) then
               if (fissrate(i)%fissparts(1) .ne. ineu) then
                    write(*,*)"Fission rate is overflowing. Fissioning nucleus is "//&
                              trim(adjustl(isotope(fissrate(i)%fissparts(1))%name))
               else
                    write(*,*)"Fission rate is overflowing. Fissioning nucleus is "//&
                                trim(adjustl(isotope(fissrate(i)%fissparts(2))%name))
               end if
            end if
            rat = rate_max_cutoff
         end if

         if(rat.ge.rate_min_cutoff) then

            if ((fissrate(i)%reac_type .eq. rrt_bf) .or. (fissrate(i)%reac_type .eq. rrt_sf)) then ! Spontaneous and beta delayed fission
                ! Check if the parent nucleus is present
                ! If not, the reaction is not possible and may not be evaluated. Note that this
                ! can cause problems as even for a zero value the jacobian will have an entry.
                ! Do this also only for implicit euler as gear calls the jacobian with a delta
                ! instead of abundances.
                if ((Y(fissrate(i)%fissparts(1)) .le. 0) .and. (solver .eq. 0)) then
                    ! Since we ordered them, the most important nuclei are at the begining of the
                    ! fissparts. We can therefore approximate the jacobian by the first few entries.
                    fl_c = min(fissrate(i)%dimens, fiss_neglect)
                else
                    ! If the parent nucleus is present, we can calculate the full jacobian
                    fl_c = fissrate(i)%dimens
                end if

                do j = 1, fl_c
                    dYdt(fissrate(i)%fissparts(j)) =     &
                        dYdt(fissrate(i)%fissparts(j)) + &
                        rat * fissrate(i)%ch_amount(j) * Y(fissrate(i)%fissparts(1))
                    vals(fissrate(i)%cscf_ind(1,j)) = vals(fissrate(i)%cscf_ind(1,j)) - rat * fissrate(i)%ch_amount(j)
                end do

            else if (fissrate(i)%reac_type .eq. rrt_nf) then ! n-induced fission
                rat = rat * rho
                ! Also check here if the parent nucleus is present
                ! If not, the reaction is not possible and may not be evaluated. Note that this
                ! can cause problems as even for a zero value the jacobian will have an entry.
                ! Therefore we set at least the first few entries.
                if ((Y(fissrate(i)%fissparts(1)) * Y(fissrate(i)%fissparts(2)) .le. 0) .and. (solver .eq. 0)) then
                    ! Since we ordered them, the most important nuclei are at the begining of the
                    ! fissparts. We can therefore approximate the jacobian by the first few entries.
                    fl_c = min(fissrate(i)%dimens, fiss_neglect)
                else
                    ! If the parent nucleus is present, we can calculate the full jacobian
                    fl_c = fissrate(i)%dimens
                end if

                do j = 1,fl_c
                    dYdt(fissrate(i)%fissparts(j)) =    &
                        dYdt(fissrate(i)%fissparts(j)) + &
                        rat * fissrate(i)%ch_amount(j) * Y(fissrate(i)%fissparts(1)) * Y(fissrate(i)%fissparts(2))
                    vals(fissrate(i)%cscf_ind(1,j)) = vals(fissrate(i)%cscf_ind(1,j)) - rat * fissrate(i)%ch_amount(j) * Y(fissrate(i)%fissparts(2))
                    vals(fissrate(i)%cscf_ind(2,j)) = vals(fissrate(i)%cscf_ind(2,j)) - rat * fissrate(i)%ch_amount(j) * Y(fissrate(i)%fissparts(1))
                end do

            end if

         end if

         if(rat.lt.rate_min_cutoff) rat=0.d0
         fissrate(i)%cached = rat
      end do fissloop
   end if


   ! switch between different solvers
   if(solver==0) then ! implicit Euler
     ! additional check for input into pardiso
      if (VERBOSE_LEVEL .ge. 2) then
         do i=1,net_size
            if (dYdt(i).ne.dYdt(i)) then
               call raise_exception("dYdt of "//trim(adjustl(net_names(i)))//" is NaN.",&
                                    "jacobi_init",280010)
            end if
         end do
      end if
      rhs = dYdt + d*(Y_p - Y)
      ! computes matrix-vector product for a sparse matrix in the CSC format:
      call mkl_dcscmv('T',net_size,net_size,1.d0,descra,vals,rows,pt_b,  &
           pt_e,Y,1.d0,rhs)
   elseif(solver==1) then ! Gear's method
      rhs = -(Y - get_predictor_Y())*(l_1/h) + (dYdt - get_predictor_dYdt()/h)
   endif
   dYdt = dYdt + d*(Y_p-Y)

   INFO_EXIT("jacobi_init")

end subroutine jacobi_init


end module jacobian_class
