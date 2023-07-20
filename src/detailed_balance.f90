!> @file detailed_balance.f90
!!
!! The error file code for this file is ***W26***.
!!
!! Contains the module \ref detailed_balance


!> @brief Module that deals with inverse reaction rates.
!!
!! This module is only used when \ref parameter_class::use_detailed_balance is
!! enabled. The detailed balance is calculated as follows:
!!
!! \f[
!! \sigma_\mathrm{inv} = N_A^{-n}
!!                       \delta
!!                       \left(\frac{\prod_i A_i}{\prod_j A_j}\right)^{1.5}
!!                       \left(\frac{\prod_i g_i}{\prod_j g_j}\right)
!!                       \left( \frac{m_u k_B T}{2\pi \hbar^2}\right)^{1.5 \cdot n}
!!                       e^{-Q/(k_B T)}
!!                       \sigma_\mathrm{forw},
!! \f]
!! where \f$n\f$ is the number of reactants minus the number of products,
!! \f$ N_A \f$ is Avogadros number, \f$A_i\f$ is the mass number of reactants of the forward reaction,
!! \f$A_j\f$ the number of products of the forward reaction, \f$g_i\f$ the spin factor \f$ 2J_i +1\f$
!! of reactants (j, for products), Q the q-value of the forward reaction, T the temperature
!! and \f$ \sigma_\mathrm{forw} \f$ the forward reaction.
!!
!! The necessity of calculating detailed balance within the network is given by the
!! tabulated rates that may not easily obey the detailed balance principle otherwise.
!! Furthermore, there exist inconsistencies between the reaclib Q-value and the
!! Q-value calculated by the winvn (version downloaded 24.07.22).
!! The Q-value given in the reaclib are mostly used for calculating the inverse rates.
!! However there seem to exist also some other inconsistencies.
!! An example of an inconsistent q-value in reaclib is given by the reaction:
!! n + k55 -> k56
!! other inconsistencies are included in, e.g.,
!! n + ge89 -> ge90
!! or p + co52 -> ni53 which deviates by a (almost) constant factor of 4 indicating
!! a different spin that was used in reaclib to calculate the inverse reaction.
!! @image html k56.png "Reaction calculated with detailed balance and taken from reaclib (25.07.22)" height=280
!! @image html ni53.png "Reaction calculated with detailed balance and taken from reaclib (25.07.22)" height=280
!! This problem was also topic of
!! [Lippuner & Roberts 2017](https://ui.adsabs.harvard.edu/abs/2017ApJS..233...18L/abstract).
!!
!! The essence of a discussion with H. Schatz is that reaclib reactions can be based on other mass models/data
!! than the winvn table. Every reaclib reaction can have its own set of mass tables. It might be better to
!! calculate inverse reaction yourself in the network as you will be more consistent, especially close to
!! equilibrium values such as NSE. However, this goes to the cost of an inconsistency between the inverse
!! and forward rate as the forward rate involves other mass models for the calculation. For (n,gamma) reactions
!! that are based on Hauser-Feshbach, Hendrik suggested that one could use a scaling relation of the forward rate
!! based on the Q-value differences.
!!
!!
!! @author M. Reichert
!! @date 22.07.22
#include "macros.h"
module detailed_balance
    use global_class,    only: reactionrate_type
    use error_msg_class, only: raise_exception, int_to_str
    use parameter_class, only: use_detailed_balance
    implicit none

    integer,private :: ninv     !< Number of inverse rates
    integer,private :: ninv_old !< Number of old inverse rates

    ! Helper variables for ignoring specific sources
    character(len=4),allocatable,dimension(:),private :: src_ignore
    integer,private                                   :: src_ignore_length
    character(len=4),allocatable,dimension(:),private :: src_q_reacl
    integer,private                                   :: src_q_reacl_length
    character(len=4),allocatable,dimension(:),private :: src_q_winvn
    integer,private                                   :: src_q_winvn_length

    !
    ! Public and private fields and methods of the module
    !
    public:: &
        calculate_detailed_balance, merge_inverse_rates, init_inverse_rates
    private:: &
        create_inverse_rate, write_reac_verbose_out, write_reac_rate


  contains



    !> Initialize everything for the inverse rates
    !!
    !! This subroutine analyzes the parameters to determine which rates have to
    !! be replaced and which Q-value has to be used.
    !!
    !! @author M. Reichert
    !! @date 04.08.22
    subroutine init_inverse_rates()
      use parameter_class, only: detailed_balance_src_ignore,detailed_balance_src_q_reac,&
                                 detailed_balance_src_q_winvn
      use nucstuff_class,  only: analyze_src_string
      implicit none

      ! Dont do anything if not enabled
      if (.not. use_detailed_balance) return

      INFO_ENTRY("init_inverse_rates")

      ! Get the ignore source string as array
      call analyze_src_string(detailed_balance_src_ignore,src_ignore,src_ignore_length)
      ! Source strings for which the Q-value of reaclib should be used
      call analyze_src_string(detailed_balance_src_q_reac,src_q_reacl,src_q_reacl_length)
      ! Source strings for which the Q-value of winvn should be used
      call analyze_src_string(detailed_balance_src_q_winvn,src_q_winvn,src_q_winvn_length)

      ! Give some verbose output
      if (VERBOSE_LEVEL .ge. 2) then
        write(*,*) "Detailed balance: ignoring source    : ",src_ignore
        write(*,*) "Detailed balance: use Q-value reaclib: ",src_q_reacl
        write(*,*) "Detailed balance: use Q-value winvn  : ",src_q_winvn
      end if

      INFO_EXIT("init_inverse_rates")

    end subroutine init_inverse_rates



    !> Delete and create new inverse rates
    !!
    !! This routine first identifies included reverse rates and deletes them.
    !! Afterwards, new reverse rates are created and attached to the rate array at the end.
    !! The position of the rates at the end of the array is important as they will use the
    !! previously calculated forward rates. Weak reactions, fission, and chapter 7
    !! reactions do not have any inverse reaction. This routine also stores
    !! the index of the forward reaction in %param(1).
    !!
    !! @author  M. Reichert
    !! @date 22.07.22
    subroutine merge_inverse_rates(rrate_array,rrate_length)
      use benam_class,     only: getcoefficients
      use parameter_class, only: use_detailed_balance_q_reac
      implicit none
      ! Declare the pass
      type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate_array  !< Large rate array, containing all reactions
      integer,intent(inout)                                          :: rrate_length !< length of rrate_array
      ! Internal varialbes
      integer                                          :: i               !< Loop variable
      type(reactionrate_type)                          :: rr_tmp          !< Temporary reaction rate
      integer                                          :: newlength       !< New rate array length
      integer                                          :: count           !< Helper variable
      type(reactionrate_type),dimension(:),allocatable :: new_rrate_array !< New rate array
      integer                                          :: rstat           !< Allocation status
      integer,dimension(11)                            :: ignore_inv      !< Ignore these chapters for inverse rates

      ! Only do something when parameter is enabled
      if (.not. use_detailed_balance) return

      INFO_ENTRY("merge_inverse_rates")

      ! For testing:
      ! Chapters that appear in the following array as number
      ! will use the reverse rates of the reaclib. E.g.,
      ! the following will ignore chapter 2 and 4:
      ! ignore_inv = (/-1,2,-1,4,-1,-1,-1,-1,-1,-1,-1/)
      ! With only "-1" included, nothing will be ignored
      ! and all reaclib reverse rates will be replaced
      ! by detailed balance ones.
      ignore_inv = (/-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/)

      ! Count the inverse rates
      ninv = 0
      ninv_old = 0
      ! Go through the rates and check if one has to add an inverse one
      do i=1,rrate_length
        ! For shorter writing
        rr_tmp = rrate_array(i)

        if (VERBOSE_LEVEL .gt. 2) then
          ! This happens in the current reaclib...
          if (rr_tmp%is_weak .and. rr_tmp%is_reverse) then
             write(*,*) "Detected rate that is weak and reverse!"
          end if
        end if

        ! Skip weak rates, neutrino rates, fission
        if (rr_tmp%is_weak) cycle
        if (rr_tmp%reac_src .eq. rrs_nu) cycle
        if (rr_tmp%reac_src .eq. rrs_fiss) cycle
        ! Ignore rates that have to be specified in the ignore list
        if (any(src_ignore .eq. adjustr(rr_tmp%source))) cycle

        ! In case chapters should get ignored
        if (any(ignore_inv==rr_tmp%group)) cycle

        ! Count old inverse rates
        if (rr_tmp%is_reverse) then
          ninv_old = ninv_old +1
          cycle
        end if

        ! Every other rate should have an inverse
        ninv = ninv+1
      end do

      ! Allocate the new array
      newlength = rrate_length-ninv_old+ninv
      allocate(new_rrate_array(newlength), stat = rstat)
      if (rstat /= 0) call raise_exception('Allocation of "new_rrate_array" failed.',&
                                           "merge_inverse_rates", 260001)

      ! Write all non reverse rates in new array
      count = 1
      do i=1,rrate_length
        ! For shorter writing
        rr_tmp = rrate_array(i)

        ! Skip reverse rates (however, not the weak ones)
        ! Ignore some chapters if you want..
        ! Ignore the list of ignore sources..
        if (rr_tmp%is_reverse .and. (.not. rr_tmp%is_weak) .and. &
            (.not. (any(ignore_inv==rr_tmp%group)))        .and. &
            (.not. (any(src_ignore == adjustr(rr_tmp%source))))) then
          cycle
        end if

        ! Add all non reverse rates here
        new_rrate_array(count) = rr_tmp

        ! next rate
        count = count+1
      end do

      ! Now create the reverse rates
      do i=1,newlength

        ! For shorter writing
        rr_tmp = new_rrate_array(i)

        ! ignore some chapters in case that it is enabled
        if (any(ignore_inv==rr_tmp%group) .or. &
            any(src_ignore == adjustr(rr_tmp%source))) then
             cycle
        end if

        ! Skip reverse rates
        if (rr_tmp%is_reverse) then
          cycle
        end if

        ! Skip weak rates, neutrino rates, fission
        if (rr_tmp%is_weak) cycle
        if (rr_tmp%reac_src .eq. rrs_nu) cycle
        if (rr_tmp%reac_src .eq. rrs_fiss) cycle

        ! Add reverse rate
        call create_inverse_rate(rr_tmp,new_rrate_array(count))
        ! index of the forward rate
        new_rrate_array(count)%param(1)=i

        ! Set the Q-value of forward and inverse reactions based on the
        ! parameter "use_detailed_balance_q_reac". If this is true,
        ! the Q-value from the reaclib will be used, otherwise
        ! it will be calculated from the winvn
        if ((use_detailed_balance_q_reac) .or. (any(src_q_reacl== rr_tmp%source)) .and. &
            (.not. any(src_q_winvn== rr_tmp%source)) ) then
          ! Set the q_value consistent to the reaclib!
          new_rrate_array(count)%q_value = -new_rrate_array(i)%q_value
        elseif ((.not. use_detailed_balance_q_reac) .or. (any(src_q_winvn== rr_tmp%source)) .and. &
                (.not. any(src_q_reacl== rr_tmp%source))) then
          ! Set the q_value consistent to the winvn as they differ for some reactions!
          new_rrate_array(i)%q_value = -new_rrate_array(count)%q_value
        end if

        ! next rate
        count = count+1
      end do

      ! Get coefficients like one_over_nfrac for inverse rates
      if (ninv .ne. 0) then
        call getcoefficients(new_rrate_array(newlength-ninv+1:newlength),ninv)
      end if

      ! Give verbose output
      if (VERBOSE_LEVEL .ge. 2) then
        call write_reac_rate(new_rrate_array,newlength,rrate_array,rrate_length)
      end if

      ! Write into the correct array and resize the old one
      deallocate(rrate_array, stat = rstat)
      if (rstat /= 0) call raise_exception('Deallocation of "rrate_array" failed.',&
                                           "merge_inverse_rates", 260002)
      allocate(rrate_array(newlength), stat = rstat)
      if (rstat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                           "merge_inverse_rates", 260001)

      rrate_array(:) = new_rrate_array(:)
      deallocate(new_rrate_array, stat = rstat)
      if (rstat /= 0) call raise_exception('Deallocation of "new_rrate_array" failed.',&
                                           "merge_inverse_rates", 260002)
      ! set the new length
      rrate_length = newlength

      ! Say how it worked out
      call write_reac_verbose_out()

      INFO_EXIT("merge_inverse_rates")

    end subroutine merge_inverse_rates


    !> Routine to create an inverse rate, based on the forward rate.
    !!
    !! Writes the correct nuclei in the parts array of each reaction.
    !! Furthermore, as the %param() array of the rates is not used, it is redefined to contain:
    !!  - 1: Index to forward rate in rate array
    !!  - 2: #Products - #Reactants.
    !!  - 3: \f$ 1.5 \cdot \log(A_{r1} \cdot ...) / (A_{p1} \cdot ...) \f$,  with A_r being
    !!       the mass number of the reactants and A_p of the products.
    !!  - 4: \f$ 1.5 \cdot \log(M_u/(2.0\cdot \pi \cdot \hbar^2)) \f$.
    !!  - 5: \f$ - \log((2J_{r1}+1) \cdot ...) / ((2J_{p1}+1) \cdot ...) \f$,  with J_r being
    !!       the spin of the reactants and J_p of the products.
    !! . Additionally, the Q-value is calculated from the mass excess.
    !!
    !! @author  M. Reichert
    !! @date 22.07.22
    subroutine create_inverse_rate(forward_rate,inverse_rate)
      use reaclib_rate_module, only: set_reaction_type
      use global_class,        only: isotope
      use parameter_class,     only: unit
      use nucstuff_class,      only: get_nr_products,get_nr_reactants,get_chapter
      implicit none
      type(reactionrate_type),intent(in)  :: forward_rate !< Temporary reaction rate
      type(reactionrate_type),intent(out) :: inverse_rate !< Temporary reaction rate
      ! Internal variables
      integer    :: nr_prod,nr_react !< Number of products/reactants
      integer    :: inv_chapter      !< Inverse reaclib chapter
      integer    :: i,j              !< Loop variable

      ! Initialize the rate
      inverse_rate%param(:)    = 0        ! will only be used for the index of forward rate
      inverse_rate%source      = "detb"   ! detailed balance
      inverse_rate%reac_src    = rrs_detb ! detailed balance
      inverse_rate%is_weak     = .false.  ! not allowed to be a weak rate
      inverse_rate%is_resonant = forward_rate%is_resonant
      inverse_rate%reac_type   = rrt_o    ! The type is set at the end of this subroutine
      inverse_rate%nu_frac     = 0        ! No weak reaction is given a inverse one, so no neutrino should be here

      ! The %parts are used to store basic things in order not to
      ! recalculate everything all over again
      ! -%param(1)  : Index of forward reaction in reaction rate array
      ! -%param(2)  : Difference of amount of reactions - products
      ! -%param(3)  : Prod(A_reac)/Prod(A_prod)
      ! -%param(4)  : Constant for calculation of inverse rate
      ! -%param(5)  : Spin factor

      ! Flag the reverse
      inverse_rate%is_reverse = .True.
      ! Get the correct chapter and put the nuclei on correct place

      ! Get number of products/educts of forward rate (#educt -> #product)
      nr_prod  = get_nr_products(forward_rate%group)
      nr_react = get_nr_reactants(forward_rate%group)

      ! Chapter 8 might be different, check manually once again for it
      if (forward_rate%group .eq. 8) then
        nr_prod = count(forward_rate%parts/=0)-nr_react
      end if

      ! Get the inverse reaclib chapter
      inv_chapter = get_chapter(nr_prod,nr_react)
      inverse_rate%group=inv_chapter

      ! Swap the parts
      inverse_rate%parts(:) = 0
      do i=1,nr_prod
        inverse_rate%parts(i) = forward_rate%parts(i+nr_react)
      end do
      do j=1,nr_react
        inverse_rate%parts(j+nr_prod) = forward_rate%parts(j)
      end do

      ! Calculate the Q-Value
      inverse_rate%q_value = 0
      do i=1,nr_prod
        inverse_rate%q_value = inverse_rate%q_value - &
                               isotope(inverse_rate%parts(i))%mass_exc
      end do
      do j=1,nr_react
        inverse_rate%q_value = inverse_rate%q_value + &
                               isotope(inverse_rate%parts(j+nr_prod))%mass_exc
      end do

      ! Save the difference between #prods and #educts
      inverse_rate%param(2) = DBLE(nr_prod-nr_react)

      ! Factor of masses
      ! Prod(A_reac)/Prod(A_prod)
      inverse_rate%param(3) = 1
      do i=1,nr_prod
        inverse_rate%param(3) = inverse_rate%param(3) * &
                                DBLE(isotope(inverse_rate%parts(i))%mass)
      end do
      do j=1,nr_react
        inverse_rate%param(3) = inverse_rate%param(3) / &
                                DBLE(isotope(inverse_rate%parts(j+nr_prod))%mass)
      end do

      ! Spin factor
      inverse_rate%param(5) = 1
      do i=1,nr_prod
        inverse_rate%param(5) = inverse_rate%param(5) * &
                                (2.0*isotope(inverse_rate%parts(i))%spin+1.0)
      end do
      do j=1,nr_react
        inverse_rate%param(5) = inverse_rate%param(5) / &
                                (2.0*isotope(inverse_rate%parts(j+nr_prod))%spin+1.0)
      end do

      ! Convert to the correct format that is necessary later
      inverse_rate%q_value  = -inverse_rate%q_value
      inverse_rate%param(2) = -inverse_rate%param(2)
      inverse_rate%param(3) = -dlog(inverse_rate%param(3)**(3.0/2.0))
      inverse_rate%param(4) =  (3./2.)*dlog((unit%mass_u/((2.0*unit%pi*(unit%hbc*1d-13)**2.))))
      inverse_rate%param(5) = -dlog(inverse_rate%param(5))

      ! Set the reaction type
      call set_reaction_type(inverse_rate)

    end subroutine create_inverse_rate



    !> Calculate the inverse rate from the forward one via detailed balance.
    !!
    !! This subroutine calculates the inverse rate via detailed balance. Detailed balance
    !! gives a direct relation between forward and inverse reaction.
    !!
    !! @note This subroutine also sets a threshold for inverse rates in case the Q-value
    !!       for a photodissociation is positive. In this case the rates can diverge towards
    !!       lower temperatures. This can mess up convergence of the system due to precision errors.
    !!       For the default FRDM winvn/reaclib q-values this is not a problem, however, when
    !!       using other mass models or when varying the masses it can become a problem.
    !!
    !! @see [Fowler, Caughlan, and Zimmerman (1967)](https://www.annualreviews.org/doi/10.1146/annurev.aa.05.090167.002521)
    !! @author  M. Reichert
    !! @date 22.07.22
    !!
    !! \b Edited:
    !!   - 18.11.22; M.R: implemented more restrictive cutoff for the rates in case of
    !!                    positive q-values
    !! .
    subroutine calculate_detailed_balance(temp,rrate_array,inverse_rate,rate)
      use parameter_class, only: unit
      use global_class,    only: nreac
      implicit none
      ! Declare the pass
      type(reactionrate_type),dimension(nreac),intent(in) :: rrate_array  !< Large rate array, containing all reactions
      type(reactionrate_type),intent(in) :: inverse_rate !< inverse reaction rate to evaluate
      real(r_kind),intent(in)            :: temp         !< Temperature [GK]
      real(r_kind),intent(out)           :: rate         !< Value for the reaction rate
      ! Internal variables
      type(reactionrate_type)   :: forward_rate    !< forward reaction rate
      real(r_kind)              :: kbT             !< kB*T
      real(r_kind)              :: delta           !< Variables with one_over_n_fac
      real(r_kind)              :: rat             !< Storage for log(forward rate)
      real(r_kind),parameter    :: cutoff_cr_h=90  !< Cutoff for the crosssection for low temperatures
      real(r_kind),parameter    :: cutoff_cr_m=70  !< Cutoff for the crosssection for low temperatures
      real(r_kind),parameter    :: cutoff_cr_l=60  !< Cutoff for the crosssection for low temperatures
      real(r_kind),parameter    :: cutoff_T9_h=1d0 !< Temperature below which the rate will be cut
      real(r_kind),parameter    :: cutoff_T9_m=8d-1!< Temperature below which the rate will be cut
      real(r_kind),parameter    :: cutoff_T9_l=5d-1!< Minimum temperature for the cut
      real(r_kind)              :: maxrat          !< helper variable to set a maximum rate

      ! Check that the code really should calculate detailed balance
      if (.not. use_detailed_balance) then
        call raise_exception("Parameter for enabling detailed balance calculation was "//&
                             "turned off, but trying to calculate it."//NEW_LINE("A"),&
                             "calculate_detailed_balance",260005)
      end if

      kbT = (unit%k_mev*1.0d9*temp) ! MeV/GK

      ! Save the forward rate
      forward_rate = rrate_array(INT(inverse_rate%param(1)))

      ! Get the forward rate, as the inverse rates are at the end of the array,
      ! there is no need to recalculate
      rat = dlog(forward_rate%cached_p)

      ! Double counting factors
      delta = dlog(forward_rate%one_over_n_fac/inverse_rate%one_over_n_fac)

      ! Calculation of the rate
      rate = inverse_rate%param(3) +\
             inverse_rate%param(2) * inverse_rate%param(4) +\
             inverse_rate%param(5) +\
             delta +\
             inverse_rate%param(2)*dlog(kbT**(3./2.))+\
             inverse_rate%q_value/(kbT)-\
             inverse_rate%param(2)*dlog(unit%n_a) +\
             rat

      ! I'm restricting the rate here since some rates close to the dripline are
      ! diverging. An example of such an diverging rate is given by
      ! the inverse rate of:
      ! k55 + n -> k56
      ! The problem with this rate is that the q-value according to the winvn is:
      ! mexc(k55) + mexc(n) - mexc(k56) = Q
      ! (-0.270)  + (8.071) - (8.747)   = -0.946 MeV
      ! The Q-value in reaclib is: 1.676 MeV. The reaclib value
      ! is the value that is used for the inverse rate in reaclib.
      ! When using a negative q-value, we get however a rate that is
      ! ~exp(q/kbT) [so positive], so it is diverging towards lower temperatures.
      ! This also makes sense as the system should be unstable when getting energy
      ! from breaking up the nucleus.
      ! As most of the rates are close to zero for low temperatures, a diverging
      ! rate can mess up the DGL system completely leading to a crash by
      ! precision errors.
      if ((temp .lt. cutoff_T9_h ) .and. (inverse_rate%q_value .gt. 0)) then

        ! Limit the rate in a smooth way and in two steps
        if ((temp .gt. cutoff_T9_l) .and. (temp .le. cutoff_T9_m)) then
          ! Linear interpolation
          maxrat = (cutoff_cr_l*(cutoff_T9_m-temp) + cutoff_cr_m*(temp-cutoff_T9_l))/&
                   (cutoff_T9_m-cutoff_T9_l)
        elseif ((temp .gt. cutoff_T9_m) .and. (temp .le. cutoff_T9_h)) then
          ! Linear interpolation
          maxrat = (cutoff_cr_m*(cutoff_T9_h-temp) + cutoff_cr_h*(temp-cutoff_T9_m))/&
                   (cutoff_T9_h-cutoff_T9_m)
        else
          maxrat = cutoff_cr_l
        end if

        ! Limit it
        if ((rate .gt. maxrat)) rate = maxrat

      end if

      ! General cutoff
      if ((rate .gt. cutoff_cr_h)) rate = cutoff_cr_h

      rate = dexp(rate)
    end subroutine calculate_detailed_balance



    !> Verbose routine to output one inverse reaction of chapter 4-9
    !!
    !! This subroutine outputs one forward, one calculated inverse, and
    !! the same inverse as given by the reaclib for chapter 4-9.
    !! All other chapters do not have forward rates so far.
    !! The calculated inverse and reaclib inverse should to large extends
    !! agree. There might be a small deviation due to different q-values
    !! used in the reaclib and in the winvn. This routine is only
    !! called when verbose_level is greater or equal two.
    !!
    !! @note For all chapter 7 reactions, the double counting factor seems
    !!       to be missing in the inverse rate of Reaclib (31.7.2022).
    !!
    !! @author M. Reichert
    !! @date 30.07.22
    subroutine write_reac_rate(rate_array,nrea,rate_array_old,nrea_old)
      use global_class,        only: ihe4, ineu, ipro
      use benam_class,         only: benam
      use nucstuff_class,      only: calc_t9_pow
      use file_handling_class, only: open_outfile, close_io_file
      use reaclib_rate_module, only: calculate_reacl_rate
      implicit none
      ! Declare the pass
      integer,intent(in)                                    :: nrea,nrea_old   !< Length of rate arrays
      type(reactionrate_type),dimension(nrea),intent(inout) :: rate_array      !< Reaction rates with db
      type(reactionrate_type),dimension(nrea),intent(inout) :: rate_array_old  !< Reaction rates from reacl.
      integer                        :: i,j                   !< Loop variables
      integer                        :: idx_ne20,idx_mg24     !< Index of ne20, mg24
      integer                        :: idx_na20,idx_c12      !< Index of na20, c12
      integer                        :: idx_d,idx_li7,idx_he3 !< Index of d, li7, he3
      real(r_kind),dimension(301)    :: temp_grid             !< Temperature grid for output
      real(r_kind),dimension(11,301) :: rat_forward           !< Storage for forward rate
      real(r_kind),dimension(11,301) :: rat_inverse           !< Storage for inverse rate (db)
      real(r_kind),dimension(11,301) :: rat_inverse_reacl     !< Storage for inverse rate (reacl.)
      real(r_kind)                   :: rat_tmp               !< Temporary rate value
      integer                        :: file_id               !< ID of the output file

      ! Logarithmic temperature grid
      temp_grid(1)=-2d0
      do i=2,301
        temp_grid(i) = temp_grid(i-1)+0.01d0
      end do
      temp_grid = 10d0**temp_grid

      ! Initialize the rates.
      rat_forward(:,:)=0d0
      rat_inverse(:,:)=0d0
      rat_inverse_reacl(:,:)=0d0

      ! search for he4 ne20 <-> mg24 as benchmark (chapter 4)
      idx_ne20 = benam(" ne20")
      idx_mg24 = benam(" mg24")
      ! search for  n na20  <->  p ne20 as benchmark (chapter 5)
      idx_na20 = benam(" na20")
      ! search for  d  li7  <->  n  he4  he4 as benchmark (chapter 6)
      idx_d    = benam("    d")
      idx_li7  = benam("  li7")
      ! search for he3  li7 <->  n    p  he4  he4 as benchmark (chapter 7)
      idx_he3    = benam("  he3")
      ! search for  a a a  <->  c12 as benchmark (chapter 8)
      idx_c12  = benam("  c12")
      ! search for n    p    p  <->   p    d  as benchmark (chapter 9)

      ! Calculate the rates at all temperature grid points
      do j=1,301
        call calc_t9_pow(temp_grid(j))
        do i=1,nrea
          ! he4 ne20 -> mg24
          if ((((rate_array(i)%parts(1) .eq. ihe4)      .and. &
                (rate_array(i)%parts(2) .eq. idx_ne20)) .or. &
               ((rate_array(i)%parts(1) .eq. idx_ne20) .and. &
                (rate_array(i)%parts(2) .eq. ihe4)))    .and. &
                (rate_array(i)%parts(3) .eq. idx_mg24) .and. &
                (rate_array(i)%group .eq. 4))    then
                call calculate_reacl_rate(rate_array(i),rat_tmp)
                rat_forward(4,j)       = rat_forward(4,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

          ! he4 he4 he4 -> c12
          if (( (rate_array(i)%parts(1) .eq. ihe4)      .and. &
                (rate_array(i)%parts(2) .eq. ihe4)      .and. &
                (rate_array(i)%parts(3) .eq. ihe4)     .and. &
                (rate_array(i)%parts(4) .eq. idx_c12))    .and. &
                (rate_array(i)%group .eq. 8))    then
                call calculate_reacl_rate(rate_array(i),rat_tmp)
                rat_forward(8,j)       = rat_forward(8,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

          ! n na20  ->  p ne20
          if ((((rate_array(i)%parts(1) .eq. ineu)      .and. &
                (rate_array(i)%parts(2) .eq. idx_na20)) .or. &
               ((rate_array(i)%parts(1) .eq. idx_na20) .and. &
                (rate_array(i)%parts(2) .eq. ineu)))    .and. &
              (((rate_array(i)%parts(3) .eq. ipro) .and. &
                (rate_array(i)%parts(4) .eq. idx_ne20)) .or. &
               ((rate_array(i)%parts(3) .eq. idx_ne20) .and. &
                (rate_array(i)%parts(4) .eq. ipro))) .and. &
                (rate_array(i)%group .eq. 5))    then
                call calculate_reacl_rate(rate_array(i),rat_tmp)
                rat_forward(5,j)       = rat_forward(5,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

          ! d  li7  ->  n   he4  he4
          if ((((rate_array(i)%parts(1) .eq. idx_d)    .and. &
                (rate_array(i)%parts(2) .eq. idx_li7)  .and. &
                (rate_array(i)%parts(3) .eq. ineu))    .and. &
                (rate_array(i)%parts(4) .eq. ihe4)     .and. &
                (rate_array(i)%parts(5) .eq. ihe4))    .and. &
                (rate_array(i)%group .eq. 6))    then
                call calculate_reacl_rate(rate_array(i),rat_tmp)
                rat_forward(6,j)       = rat_forward(6,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

          ! he3  li7 <->  n    p  he4  he4
          if ((((rate_array(i)%parts(1) .eq. idx_he3)    .and. &
                (rate_array(i)%parts(2) .eq. idx_li7)  .and. &
                (rate_array(i)%parts(3) .eq. ineu))    .and. &
                (rate_array(i)%parts(4) .eq. ipro)     .and. &
                (rate_array(i)%parts(5) .eq. ihe4)     .and. &
                (rate_array(i)%parts(6) .eq. ihe4))    .and. &
                (rate_array(i)%group .eq. 7))    then
                call calculate_reacl_rate(rate_array(i),rat_tmp)
                rat_forward(7,j)       = rat_forward(7,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

          ! n    p    p  <->   p    d
          if ((((rate_array(i)%parts(1) .eq. ineu)    .and. &
                (rate_array(i)%parts(2) .eq. ipro)    .and. &
                (rate_array(i)%parts(3) .eq. ipro))   .and. &
                (rate_array(i)%parts(4) .eq. ipro)    .and. &
                (rate_array(i)%parts(5) .eq. idx_d))  .and. &
                (rate_array(i)%group .eq. 9))    then
                call calculate_reacl_rate(rate_array(i),rat_tmp)
                rat_forward(9,j)       = rat_forward(9,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if
        end do

        do i=1,nrea
          ! he4 ne20 <- mg24
          if ((((rate_array(i)%parts(3) .eq. ihe4)       .and. &
                (rate_array(i)%parts(2) .eq. idx_ne20))  .or.  &
               ((rate_array(i)%parts(2) .eq. ihe4)       .and. &
                (rate_array(i)%parts(3) .eq. idx_ne20))) .and. &
                (rate_array(i)%parts(1) .eq. idx_mg24)   .and. &
                (rate_array(i)%group .eq. 2)             .and. &
                (rate_array(i)%source .eq. "detb")) then
                call calculate_detailed_balance(temp_grid(j),rate_array,rate_array(i),rat_tmp)
                rat_inverse(4,j)       = rat_inverse(4,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

          ! he4 he4 he4 <- c12
          if (( (rate_array(i)%parts(2) .eq. ihe4)        .and. &
                (rate_array(i)%parts(3) .eq. ihe4)        .and. &
                (rate_array(i)%parts(4) .eq. ihe4)        .and. &
                (rate_array(i)%parts(1) .eq. idx_c12))    .and. &
                (rate_array(i)%group .eq. 3)              .and. &
                (rate_array(i)%source .eq. "detb")) then
                call calculate_detailed_balance(temp_grid(j),rate_array,rate_array(i),rat_tmp)
                rat_inverse(8,j)       = rat_inverse(8,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

          ! n na20  <-  p ne20
          if ((((rate_array(i)%parts(3) .eq. ineu)       .and. &
                (rate_array(i)%parts(4) .eq. idx_na20))  .or.  &
               ((rate_array(i)%parts(3) .eq. idx_na20)   .and. &
                (rate_array(i)%parts(4) .eq. ineu)))     .and. &
              (((rate_array(i)%parts(1) .eq. ipro)       .and. &
                (rate_array(i)%parts(2) .eq. idx_ne20))  .or.  &
               ((rate_array(i)%parts(1) .eq. idx_ne20)   .and. &
                (rate_array(i)%parts(2) .eq. ipro)))     .and. &
                (rate_array(i)%group .eq. 5)             .and. &
                (rate_array(i)%source .eq. "detb")) then
                call calculate_detailed_balance(temp_grid(j),rate_array,rate_array(i),rat_tmp)
                rat_inverse(5,j)       = rat_inverse(5,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

          ! d  li7  <-  n  he4  he4
          if ((((rate_array(i)%parts(4) .eq. idx_d)     .and. &
                (rate_array(i)%parts(5) .eq. idx_li7))   .or. &
               ((rate_array(i)%parts(5) .eq. idx_d)     .and. &
                (rate_array(i)%parts(4) .eq. idx_li7))  .and. &
               ((rate_array(i)%parts(1) .eq. ineu)     .and. &
                (rate_array(i)%parts(2) .eq. ihe4)     .and. &
                (rate_array(i)%parts(3) .eq. ihe4))    .or.  &
               ((rate_array(i)%parts(1) .eq. ihe4)     .and. &
                (rate_array(i)%parts(2) .eq. ineu)     .and. &
                (rate_array(i)%parts(3) .eq. ihe4))    .or.  &
               ((rate_array(i)%parts(1) .eq. ihe4)     .and. &
                (rate_array(i)%parts(2) .eq. ihe4)     .and. &
                (rate_array(i)%parts(3) .eq. ineu)))   .and. &
                (rate_array(i)%group .eq. 9)           .and. &
                (rate_array(i)%source .eq. "detb")) then
                call calculate_detailed_balance(temp_grid(j),rate_array,rate_array(i),rat_tmp)
                rat_inverse(6,j)       = rat_inverse(6,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

          ! he3  li7 <-  n    p  he4  he4
          if ((((rate_array(i)%parts(5) .eq. idx_he3)    .and. &
                (rate_array(i)%parts(6) .eq. idx_li7)  .and. &
                (rate_array(i)%parts(1) .eq. ineu))    .and. &
                (rate_array(i)%parts(2) .eq. ipro)     .and. &
                (rate_array(i)%parts(3) .eq. ihe4)     .and. &
                (rate_array(i)%parts(4) .eq. ihe4))    .and. &
                (rate_array(i)%group .eq. 10)          .and. &
                (rate_array(i)%source .eq. "detb")) then
                call calculate_detailed_balance(temp_grid(j),rate_array,rate_array(i),rat_tmp)
                rat_inverse(7,j)       = rat_inverse(7,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

          ! n    p    p  <->   p    d
          if ((((rate_array(i)%parts(3) .eq. ineu)    .and. &
                (rate_array(i)%parts(4) .eq. ipro)    .and. &
                (rate_array(i)%parts(5) .eq. ipro))   .and. &
                (rate_array(i)%parts(1) .eq. ipro)    .and. &
                (rate_array(i)%parts(2) .eq. idx_d))  .and. &
                (rate_array(i)%group .eq. 6)          .and. &
                (rate_array(i)%source .eq. "detb")) then
                call calculate_detailed_balance(temp_grid(j),rate_array,rate_array(i),rat_tmp)
                rat_inverse(9,j)       = rat_inverse(9,j)+rat_tmp
                rate_array(i)%cached_p = rat_tmp
          end if

        end do

        do i=1,nrea_old
          if ((((rate_array_old(i)%parts(3) .eq. ihe4)      .and. &
                (rate_array_old(i)%parts(2) .eq. idx_ne20)) .or.  &
               ((rate_array_old(i)%parts(2) .eq. ihe4)      .and. &
                (rate_array_old(i)%parts(3) .eq. idx_ne20))) .and. &
                (rate_array_old(i)%parts(1) .eq. idx_mg24)  .and. &
                (rate_array_old(i)%group .eq. 2)) then
                call calculate_reacl_rate(rate_array_old(i),rat_tmp)
                rat_inverse_reacl(4,j) = rat_inverse_reacl(4,j)+rat_tmp
                rate_array_old(i)%cached_p = rat_tmp
          end if

          ! he4 he4 he4 <- c12
          if (( (rate_array_old(i)%parts(2) .eq. ihe4)      .and. &
                (rate_array_old(i)%parts(3) .eq. ihe4)      .and. &
                (rate_array_old(i)%parts(4) .eq. ihe4)     .and. &
                (rate_array_old(i)%parts(1) .eq. idx_c12))    .and. &
                (rate_array_old(i)%group .eq. 3))    then
                call calculate_reacl_rate(rate_array_old(i),rat_tmp)
                rat_inverse_reacl(8,j) = rat_inverse_reacl(8,j)+rat_tmp
                rate_array_old(i)%cached_p = rat_tmp
          end if

          if ((((rate_array_old(i)%parts(3) .eq. ineu)      .and. &
                (rate_array_old(i)%parts(4) .eq. idx_na20)) .or. &
               ((rate_array_old(i)%parts(3) .eq. idx_na20) .and. &
                (rate_array_old(i)%parts(4) .eq. ineu)))    .and. &
              (((rate_array_old(i)%parts(1) .eq. ipro) .and. &
                (rate_array_old(i)%parts(2) .eq. idx_ne20)) .or. &
               ((rate_array_old(i)%parts(1) .eq. idx_ne20) .and. &
                (rate_array_old(i)%parts(2) .eq. ipro))) .and. &
                (rate_array_old(i)%group .eq. 5))    then
                call calculate_reacl_rate(rate_array_old(i),rat_tmp)
                rat_inverse_reacl(5,j)     = rat_inverse_reacl(5,j)+rat_tmp
                rate_array_old(i)%cached_p = rat_tmp
          end if

          ! d  li7  <-   n  he4  he4
          if ((((rate_array_old(i)%parts(4) .eq. idx_d)    .and. &
                (rate_array_old(i)%parts(5) .eq. idx_li7)  .and. &
                (rate_array_old(i)%parts(1) .eq. ineu))    .and. &
                (rate_array_old(i)%parts(2) .eq. ihe4)     .and. &
                (rate_array_old(i)%parts(3) .eq. ihe4))    .and. &
                (rate_array_old(i)%group .eq. 9))    then
                call calculate_reacl_rate(rate_array_old(i),rat_tmp)
                rat_inverse_reacl(6,j)     = rat_inverse_reacl(6,j)+rat_tmp
                rate_array_old(i)%cached_p = rat_tmp
          end if

          ! he3  li7 <-  n    p  he4  he4
          if ((((rate_array_old(i)%parts(5) .eq. idx_he3)    .and. &
                (rate_array_old(i)%parts(6) .eq. idx_li7)  .and. &
                (rate_array_old(i)%parts(1) .eq. ineu))    .and. &
                (rate_array_old(i)%parts(2) .eq. ipro)     .and. &
                (rate_array_old(i)%parts(3) .eq. ihe4)     .and. &
                (rate_array_old(i)%parts(4) .eq. ihe4))    .and. &
                (rate_array_old(i)%group .eq. 10))    then
                call calculate_reacl_rate(rate_array_old(i),rat_tmp)
                rat_inverse_reacl(7,j)     = rat_inverse_reacl(7,j)+rat_tmp
                rate_array_old(i)%cached_p = rat_tmp
          end if

          ! n    p    p  <->   p    d
          if ((((rate_array_old(i)%parts(3) .eq. ineu)    .and. &
                (rate_array_old(i)%parts(4) .eq. ipro)  .and. &
                (rate_array_old(i)%parts(5) .eq. ipro))    .and. &
                (rate_array_old(i)%parts(1) .eq. ipro)     .and. &
                (rate_array_old(i)%parts(2) .eq. idx_d))     .and. &
                (rate_array_old(i)%group .eq. 6))    then
                call calculate_reacl_rate(rate_array_old(i),rat_tmp)
                rat_inverse_reacl(9,j)     = rat_inverse_reacl(9,j)+rat_tmp
                rate_array_old(i)%cached_p = rat_tmp
          end if

        end do

        ! Ensure the format works
        do i=1,11
          rat_inverse_reacl(i,j) = max(rat_inverse_reacl(i,j),1d-99)
          rat_inverse(i,j)       = max(rat_inverse(i,j),1d-99)
          rat_forward(i,j)       = max(rat_forward(i,j),1d-99)
        end do
      end do

      file_id = open_outfile("debug_detailed_balance.dat")
      write(file_id,"(A)")"# File to check inverse rates. The inverse rates calculated (I) and "
      write(file_id,"(A)")"# from the reaclib (IR) should agree. There might be slight deviations"
      write(file_id,"(A)")"# due to different Q-values in the Reaclib and winvn file."
      write(file_id,"(A)")"# Tested are a couple of chapters, e.g.,"
      write(file_id,"(A)")"# Chapter 4: Ne20 + He4      ->  Mg24"
      write(file_id,"(A)")"# Chapter 5: n + Na20        ->  p Ne20"
      write(file_id,"(A)")"# Chapter 6: d + Li7         ->  n + He4 + He4"
      write(file_id,"(A)")"# Chapter 7: He3 + Li7       ->  n + p + He4 + He4"
      write(file_id,"(A)")"# Chapter 8: He4 + He4 + He4 ->  C12"
      write(file_id,"(A)")"# Chapter 9: n + p + p       ->  p + d"
      write(file_id,"(A)")"###############################################"
      write(file_id,"(A)")"#  T [GK],      F (4),       I (4),       IR (4),"//&
                                    "      F (5),       I (5),       IR (5),"//&
                                    "      F (6),       I (6),       IR (6),"//&
                                    "      F (7),       I (7),       IR (7),"//&
                                    "      F (8),       I (8),       IR (8),"//&
                                    "      F (9),       I (9),       IR (9)"
      do i=1,301
        write(file_id,"(*(1pE13.4))")temp_grid(i),rat_forward(4,i),rat_inverse(4,i),rat_inverse_reacl(4,i),&
                                                  rat_forward(5,i),rat_inverse(5,i),rat_inverse_reacl(5,i),&
                                                  rat_forward(6,i),rat_inverse(6,i),rat_inverse_reacl(6,i),&
                                                  rat_forward(7,i),rat_inverse(7,i),rat_inverse_reacl(7,i),&
                                                  rat_forward(8,i),rat_inverse(8,i),rat_inverse_reacl(8,i),&
                                                  rat_forward(9,i),rat_inverse(9,i),rat_inverse_reacl(9,i)
      end do

      call close_io_file(file_id,"debug_detailed_balance.dat")

    end subroutine write_reac_rate


    !> Write the verbose output of the reaction rates
    !!
    !! The rates are always counted, for a certain verbose level they
    !! are also printed to the OUT file
    !!
    !! @author M. Reichert
    !! @date 21.07.22
    subroutine write_reac_verbose_out()
       use error_msg_class, only: int_to_str,write_data_to_std_out
       implicit none

       if (VERBOSE_LEVEL .ge. 1) then
          call write_data_to_std_out("Removing inverse rates",int_to_str(ninv_old))
          call write_data_to_std_out("Adding inverse rates",int_to_str(ninv))
       end if
    end subroutine write_reac_verbose_out


end module detailed_balance
