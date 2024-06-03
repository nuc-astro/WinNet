!> @file fission_rate_module.f90
!!
!! The error file code for this file is ***W19***.
!! Contains the module \ref fission_rate_module

!> Module to deal with fission reactions
!!
!! This module contains subroutines to read fission reactions
!! and fission fragment distributions. Furthermore, it counts the
!! amount of reactions. Many subroutines here were originally implemented
!! by M. Eichler and C. Winteler.
!!
!! @remark All subroutines have often been existed at other places in
!!         the network and possibly as parts of other subroutines.
#include "macros.h"
module fission_rate_module
   implicit none


   !> fission rate type, designed to save fragment distribution at initialization
   type,public                     :: fissionrate_type
      integer                      :: fissnuc_index        !< index of parent nucleus
      integer                      :: channels             !< number of fragment pairs for each fission rate
      integer                      :: mode                 !< fission mode; 1:n-induced fission 2: spontaneous and beta-delayed fission
      integer                      :: released_n           !< Amount of released neutrons for beta-delayed fission
      integer                      :: dimens               !< dimension of fissparts array
      real(r_kind)                 :: cached               !< computed rate
      real(r_kind)                 :: averageQ             !< average Q value, weighted over all fragment channels
      character(4)                 :: src                  !< Source of the reaction
      integer                      :: reac_type            !< "rrt_sf": spontaneous fission, "rrt_bf": beta-delayed fission, "rrt_nf": neutron-induced fission
      integer,dimension(:),allocatable      :: fissparts   !< corresponds to rrate()%parts() array
      integer,dimension(:,:),allocatable    :: cscf_ind    !< cscf_ind (i,j) gives the position of the entry \f$\frac{d\dot{Y}_{j}}{dY_{i}}\f$ in the cscf data array
      real(r_kind),dimension(9)             :: param       !< rate parameters a0-a9, as given in reaclib file
      real(r_kind),dimension(:),allocatable :: q_value     !< Q value of each fission channel
      real(r_kind),dimension(:),allocatable :: channelprob !< probability for each fragment pair
      real(r_kind),dimension(:),allocatable :: ch_amount   !< corresponds to rrate(i)%ch_amount(j); for fragment i ch_amount(i) = +1 * channelprob(i)
   end type fissionrate_type

   type(fissionrate_type),dimension(:),allocatable,public   :: fissrate !< Array storing fission reactions

   real(r_kind),dimension(:,:),allocatable,private :: beta_delayed_fiss_probs
   integer,private  :: amount_cols !< Amount of columns in the beta-delayed fission file


   character(len=*), private, parameter       :: fiss_binary_name='fiss_rates.windat' !< Name of the binary file containing the fission rates


    type,private                              :: fragtype
        integer                               :: Zp       !< Z of fissioning nucleus
        integer                               :: Ap       !< A of fissioning nucleus
        integer                               :: neutrons !< number of neutrons emitted
        integer                               :: nr_frags !< Number of fragments
        real(r_kind)                          :: Yn       !< Yield for neutrons
        integer,dimension(:),allocatable      :: net_idx  !< index of fragment in network
        integer,dimension(:),allocatable      :: Z        !< Z of fragment
        integer,dimension(:),allocatable      :: A        !< A of fragment
        real(r_kind),dimension(:),allocatable :: Y        !< Yields of fragment
    end type fragtype
    type(fragtype),dimension(:),allocatable,private   :: fissfrags_n_induced    !< Array storing fragment distributions of neutron-induced fission from file
    type(fragtype),dimension(:),allocatable,private   :: fissfrags_spontaneous  !< Array storing fragment distributions of spontaneous fission from file
    type(fragtype),dimension(:),allocatable,private   :: fissfrags_beta_delayed !< Array storing fragment distributions of beta-delayed fission from file




   integer, public, parameter    :: fiss_neglect=5 !< Amount of fission fragments not to be neglected
                                                   !  in the jacobian in case of vanishing parent
   integer, public               :: nfiss        !< Amount of fission rates
   integer, private              :: nfiss_spont  !< Amount of spontaneous fission rates
   integer, private              :: nfiss_n_ind  !< Amount of neutron induced fission rates
   integer, private              :: nfiss_bdel   !< Amount of beta-delayed fission rates
   integer, private              :: nufiss       !< total number of fission neutrons

   integer, private  :: n_nf,n_bf,n_sf !< Amount of individual reactions
   !
   ! Public and private fields and methods of the module
   !
   public:: &
      init_fission_rates, merge_fission_rates, output_binary_fission_reaction_data
   private:: &
      count_fission_rates, read_fission_rates, fiss_dist, kodtakdist,&
      abla_nfiss, abla_betafiss, read_binary_fission_reaction_data,&
      reorder_fragments, read_fission_rates_halflife_format, read_fission_rates_probability_format,&
      read_fission_rates_reaclib_format, count_fission_rates_halflife_format,&
      count_fission_rates_probability_format, count_fission_rates_reaclib_format,&
      modify_halflifes, write_reac_verbose_out
contains


   !> Initialize the fission reactions
   !!
   !! This subroutine counts and reads fission reactions.
   !! After calling it \ref fissrate will be filled as well as the integer
   !! \ref nfiss will store the amount of fission reactions.
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine init_fission_rates()
      use error_msg_class, only: raise_exception
      use parameter_class, only: fissflag, use_prepared_network,&
                                 prepared_network_path
      implicit none
      integer   ::  alloc_stat !< Allocation status flag

      n_nf=0; n_bf=0; n_sf=0; ! Set counters to 0
      nfiss=0

      if (fissflag.ne.0) then
        if (use_prepared_network) then
          call read_binary_fission_reaction_data(prepared_network_path)
        else
         !-- Count the amount of fission rates
         call count_fission_rates()

         !-- Allocate the fission rate arrays
         allocate(fissrate(nfiss),stat=alloc_stat)
         if ( alloc_stat /= 0) call raise_exception('Allocation of "fissrate" failed.',&
                                                    "init_fission_rates",190001)

         !-- Read the fission rates into both arrays
         call read_fission_rates()
         call add_fission_fragments()

         !-- Reorder the fragments, having the ones with the largest
         !!  ch_amount at the beginning.
         call reorder_fragments()

        end if
        !-- Write an overview of the reactions
        call write_reac_verbose_out()
      end if


   end subroutine init_fission_rates



   !> Reorder the fission fragments
   !!
   !! Subroutine to reorder the fragments, having the ones with the largest
   !! ch_amount at the beginning. This is benificial as we can ignore
   !! some of the entries in the jacobian in case the parent abundance
   !! is zero.
   !!
   !! @author M. Reichert
   !! @date 22.02.23
   subroutine reorder_fragments()
    use mergesort_module, only: bubblesort, reorder_int
    implicit none
    integer                          :: i !< loop variable
    integer,dimension(:),allocatable :: indices

    INFO_ENTRY("reorder_fragments")

    do i=1,nfiss
        ! Dont do this for negligible amount of fragments
        if (fissrate(i)%dimens .le. fiss_neglect+2) cycle

        ! Allocate indice array
        if (allocated(indices)) deallocate(indices)
        allocate(indices(fissrate(i)%dimens-2))

        ! Order descending according to ch_amount,
        ! start with the 3rd entry, since the first two are the fissioning nuclei
        call bubblesort(0,fissrate(i)%dimens-2,fissrate(i)%ch_amount(3:),indices)
        ! Also reorder the fissparts array
        call reorder_int(fissrate(i)%fissparts(3:),indices,fissrate(i)%dimens-2)

    end do

    INFO_EXIT("reorder_fragments")

   end subroutine reorder_fragments


   !> Write the amount of individual reactions to the out
   !!
   !! The rates are always counted, for a certain verbose level they
   !! are also printed to the OUT file
   !!
   !! @author M. Reichert
   !! @date 27.01.21
   subroutine write_reac_verbose_out()
      use error_msg_class, only: int_to_str,write_data_to_std_out
      implicit none
      character(len=7) :: tmp !< temporary character for pretty output

      if (VERBOSE_LEVEL .ge. 1) then
         call write_data_to_std_out("Amount fission rates",int_to_str(nfiss))
         call write_data_to_std_out("Total number of fission neutrons",int_to_str(nufiss))
      elseif (VERBOSE_LEVEL .ge. 2) then
         if (nfiss .gt. 0) write(*,"(A)") ""
         if (nfiss .gt. 0) write(*,"(A)") "    Fission rates:  "
         if (nfiss .gt. 0) write(*,"(A)") "   |----------------------|"
         tmp = int_to_str(nfiss)
         if (nfiss .gt. 0) write(*,"(A)") "   | Total       :"//adjustr(tmp)//" |"
         tmp = int_to_str(n_nf)
         if (n_nf .gt. 0) write(*,"(A)")  "   | n-induced   :"//adjustr(tmp)//" |"
         tmp = int_to_str(n_bf)
         if (n_bf .gt. 0) write(*,"(A)")  "   | b-delayed   :"//adjustr(tmp)//" |"
         tmp = int_to_str(n_sf)
         if (n_sf .gt. 0) write(*,"(A)")  "   | spontaneous :"//adjustr(tmp)//" |"
         if (nfiss .gt. 0) write(*,"(A)") "   |----------------------|"
         if (nfiss .gt. 0) write(*,"(A)") ""
      end if

   end subroutine write_reac_verbose_out



   !> Read the fission reactions and fragment distributions from a binary file
   !!
   !! This subroutine reads the fission reactions and fragment distributions
   !! from a binary, unformatted file. The file is assumed to be created
   !! beforehand. It is only read when use_prepared_network is set to true.
   !!
   !! @author M. Reichert
   !! @date 21.07.23
   subroutine read_binary_fission_reaction_data(path)
    use file_handling_class, only: open_unformatted_infile
    use parameter_class,     only: max_fname_len, fissflag
    use error_msg_class,     only: raise_exception
    implicit none
    character(len=*), intent(in) :: path         !< Path to folder with binary files
    integer                      :: i            !< Loop variable
    integer                      :: file_id      !< File id
    integer                      :: alloc_stat   !< Allocation status
    logical                      :: is_allocated !< Logical to check if array is allocated

    if (fissflag .eq. 0) return

    ! Open an unformatted file
    file_id = open_unformatted_infile(trim(adjustl(path))//trim(adjustl(fiss_binary_name)))

    read(file_id) nfiss
    read(file_id) nufiss
    read(file_id) n_nf
    read(file_id) n_bf
    read(file_id) n_sf

    ! Allocate the fission rate array
    allocate(fissrate(nfiss),stat=alloc_stat)
    if ( alloc_stat /= 0) call raise_exception('Allocation of "fissrate" failed.',&
                                               "read_binary_fission_reaction_data",190001)

    do i=1,nfiss
      ! Write all fission rates
      read(file_id) fissrate(i)%fissnuc_index
      read(file_id) fissrate(i)%channels
      read(file_id) fissrate(i)%mode
      read(file_id) fissrate(i)%src
      read(file_id) fissrate(i)%dimens
      read(file_id) fissrate(i)%cached
      read(file_id) fissrate(i)%averageQ
      read(file_id) fissrate(i)%reac_type

      read(file_id) is_allocated
      if (is_allocated) then
        ! Allocate fissparts, cscf_ind, q_value, channelprob, and ch_amount array
        allocate(fissrate(i)%fissparts(fissrate(i)%dimens),stat=alloc_stat)
        if ( alloc_stat /= 0) call raise_exception('Allocation of "fissrate(i)%fissparts" failed.',&
                                                    "read_binary_fission_reaction_data",190001)
        read(file_id) fissrate(i)%fissparts   ! has dimension dimens
      end if

      read(file_id) is_allocated
      if (is_allocated) then
          allocate(fissrate(i)%cscf_ind(fissrate(i)%dimens,fissrate(i)%dimens),stat=alloc_stat)
          if ( alloc_stat /= 0) call raise_exception('Allocation of "fissrate(i)%cscf_ind" failed.',&
                                                      "read_binary_fission_reaction_data",190001)
          read(file_id) fissrate(i)%cscf_ind    ! has dimension (dimens,dimens)
      end if

      read(file_id) fissrate(i)%param       ! has dimension 9

      read(file_id) is_allocated
      if (is_allocated) then
          allocate(fissrate(i)%q_value(fissrate(i)%channels),stat=alloc_stat)
          if ( alloc_stat /= 0) call raise_exception('Allocation of "fissrate(i)%q_value" failed.',&
                                                      "read_binary_fission_reaction_data",190001)
          read(file_id) fissrate(i)%q_value     ! has dimension channels
      end if

      read(file_id) is_allocated
      if (is_allocated) then
          allocate(fissrate(i)%channelprob(fissrate(i)%channels),stat=alloc_stat)
          if ( alloc_stat /= 0) call raise_exception('Allocation of "fissrate(i)%channelprob" failed.',&
                                                      "read_binary_fission_reaction_data",190001)
          read(file_id) fissrate(i)%channelprob ! has dimension channels
      end if

      read(file_id) is_allocated
      if (is_allocated) then
          allocate(fissrate(i)%ch_amount(fissrate(i)%dimens),stat=alloc_stat)
          if ( alloc_stat /= 0) call raise_exception('Allocation of "fissrate(i)%ch_amount" failed.',&
                                                      "read_binary_fission_reaction_data",190001)
          read(file_id) fissrate(i)%ch_amount   ! has dimension dimens
      end if


    end do

    close(file_id)


    end subroutine read_binary_fission_reaction_data


   !> Save the fission data to a unformatted binary file
   !!
   !! This subroutine saves the fission data to a unformatted binary file.
   !!
   !! @author M. Reichert
   !! @date 21.07.23
   subroutine output_binary_fission_reaction_data(path)
    use file_handling_class, only: open_unformatted_outfile
    use parameter_class,     only: fissflag
    implicit none
    character(len=*), intent(in) :: path
    integer                      :: i
    integer                      :: file_id

    if (fissflag .eq. 0) return
    ! Open an unformatted file
    file_id = open_unformatted_outfile(trim(adjustl(path))//trim(adjustl(fiss_binary_name)))
    ! Write the fission rates to binary file
    write(file_id) nfiss
    write(file_id) nufiss
    write(file_id) n_nf
    write(file_id) n_bf
    write(file_id) n_sf


    do i=1,nfiss
      ! Write all fission rates
      write(file_id) fissrate(i)%fissnuc_index
      write(file_id) fissrate(i)%channels
      write(file_id) fissrate(i)%mode
      write(file_id) fissrate(i)%src
      write(file_id) fissrate(i)%dimens
      write(file_id) fissrate(i)%cached
      write(file_id) fissrate(i)%averageQ
      write(file_id) fissrate(i)%reac_type
      write(file_id) allocated(fissrate(i)%fissparts)
      if (allocated(fissrate(i)%fissparts)) then
        write(file_id) fissrate(i)%fissparts   ! has dimension dimens
      end if
      write(file_id) allocated(fissrate(i)%cscf_ind)
      if (allocated(fissrate(i)%cscf_ind)) then
        write(file_id) fissrate(i)%cscf_ind    ! has dimension (dimens,dimens)
      end if
      write(file_id) fissrate(i)%param       ! has dimension 9
      write(file_id) allocated(fissrate(i)%q_value)
      if (allocated(fissrate(i)%q_value)) then
        write(file_id) fissrate(i)%q_value     ! has dimension channels
      end if
      write(file_id) allocated(fissrate(i)%channelprob)
      if (allocated(fissrate(i)%channelprob)) then
          write(file_id) fissrate(i)%channelprob ! has dimension channels
      end if
      write(file_id) allocated(fissrate(i)%ch_amount)
      if (allocated(fissrate(i)%ch_amount)) then
        write(file_id) fissrate(i)%ch_amount   ! has dimension dimens
      end if
    end do

    ! Close the file again
    close(file_id)

   end subroutine output_binary_fission_reaction_data


   !> Merge fission rates with larger array.
   !!
   !! This subroutine will merge the fission rates.
   !! Since fission rates are treated separately, there is not really much
   !! to do. However, in case of beta-delayed fission, the half lifes
   !! of the beta decays have to be modified to include the fission rates.
   !!
   !! @note The other fission rate array \ref fissrate will exist independently
   !!       and are used to calculate the correct equations for the
   !!       fission fragments in \ref jacobian_class::jacobi_init,
   !!       \ref jacobian_class::abchange.
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine merge_fission_rates(rrate_array,rrate_length,fiss_count)
      use error_msg_class, only: raise_exception
      use parameter_class, only: use_prepared_network, fission_format_beta_delayed
      use global_class,    only: reactionrate_type
      implicit none
      type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate_array  !< Large rate array, containing all reactions
      integer,intent(inout)                                          :: rrate_length !< length of rrate_array
      integer,intent(out)                                            :: fiss_count   !< Amount of fission rates
      integer                                                        :: new_length   !< New length of rrate_array

      !-- Only do something if there are fission rates
      if (allocated(fissrate) .and. (.not. use_prepared_network)) then
         ! New length of the array
         new_length = rrate_length
         if (nfiss .ne. 0) then
            ! Prepare something for the case that fission_format_beta_delayed is
            ! Set to probability format. In this case, the halflifes of the beta decays
            ! have to be modified.
            if (fission_format_beta_delayed .eq. 3) then
                call modify_halflifes(rrate_array,rrate_length,fissrate,nfiss)
                ! The length may have gotten modified
                new_length = rrate_length
            end if
           !-- Output the new length
           rrate_length = new_length
         end if
      end if
      fiss_count = nfiss
   end subroutine merge_fission_rates


   !> Modifies half lifes of beta decays to include fission
   !!
   !! In case the beta delayed fission rates are given as a probability,
   !! the beta decays are rescaled to include the fission rates without
   !! changing the total half life.
   !! The fission rates are given by:
   !! \f[
   !! \lambda_{\text{fiss}} =P_{\text{fiss}} * \lambda_{\beta}
   !! \f]
   !! where \f$P_{\text{fiss}}\f$ is the probability of fission and
   !! \f$\lambda_{\beta}\f$ is the beta decay rate as given initially
   !! in the Reaclib. The beta decay rate is then modified to
   !! \f[
   !! \lambda_{\beta, new} = \left(1-P_{\text{fiss}}\right) * \lambda_{\beta}
   !! \f]
   !! There are special cases to consider, namely \f$P_{\text{fiss}}=1\f$ in
   !! which case the reaclib rate has to be removed completely and
   !! \f$\lambda_{\beta}\f = 0$ in which case the fission rate has to be removed.
   !! Note also that fission can have multiple channels, which are by emitting
   !! different amount of neutrons.
   !! The fission rates in form of probabilities are given, e.g., in
   !! [Panov et al. 2005](https://ui.adsabs.harvard.edu/abs/2005NuPhA.747..633P/abstract) or
   !! [Mumpower et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...869...14M/abstract).
   !!
   !! @author M. Reichert
   !! @date 30.05.24
   subroutine modify_halflifes(rrate,rrate_length,fissrate_in,nfiss_in)
    use global_class,    only: reactionrate_type, net_size, ineu, ipro, net_names
    use nucstuff_class,  only: get_nr_reactants
    use error_msg_class, only: int_to_str, raise_exception
    implicit none
    type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate            !< Large rate array, containing all reactions
    integer,intent(inout)                                          :: rrate_length     !< length of rrate_array
    type(fissionrate_type),dimension(:),allocatable,intent(inout)  :: fissrate_in      !< Fission rate array
    integer,intent(inout)                                          :: nfiss_in         !< Amount of fission rates
    integer                                            :: i                !< Loop variable
    real(r_kind),dimension(net_size)                   :: lambdas          !< Lambdas of the beta decay
    integer                                            :: j                !< Loop variable
    real(r_kind)                                       :: tot_fiss_prob    !< Total fission probability
    integer                                            :: parent_idx       !< Index of fissioning nucleus
    integer                                            :: Px_idx           !< Channel of the fission rate
    real(r_kind)                                       :: Px               !< Probability of the fission rate
    logical,dimension(rrate_length)                    :: rrate_mask       !< Mask for the rrate array
    logical,dimension(nfiss_in)                        :: fissrate_mask    !< Mask for the fissrate array
    type(reactionrate_type),dimension(:),allocatable   :: rrate_tmp        !< Temporary array to store the rates
    type(fissionrate_type),dimension(:),allocatable    :: fissrate_tmp     !< Temporary array to store the fission rates
    integer                                            :: new_rrate_length !< New length of the rrate array
    integer                                            :: new_nfiss_in     !< New amount of fission rates
    integer                                            :: alloc_stat       !< Allocation status
    real(r_kind)                                       :: lambda_rate      !< Rate of the beta decay

    INFO_ENTRY("modify_halflifes")

    lambdas(:) = 0d0

    fissrate_mask(:) = .True.

    ! Ensure that there are reactions in rrate
    if (allocated(rrate)) then
        rrate_mask(:) = .True.
        ! Modify beta decays in rrate array and get the half lifes of the beta decays
        do i=1,rrate_length

            ! Ignore things that are definitely not beta decays
            ! and also theoretical weak rates
            if ((rrate(i)%reac_src .eq. rrs_twr)   .or. &
                (rrate(i)%reac_src .eq. rrs_nu)    .or. &
                (rrate(i)%reac_src .eq. rrs_fiss)  .or. &
                (rrate(i)%reac_src .eq. rrs_aext)  .or. &
                (rrate(i)%reac_src .eq. rrs_detb)) cycle

            if (rrate(i)%reac_type.eq.rrt_betm) then
                part_loop: do j=1,get_nr_reactants(rrate(i)%group)
                    if ((rrate(i)%parts(j).ne.ineu) .and. (rrate(i)%parts(j).ne.ipro)) then
                        ! Get the total fission fraction over all channels (P0n, P1n, ...)
                        tot_fiss_prob = sum(beta_delayed_fiss_probs(rrate(i)%parts(j),:))

                        ! Take care of reaclib type reactions
                        if ((rrate(i)%reac_src .eq. rrs_reacl) .or. &
                            (rrate(i)%reac_src .eq. rrs_wext)) then
                            lambda_rate =  dexp(rrate(i)%param(1))
                            lambdas(rrate(i)%parts(j)) = lambdas(rrate(i)%parts(j))+lambda_rate
                            if ((1d0-tot_fiss_prob) .ne. 0d0) then
                                rrate(i)%param(1) = rrate(i)%param(1)+dlog(1d0-tot_fiss_prob)
                            else
                                ! Remove the rate, it will become a fission rate
                                rrate_mask(i) = .False.
                            end if
                        ! Take care of tabulated rates
                        else if (rrate(i)%reac_src .eq. rrs_tabl)  then
                            lambda_rate = rrate(i)%tabulated(1)
                            lambdas(rrate(i)%parts(j)) = lambdas(rrate(i)%parts(j))+lambda_rate
                            if ((1d0-tot_fiss_prob) .ne. 0d0) then
                                rrate(i)%tabulated(:) = rrate(i)%tabulated(:)*(1d0-tot_fiss_prob)
                            else
                                ! Remove the rate, it will become a fission rate
                                rrate_mask(i) = .False.
                            end if
                        else
                            call raise_exception('Unknown reac_src in modify_halflifes.',&
                                                 "modify_halflifes",190014)
                        end if
                        ! Go out of the loop
                        exit part_loop
                    end if
                end do part_loop
            end if
        end do

        ! Reduce the rrate array to the correct size
        new_rrate_length = count(rrate_mask)
        if (new_rrate_length .ne. rrate_length) then

            ! Say something if the verbose level is high enough
            if (VERBOSE_LEVEL .ge. 2) then
                write(*,*) "Reducing the rate array to the correct size by removing "//&
                           int_to_str(rrate_length-new_rrate_length)//&
                           " beta-decays and putting it into fission rates."
            end if

            ! Allocate a temporary array to store the content of rrate
            allocate(rrate_tmp(new_rrate_length),stat=alloc_stat)
            if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_tmp" failed.',&
                                                       "modify_halflifes",190001)

            rrate_tmp(1:new_rrate_length) = pack(rrate,rrate_mask)
            deallocate(rrate,stat=alloc_stat)
            if ( alloc_stat /= 0) call raise_exception('Deallocation of "rrate" failed.',&
                                                       "modify_halflifes",190002)
            allocate(rrate(new_rrate_length),stat=alloc_stat)
            if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate" failed.',&
                                                       "modify_halflifes",190001)

            rrate(1:new_rrate_length) = rrate_tmp(1:new_rrate_length)
            deallocate(rrate_tmp,stat=alloc_stat)
            if ( alloc_stat /= 0) call raise_exception('Deallocation of "rrate_tmp" failed.',&
                                                       "modify_halflifes",190002)
            rrate_length = new_rrate_length
        end if
    end if

    ! Make the correct lambdas in the fission array
    do i=1,nfiss_in
        if (fissrate_in(i)%reac_type.eq.rrt_bf) then
            ! Beta-delayed fission
            parent_idx = fissrate_in(i)%fissnuc_index
            Px_idx = int(fissrate_in(i)%released_n+1d0)
            Px = beta_delayed_fiss_probs(parent_idx,Px_idx)
            if (lambdas(parent_idx) .eq. 0d0) then
                ! Remove the rate, the parent is not beta-decaying
                fissrate_mask(i) = .False.
            end if
            ! Fissrate
            fissrate_in(i)%param(:) = 0d0
            fissrate_in(i)%param(1) = dlog(Px * lambdas(parent_idx))
        end if
    end do

    ! Reduce the fission rate array to the correct size
    new_nfiss_in = count(fissrate_mask)
    if (new_nfiss_in .ne. nfiss_in) then
        ! Say something if the verbose level is high enough
        if (VERBOSE_LEVEL .ge. 2) then
            write(*,*) "Reducing the fission rate array to the correct size by removing "//&
                       int_to_str(nfiss_in-new_nfiss_in)//" beta-delayed fission rates."
        end if
        ! Allocate a temporary array to store the content of fissrate
        allocate(fissrate_tmp(new_nfiss_in),stat=alloc_stat)
        if ( alloc_stat /= 0) call raise_exception('Allocation of "fissrate_tmp" failed.',&
                                                   "modify_halflifes",190001)

        fissrate_tmp(1:new_nfiss_in) = pack(fissrate_in,fissrate_mask)
        deallocate(fissrate_in,stat=alloc_stat)
        if ( alloc_stat /= 0) call raise_exception('Deallocation of "fissrate_in" failed.',&
                                                   "modify_halflifes",190002)
        allocate(fissrate_in(new_nfiss_in),stat=alloc_stat)
        if ( alloc_stat /= 0) call raise_exception('Allocation of "fissrate_in" or failed.',&
                                                   "modify_halflifes",190001)

        fissrate_in(1:new_nfiss_in) = fissrate_tmp(1:new_nfiss_in)
        deallocate(fissrate_tmp,stat=alloc_stat)
        if ( alloc_stat /= 0) call raise_exception('Deallocation of "fissrate_tmp" or failed.',&
                                                   "modify_halflifes",190002)
        nfiss_in = new_nfiss_in
    end if

    ! Deallocate the beta-delayed fission probabilities
    ! as everything is now in the fission rate array
    deallocate(beta_delayed_fiss_probs,stat=alloc_stat)
    if ( alloc_stat /= 0) call raise_exception('Deallocation of "beta_delayed_fiss_probs" failed.',&
                                               "modify_halflifes",190002)

    INFO_EXIT("modify_halflifes")
   end subroutine modify_halflifes



   !> Count the amount of fission rates
   !!
   !! This subroutine counts the amount of fission rates and
   !! stores the result in \ref nfiss.
   !! In case that a new fission type will be implemented,
   !! the reaction rates should be counted here.
   !! The different formats for the different types of fission reactions are:
   !! <table>
   !! <caption id="multi_row">Fission formats</caption>
   !! <tr><th> Type                <th> Value <th> Description
   !! <tr><td> Spontaneous         <td> 0     <td> No rates are read
   !! <tr><td> Spontaneous         <td> 1     <td> Reaclib format
   !! <tr><td> Spontaneous         <td> 2     <td> Half life format
   !! <tr><td> n-induced           <td> 0     <td> No rates are read
   !! <tr><td> n-induced           <td> 1     <td> Reaclib format
   !! <tr><td> \f$\beta\f$-delayed <td> 0     <td> No rates are read
   !! <tr><td> \f$\beta\f$-delayed <td> 1     <td> Reaclib format
   !! <tr><td> \f$\beta\f$-delayed <td> 2     <td> Half life format
   !! <tr><td> \f$\beta\f$-delayed <td> 3     <td> Probability format
   !! </table>
   !!
   !! @see \ref read_fission_rates
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine count_fission_rates()
    use parameter_class, only: fission_rates_spontaneous, fission_format_spontaneous, &
                               fission_rates_n_induced, fission_format_n_induced, &
                               fission_rates_beta_delayed, fission_format_beta_delayed
    use global_class,    only: net_size
    use error_msg_class, only: int_to_str, raise_exception
    implicit none
    integer :: alloc_stat !< Allocation status

    INFO_ENTRY("count_fission_rates")

    ! Spontaneous fission rates
    select case(fission_format_spontaneous)
    case(0)
        nfiss_spont = 0
    case(1)
        call count_fission_rates_reaclib_format(fission_rates_spontaneous,nfiss_spont)
    case(2)
        call count_fission_rates_halflife_format(fission_rates_spontaneous,nfiss_spont)
    case default
        call raise_exception("Fission format for spontaneous fission rates not implemented, got '"//&
                              int_to_str(fission_format_spontaneous)//"'. Change the "//&
                              "'fission_format_spontaneous' parameter.",&
                              "count_fission_rates",190010)
    end select

    ! Neutron induced fission rates
    select case(fission_format_n_induced)
    case(0)
        nfiss_n_ind = 0
    case(1)
        call count_fission_rates_reaclib_format(fission_rates_n_induced,nfiss_n_ind)
    case default
        call raise_exception("Fission format for neutron-induced fission rates not implemented, got '"//&
                              int_to_str(fission_format_n_induced)//"'. Change the "//&
                              "'fission_format_n_induced' parameter.",&
                              "count_fission_rates",190010)
    end select

    ! Beta-delayed fission rates
    select case(fission_format_beta_delayed)
    case(0)
        nfiss_bdel = 0
    case(1)
        call count_fission_rates_reaclib_format(fission_rates_beta_delayed,nfiss_bdel)
    case(2)
        call count_fission_rates_halflife_format(fission_rates_beta_delayed,nfiss_bdel)
    case(3)
        call count_fission_rates_probability_format(fission_rates_beta_delayed,nfiss_bdel,amount_cols)
        ! Initialize an array to store the beta-delayed fission probabilities
        allocate(beta_delayed_fiss_probs(net_size,amount_cols),stat=alloc_stat)
        if ( alloc_stat /= 0) call raise_exception('Allocation of "beta_delayed_fiss_probs" failed.',&
                                                   "count_fission_rates",190001)
        beta_delayed_fiss_probs(:,:) = 0d0
    case default
        call raise_exception("Fission format for beta-delayed fission rates not implemented, got '"//&
                              int_to_str(fission_format_beta_delayed)//"'. Change the "//&
                              "'fission_format_beta_delayed' parameter.",&
                              "count_fission_rates",190010)
    end select

    ! Calculate the total amount of fission rates
    nfiss = nfiss_spont + nfiss_n_ind + nfiss_bdel

    INFO_EXIT("count_fission_rates")
   end subroutine count_fission_rates


   !> Count the amount of fission rates in reaclib format
   !!
   !! This subroutine counts the amount of fission rates in reaclib format
   !! and stores the result in count_rates. This is necessary to allocate
   !! the fission array with the correct length. The subroutine can be called
   !! for any kind of fission (spontaneous, n-induced, beta-delayed).
   !!
   !! \b Edited:
   !!    - 30.05.2024 (MR): Made this its own subroutine to add more file formats
   !! .
   subroutine count_fission_rates_reaclib_format(fission_rate_file,count_rates)
    use file_handling_class
    use benam_class,     only: benam
    use format_class
    implicit none
    character(len=*), intent(in)  :: fission_rate_file !< File containing fission rates
    integer, intent(out)          :: count_rates       !< Amount of fission rates
    integer                    :: j           !< Loop variable
    integer                    :: read_stat   !< Read status
    integer                    :: fissionlib  !< File id of fisionrates
    character(5), dimension(6) :: parts       !< Participating nuclei names
    integer, dimension(6)      :: parts_index !< Participating nuclei indices
    real(r_kind), dimension(9) :: params      !< Reaclib fitting parameter
    character(4)               :: src         !< Reaclib source label
    character(1)               :: res         !< Reaclib weak flag label
    character(1)               :: rev         !< Reaclib reverse label
    real(r_kind)               :: q           !< Reaclib Q-Value
    integer                    :: grp         !< Reaclib chapter
    integer                    :: ifiss       !< reaction rate counter

    ! Counter for fission rates
    ifiss = 0

    ! Open fission file
    fissionlib = open_infile(fission_rate_file)

    ! Count the rates
    fission: do
        read(fissionlib,my_format(1), iostat = read_stat)  &
                grp, parts(1:6), src, res, rev, q
        if (read_stat /= 0) exit
        read(fissionlib,"(4E13.6)") params(1:4)
        read(fissionlib,"(5E13.6)") params(5:9)
        if (grp.ne.0) cycle fission
        parts_index = 0
        indices: do j=1,6
            if (parts(j) .eq. '     ') exit indices
            parts_index(j) = benam(parts(j))
            ! parts_index(i) = 0 if nuclide i is not in network
            if (parts_index(j) .eq. 0) cycle fission
        end do indices
        ifiss = ifiss + 1
    end do fission

    ! Close the fission file again
    close(fissionlib)

    count_rates = ifiss
   end subroutine count_fission_rates_reaclib_format


   !> Count the amount of fission rates in half-life format
   !!
   !! This subroutine counts the amount of fission rates in half-life format
   !! and stores the result in count_rates. This is necessary to allocate
   !! the fission array with the correct length. The subroutine can be called
   !! for spontaneous and beta-delayed fission.
   !! An example of an entry can be:
   !! \file{
   !!  th198    4.760787e-09    kh20
   !!  th199    2.808520e-08    kh20
   !!  th200    2.841403e-08    kh20
   !! ...}
   !! where the first columns gives the nucleus, the second the halflife in seconds
   !! and the third the source. Note that in this format, no beta-delayed neutrons
   !! can be specified.
   !!
   !! @author M. Reichert
   !! @date 30.05.24
   subroutine count_fission_rates_halflife_format(fission_rate_file,count_rates)
    use file_handling_class
    use benam_class,     only: benam
    use format_class
    implicit none
    character(len=*), intent(in)  :: fission_rate_file !< File containing fission rates
    integer, intent(out)          :: count_rates       !< Amount of fission rates
    integer                    :: j           !< Loop variable
    integer                    :: read_stat   !< Read status
    integer                    :: fissionlib  !< File id of fisionrates
    character(5), dimension(6) :: parts       !< Participating nuclei names
    integer, dimension(6)      :: parts_index !< Participating nuclei indices
    integer                    :: ifiss       !< reaction rate counter
    real(r_kind), dimension(9) :: params      !< Reaclib fitting parameter
    character(4)               :: src         !< Reaclib source label

    ! Counter for fission rates
    ifiss = 0

    ! Open fission file
    fissionlib = open_infile(fission_rate_file)

    ! Count the rates
    fission_hl: do
        ! Set parts to empty strings
        parts(1:6) = '     '
        params(:)  = 0d0

        ! An example line could be: pa231    3.446020e+23    kh20
        read(fissionlib, *, iostat = read_stat)  &
                parts(1), params(1), src
        if (read_stat /= 0) exit

        parts_index = 0
        indices_hl: do j=1,6
            if (parts(j) .eq. '     ') exit indices_hl
            parts_index(j) = benam(parts(j))
            ! parts_index(i) = 0 if nuclide i is not in network
            if (parts_index(j) .eq. 0) cycle fission_hl
        end do indices_hl
        ifiss = ifiss + 1
    end do fission_hl

    ! Close the fission file again
    close(fissionlib)

    ! Store the amount of fission rates
    count_rates = ifiss
   end subroutine count_fission_rates_halflife_format


   !> Count the amount of fission rates in probability format
   !!
   !! This subroutine counts the amount of fission rates in probability format
   !! and stores the result in count_rates. This is necessary to allocate
   !! the fission array with the correct length. The subroutine can be called
   !! for beta-delayed fission.
   !! An example of an entry can be:
   !! \file{
   !! pa291        mp18
   !! 0.000    0.008    0.000    0.015    0.000    0.018    0.000    0.025    0.000    0.003
   !! pa292        mp18
   !! 0.000    0.001    0.002    0.000    0.004    0.000    0.006    0.000    0.018    0.000
   !! pa293        mp18
   !! 0.239    0.001    0.000    0.009    0.000    0.010    0.000    0.014    0.000    0.009
   !! ...}
   !! where the first line of an entry gives nucleus name and the source. The second line
   !! gives the probabilities of the fission channels (i.e., P0nf, P1nf, ...).
   !! Note that there is a maximum number of columns (i.e., neutron emission) set to 30 here.
   !! Throughout the file, the amount of possible channels have to stay the same, but they
   !! can be different for two different files. For example, for Panov et al. 2005 the file
   !! could look like:
   !! \file{
   !!  u259        pa05
   !! 0.990
   !!  u260        pa05
   !! 0.020
   !!  u261        pa05
   !! 0.990
   !! ...}.
   !!
   !! @author M. Reichert
   !! @date 30.05.24
   subroutine count_fission_rates_probability_format(fission_rate_file,count_rates,nr_cols)
    use error_msg_class, only: raise_exception
    use benam_class,     only: benam
    use file_handling_class
    use format_class
    implicit none
    character(len=*), intent(in)  :: fission_rate_file !< File containing fission rates
    integer, intent(out)          :: count_rates       !< Amount of fission rates
    integer, intent(out)          :: nr_cols           !< Number of columns in the file
    integer,parameter               :: maxcols = 30!< Maximum number of columns
    integer                         :: j           !< Loop variable
    integer                         :: read_stat   !< Read status
    integer                         :: fissionlib  !< File id of fisionrates
    character(5)                    :: parent      !< Participating nuclei names
    integer                         :: idx_parent  !< Participating nuclei indices
    integer                         :: ifiss       !< reaction rate counter
    real(r_kind), dimension(maxcols):: buff        !< Reaclib fitting parameter
    character(500)                  :: dummy       !< Reaclib source label
    character(4)                    :: src         !< Reaclib source label
    integer                         :: fiss_probs  !< Number of fission probabilities not equal zero

    ! Counter for fission rates
    ifiss = 0

    ! Open fission file
    fissionlib = open_infile(fission_rate_file)

    ! Count the rates
    fission_prob: do
        buff(:) = -1d0

        ! An example entry could be: pa250        pa05
        !                            0.010
        read(fissionlib, *   , iostat = read_stat)  parent, src

        if (read_stat /= 0) exit

        ! Count columns, maximum allowed ones are 30
        if (ifiss .eq. 0) then
            read(fissionlib,'(A)', iostat = read_stat)  dummy
            do j=1,maxcols
                read(dummy,*,iostat=read_stat) buff(1:j)
                if (read_stat .ne. 0) exit
            end do
            ! Count the number of columns
            nr_cols = count(buff .ne. -1d0)
        else
            read(fissionlib,*, iostat = read_stat)  buff(1:nr_cols)
            if (read_stat .ne. 0) call raise_exception('Reading the fission probabilities failed.',&
                                                       "count_fission_rates_probability_format",190012)
        end if

        ! Count the number of fission probabilities not equal zero
        fiss_probs = count((buff .ne. -1d0) .and. (buff .ne. 0d0))

        ! Check if the parent is included in the network
        idx_parent = benam(parent)

        if (idx_parent .eq. 0) cycle fission_prob

        ! Count the amount of fission rates that have a
        ! probability not equal zero
        ifiss = ifiss + fiss_probs
    end do fission_prob

    ! Close the fission file again
    close(fissionlib)

    ! Store the amount of fission rates
    count_rates = ifiss
   end subroutine count_fission_rates_probability_format


   !> Initializes the rates with fragments
   !!
   !! This routine initializes the fission rates with the specified fragments.
   !! The fragments can differ for the types of fission.
   !!
   !!
   !!
   !!
   !! \b Edited:
   !!     - 31.05.2024 (MR): Created this subroutine from code parts contained in other routines
   !! .
   subroutine add_fission_fragments()
    use parameter_class,     only: fissflag, &
                                   nfission_file, bfission_file, sfission_file, &
                                   fission_frag_n_induced, fission_frag_beta_delayed,&
                                   fission_frag_spontaneous, fission_frag_missing, &
                                   fission_format_beta_delayed, fission_format_n_induced, &
                                   fission_format_spontaneous
    use global_class,        only: isotope
    use error_msg_class,     only: raise_exception, int_to_str
    use file_handling_class, only: open_infile, close_io_file
    implicit none
    integer                    :: neutfission=0 !< File IDs for abla fiss. distr. file
    integer                    :: betafission=0 !< File IDs for abla fiss. distr. file
    integer                    :: fissindex     !< Index for current fission rate
    integer                    :: astat         !< Allocation status
    real(r_kind)               :: q             !< Reaclib Q-Value
    integer                    :: nfrac_distr   !< Number of fragment distributions in file


    !-- Open the abla files, containing probabilities of the fiss. fragments
    if (fissflag .eq. 3) then
        call read_fiss_fragment_file(nfission_file, fissfrags_n_induced)
        ! Create the same for beta-delayed fission without explicitely reading the file again
        nfrac_distr = size(fissfrags_n_induced)
        allocate(fissfrags_beta_delayed(nfrac_distr),stat=astat)
        if (astat /= 0) call raise_exception('Allocation of "fissfrags_bdel" failed.',&
                                             "add_fission_fragments",190001)
        fissfrags_beta_delayed(:) = fissfrags_n_induced(:)
     else if (fissflag .eq. 4) then
        ! Read the fission fragment files if necessary
        ! n_induced
        if ((fission_frag_n_induced .eq. 3) .and. (fission_format_n_induced .ne. 0)) then
            call read_fiss_fragment_file(nfission_file, fissfrags_n_induced)
        end if
        ! beta-delayed
        if ((fission_frag_beta_delayed .eq. 3)  .and. (fission_format_beta_delayed .ne. 0)) then
            call read_fiss_fragment_file(bfission_file, fissfrags_beta_delayed)
        end if
        ! spontaneous
        if ((fission_frag_spontaneous .eq. 3)  .and. (fission_format_spontaneous .ne. 0)) then
            call read_fiss_fragment_file(sfission_file, fissfrags_spontaneous)
        end if
     else if (fissflag .eq. 5) then
        neutfission= open_infile(nfission_file)
        betafission= open_infile(bfission_file)
     end if

     ! Loop through reactions and add fission fragments
     do fissindex=1,nfiss

        if (fissrate(fissindex)%reac_type.eq.rrt_nf) then      ! neutron-induced fission
            select case (fissflag)
            case(1)
               call fiss_dist(fissrate(fissindex))
            case(2)
               call kodtakdist(fissrate(fissindex))
            case(3)
               call file_fiss_frag(fissrate(fissindex),2)
            case(4)
                select case (fission_frag_n_induced)
                case(1)
                    call fiss_dist(fissrate(fissindex))
                case(2)
                    call kodtakdist(fissrate(fissindex))
                case(3)
                    call file_fiss_frag(fissrate(fissindex),fission_frag_missing)
                case default
                    call raise_exception("Fission fragment flag ("//trim(adjustl(int_to_str(fission_frag_n_induced)))&
                                         //") not implemented yet. "//&
                                         "Set it to a supported value!",&
                                         "add_fission_fragments",&
                                         190003)
                end select
            case(5)
               call abla_nfiss(fissindex,isotope(fissrate(fissindex)%fissnuc_index)%mass,  &
                    isotope(fissrate(fissindex)%fissnuc_index)%p_nr,neutfission,q)
             case default
                call raise_exception("Fission flag ("//trim(adjustl(int_to_str(fissflag)))&
                                     //") not implemented yet. "//&
                                     "Set it to a supported value!",&
                                     "add_fission_fragments",&
                                     190003)
            end select
        else if (fissrate(fissindex)%reac_type.eq.rrt_sf) then  ! spontaneous fission
            select case (fissflag)
            case(1)
               call fiss_dist(fissrate(fissindex))
            case(2)
               call kodtakdist(fissrate(fissindex))
            case(3)
               call kodtakdist(fissrate(fissindex))
            case(4)
                select case (fission_frag_spontaneous)
                case(1)
                    call fiss_dist(fissrate(fissindex))
                case(2)
                    call kodtakdist(fissrate(fissindex))
                case(3)
                    call file_fiss_frag(fissrate(fissindex),fission_frag_missing)
                case default
                    call raise_exception("Fission fragment flag ("//trim(adjustl(int_to_str(fission_frag_spontaneous)))&
                                         //") not implemented yet. "//&
                                         "Set it to a supported value!",&
                                         "add_fission_fragments",&
                                         190003)
                end select
            case(5)
               call abla_betafiss(fissindex,isotope(fissrate(fissindex)%fissnuc_index)%mass,  &
                    isotope(fissrate(fissindex)%fissnuc_index)%p_nr,betafission,q)
            case default
                call raise_exception("Fission flag ("//trim(adjustl(int_to_str(fissflag)))&
                                        //") not implemented yet. "//&
                                        "Set it to a supported value!",&
                                        "add_fission_fragments",&
                                        190003)
            end select
        else if (fissrate(fissindex)%reac_type.eq.rrt_bf) then  ! beta-delayed fission
            select case (fissflag)
            case(1)
               call fiss_dist(fissrate(fissindex))
            case(2)
               call kodtakdist(fissrate(fissindex))
            case(3)
               call file_fiss_frag(fissrate(fissindex),2)
            case(4)
                select case (fission_frag_beta_delayed)
                case(1)
                    call fiss_dist(fissrate(fissindex))
                case(2)
                    call kodtakdist(fissrate(fissindex))
                case(3)
                    call file_fiss_frag(fissrate(fissindex),fission_frag_missing)
                case default
                    call raise_exception("Fission fragment flag ("//trim(adjustl(int_to_str(fission_frag_beta_delayed)))&
                                         //") not implemented yet. "//&
                                         "Set it to a supported value!",&
                                         "add_fission_fragments",&
                                         190003)
                end select
            case(5)
               call abla_betafiss(fissindex,isotope(fissrate(fissindex)%fissnuc_index)%mass,  &
                    isotope(fissrate(fissindex)%fissnuc_index)%p_nr,betafission,q)
            case default
                call raise_exception("Fission flag ("//trim(adjustl(int_to_str(fissflag)))&
                                        //") not implemented yet. "//&
                                        "Set it to a supported value!",&
                                        "add_fission_fragments",&
                                        190003)
            end select
        end if


     end do

    ! Clean up
    ! Close the neutron and betafission files again
    if (neutfission.gt.0) call close_io_file(neutfission,nfission_file)
    if (betafission.gt.0) call close_io_file(betafission,bfission_file)
    ! Deallocate the fission fragment arrays in case they were read from file
    if (allocated(fissfrags_beta_delayed)) then
       deallocate(fissfrags_beta_delayed,stat=astat)
       if (astat.ne.0) call raise_exception("Could not deallocate fissfrags_beta_delayed array.",&
                                            "add_fission_fragments",190002)
    end if
    if (allocated(fissfrags_n_induced)) then
       deallocate(fissfrags_n_induced,stat=astat)
       if (astat.ne.0) call raise_exception("Could not deallocate fissfrags_n_induced array.",&
                                            "add_fission_fragments",190002)
    end if
    if (allocated(fissfrags_spontaneous)) then
       deallocate(fissfrags_spontaneous,stat=astat)
       if (astat.ne.0) call raise_exception("Could not deallocate fissfrags_spontaneous array.",&
                                            "add_fission_fragments",190002)
    end if
   end subroutine add_fission_fragments




   !> Read the fission rates, splitted into different types of fission
   !!
   !! This routine serves as interface between the reading routines and
   !! the different fission formats and types. In case a new type of
   !! fission will be implemented, the reaction rates should be read here.
   !! The different formats for the different types of fission reactions are:
   !! <table>
   !! <caption id="multi_row">Fission formats</caption>
   !! <tr><th> Type                <th> Value <th> Description
   !! <tr><td> Spontaneous         <td> 0     <td> No rates are read
   !! <tr><td> Spontaneous         <td> 1     <td> Reaclib format
   !! <tr><td> Spontaneous         <td> 2     <td> Half life format
   !! <tr><td> n-induced           <td> 0     <td> No rates are read
   !! <tr><td> n-induced           <td> 1     <td> Reaclib format
   !! <tr><td> \f$\beta\f$-delayed <td> 0     <td> No rates are read
   !! <tr><td> \f$\beta\f$-delayed <td> 1     <td> Reaclib format
   !! <tr><td> \f$\beta\f$-delayed <td> 2     <td> Half life format
   !! <tr><td> \f$\beta\f$-delayed <td> 3     <td> Probability format
   !! </table>
   !!
   !! @author M. Reichert
   !! @date 31.05.2024
   subroutine read_fission_rates()
    use parameter_class, only: fission_rates_spontaneous, fission_rates_n_induced, &
                               fission_rates_beta_delayed, fission_format_spontaneous,&
                               fission_format_beta_delayed, fission_format_n_induced
    use error_msg_class, only: raise_exception, int_to_str
    implicit none
    integer :: fisscount !< Amount of fission rates

    ! Keep track of fission rates
    fisscount = 1

    ! Read the spontaneous fission rates
    select case (fission_format_spontaneous)
    case(0)
        ! Dont read anything
        continue
    case(1)
        ! Read the spontaneous fission rates in reaclib format
        call read_fission_rates_reaclib_format(fission_rates_spontaneous,rrt_sf,fisscount)
    case(2)
        ! Read the spontaneous fission rates in halflife format
        call read_fission_rates_halflife_format(fission_rates_spontaneous,rrt_sf,fisscount)
    case default
        call raise_exception("Fission format for spontaneous fission rates not implemented, got '"//&
                              int_to_str(fission_format_spontaneous)//"'. Change the "//&
                              "'fission_format_spontaneous' parameter.",&
                              "read_fission_rates",190010)
    end select

    ! Read the neutron-induced fission rates
    select case (fission_format_n_induced)
    case(0)
        ! Dont read anything
        continue
    case(1)
        call read_fission_rates_reaclib_format(fission_rates_n_induced,rrt_nf,fisscount)
    case default
        call raise_exception("Fission format for neutron-induced fission rates not implemented, got '"//&
                              int_to_str(fission_format_n_induced)//"'. Change the "//&
                              "'fission_format_n_induced' parameter.",&
                              "read_fission_rates",190010)
    end select

    ! Read the beta-delayed fission rates
    select case (fission_format_beta_delayed)
    case(0)
        ! Dont read anything
        continue
    case(1)
        ! Read the beta-delayed fission rates in reaclib format
        call read_fission_rates_reaclib_format(fission_rates_beta_delayed,rrt_bf,fisscount)
    case(2)
        ! Read the beta-delayed fission rates in halflife format
        call read_fission_rates_halflife_format(fission_rates_beta_delayed,rrt_bf,fisscount)
    case(3)
        ! Read the beta-delayed fission rates in probability format
        call read_fission_rates_probability_format(fission_rates_beta_delayed,rrt_bf,fisscount)
    case default
        call raise_exception("Fission format for beta-delayed fission rates not implemented, got '"//&
                              int_to_str(fission_format_beta_delayed)//"'. Change the "//&
                              "'fission_format_beta_delayed' parameter.",&
                              "read_fission_rates",190010)
    end select

   end subroutine read_fission_rates


   !> Reads fission rates in Reaclib format
   !!
   !! This subroutine reads the fission rates that are specified in the file
   !! parameter_class::fission_rates. Fission rates are implemented with the same
   !! fits and a similar format as usual reaclib rates (see reaclib_rate_module::read_reaclib).
   !! An example looks like:
   !! \file{
   !!    th256                                 ms99w     0.00000E+00
   !! 1.988512E-02 0.000000E+00 0.000000E+00 0.000000E+00
   !! 0.000000E+00 0.000000E+00 0.000000E+00
   !!    th257                                 ms99w     0.00000E+00
   !! -1.483947E+00 0.000000E+00 0.000000E+00 0.000000E+00
   !! 0.000000E+00 0.000000E+00 0.000000E+00
   !! ...}
   !! The chapter structure is the same as in reaclib, but only the parent nucleus
   !! has an entry. The products are given by the respective fission fragment
   !! distribution.
   !!
   !! @returns \ref fissrate filled with fission rates
   !! @author M. Eichler
   !!
   !! \b Edited:
   !!   - 25.01.21 (MR): Moved from init_network to this new independent subroutine
   !!   - 27.05.24 (MR): Moved initializing fragments to independent subroutine
   !! -
   subroutine read_fission_rates_reaclib_format(fission_path,reac_type,start_idx)
      use parameter_class, only: fissflag
      use benam_class,     only: benam
      use global_class,    only: isotope, ineu
      use file_handling_class
      use format_class
      use error_msg_class
      implicit none
      character(len=*), intent(in) :: fission_path  !< Path to fission rates file
      integer, intent(in)          :: reac_type     !< Reaction type
      integer, intent(inout)       :: start_idx     !< Start index for fission rates
      integer                      :: read_stat     !< Read status
      integer                      :: fissionlib    !< File id of fisionrates
      character(5), dimension(6)   :: parts         !< Participating nuclei names
      integer, dimension(6)        :: parts_index   !< Participating nuclei indices
      real(r_kind), dimension(9)   :: params        !< Reaclib fitting parameter
      character(4)                 :: src           !< Reaclib source label
      character(1)                 :: res           !< Reaclib weak flag label
      character(1)                 :: rev           !< Reaclib reverse label
      real(r_kind)                 :: q             !< Reaclib Q-Value
      integer                      :: grp           !< Reaclib chapter
      integer                      :: group_index   !< Storage for the current chapter
      integer                      :: fissindex     !< Index for current fission rate
      integer                      :: j             !< Loop variable



      !-- Open the fission library file
      fissionlib = open_infile(fission_path)

      ! make sure counting index for reaction rates continues for fission rates
      fissindex = start_idx
      fission_loop: do
         read(fissionlib,my_format(1), iostat = read_stat)  &
              grp, parts(1:6), src, res, rev, q
         if (read_stat /= 0) exit fission_loop
         read(fissionlib,"(4E13.6)") params(1:4)
         read(fissionlib,"(5E13.6)") params(5:9)

         select case (grp)
         case (1:11)
            group_index = grp
            cycle fission_loop
         case default
            parts_index = 0
            inner_fission: do j=1,6
               if (parts(j) .eq. '     ') exit inner_fission
               parts_index(j) = benam(parts(j))
               ! parts_index(i) = 0 if nuclide i is not in network
               if (parts_index(j) .eq. 0) cycle fission_loop
            end do inner_fission
         end select

         ! Keep track for debugging
         if (reac_type .eq. rrt_nf) n_nf=n_nf+1 !< count neutron induced fission
         if (reac_type .eq. rrt_sf) n_sf=n_sf+1 !< count spontaneous fission
         if (reac_type .eq. rrt_bf) n_bf=n_bf+1 !< count beta delayed fission

         ! Set the amount of released neutrons per default to 0
         fissrate(fissindex)%released_n = 0

         if (reac_type.eq.rrt_nf) then      ! neutron-induced fission
            fissrate(fissindex)%reac_type = reac_type
            fissrate(fissindex)%param = params
            ! Save fission nucleus index
            if (parts_index(1) .eq. ineu) then
               fissrate(fissindex)%fissnuc_index = parts_index(2)
            else
               fissrate(fissindex)%fissnuc_index = parts_index(1)
            end if

         else if (reac_type.eq.rrt_sf) then  ! spontaneous fission
            fissrate(fissindex)%reac_type = reac_type
            fissrate(fissindex)%param = params
            ! Save fission nucleus index
            fissrate(fissindex)%fissnuc_index = parts_index(1)

         else if (reac_type.eq.rrt_bf) then  ! beta-delayed fission
            ! Check if there should be neutrons released
            fissrate(fissindex)%released_n = count(parts_index(:) .eq. ineu)
            fissrate(fissindex)%reac_type = reac_type
            fissrate(fissindex)%param = params
            ! Save fission nucleus index
            fissrate(fissindex)%fissnuc_index = parts_index(1)


         end if
         fissrate(fissindex)%src      = src


         fissindex = fissindex + 1
      end do fission_loop

      ! Return the number of fission rates that have been added
      start_idx = fissindex

      ! Close the fission file again
      close(fissionlib)

   end subroutine read_fission_rates_reaclib_format



   !> Read the fission rates in probability format
   !!
   !! This subroutine reads fission rates in probability format
   !! and stores the result in count_rates. The subroutine can be called
   !! for beta-delayed fission and will throw an error for other types of fission.
   !! An example of an entry can be:
   !! \file{
   !! pa291        mp18
   !! 0.000    0.008    0.000    0.015    0.000    0.018    0.000    0.025    0.000    0.003
   !! pa292        mp18
   !! 0.000    0.001    0.002    0.000    0.004    0.000    0.006    0.000    0.018    0.000
   !! pa293        mp18
   !! 0.239    0.001    0.000    0.009    0.000    0.010    0.000    0.014    0.000    0.009
   !! ...}
   !! where the first line of an entry gives nucleus name and the source. The second line
   !! gives the probabilities of the fission channels (i.e., P0nf, P1nf, ...).
   !! Note that there is a maximum number of columns (i.e., neutron emission) set to 30 here.
   !! Throughout the file, the amount of possible channels have to stay the same, but they
   !! can be different for two different files. For example, for Panov et al. 2005 the file
   !! could look like:
   !! \file{
   !!  u259        pa05
   !! 0.990
   !!  u260        pa05
   !! 0.020
   !!  u261        pa05
   !! 0.990
   !! ...}.
   !!
   !! @author M. Reichert
   !! @date 31.05.24
   subroutine read_fission_rates_probability_format(fission_path,reac_type,start_idx)
    use benam_class,     only: benam
    use global_class,    only: isotope, ineu
    use error_msg_class, only: raise_exception
    use file_handling_class
    use format_class
    implicit none
    character(len=*), intent(in)       :: fission_path!< Path to fission rates file
    integer, intent(in)                :: reac_type   !< Reaction type
    integer, intent(inout)             :: start_idx   !< Start index for fission rates
    ! Internal variables
    integer                            :: read_stat   !< Read status
    integer                            :: fissionlib  !< File id of fisionrates
    character(5), dimension(6)         :: parts       !< Participating nuclei names
    integer, dimension(6)              :: parts_index !< Participating nuclei indices
    real(r_kind), dimension(9)         :: params      !< Reaclib fitting parameter
    character(4)                       :: src         !< Reaclib source label
    integer                            :: group_index !< Storage for the current chapter
    integer                            :: fissindex   !< Index for current fission rate
    integer                            :: j           !< Loop variable
    integer                            :: idx_parent  !< Index of the parent nucleus
    character(5)                       :: parent      !< Parent nucleus name
    real(r_kind),dimension(amount_cols):: Pnf         !< Probability of beta-delayed neutron emission fission

    ! Counter for fission rates
    fissindex = start_idx

    ! Open fission file
    fissionlib = open_infile(fission_path)

    ! Count the rates
    fission_prob: do
        ! An example entry could be: pa250        pa05
        !                            0.010
        read(fissionlib,*, iostat = read_stat)  parent, src
        if (read_stat /= 0) exit
        read(fissionlib,*, iostat = read_stat)  Pnf
        ! Throw error if the file ends with a parent and a source
        if (read_stat /= 0) call raise_exception("Could not read fission probabilities.",&
                                                 "read_fission_rates_probability_format",&
                                                 190012)


        ! Check if the parent is included in the network
        idx_parent = benam(parent)
        if (idx_parent .eq. 0) cycle fission_prob


        parts_index(:) = 0
        parts_index(1) = idx_parent

        do j=1,amount_cols
            if (Pnf(j) .eq. 0d0) cycle

            params(:) = 0d0
            params(1) = Pnf(j) ! This will be changed once the reactions are merged
            params(2) = j-1    ! Gives the amount of released neutrons

            if (reac_type.eq.rrt_bf) then  ! beta-delayed fission
                group_index = 2
                fissrate(fissindex)%reac_type = reac_type
                fissrate(fissindex)%param = params
                ! Save fission nucleus index
                fissrate(fissindex)%fissnuc_index = parts_index(1)
                fissrate(fissindex)%released_n = j-1
            else
                call raise_exception("Fission type not implemented for probability format!",&
                                     "read_fission_rates_probability_format",&
                                     190011)
            end if
            ! Put also the source there
            fissrate(fissindex)%src   = src

            ! Also save it in a separate array that is more easier to
            ! access by index
            beta_delayed_fiss_probs(idx_parent,j) = Pnf(j)

            ! Keep track for debugging
            if (reac_type .eq. rrt_bf) n_bf=n_bf+1 ! count beta delayed fission

            fissindex = fissindex + 1
        end do

    end do fission_prob

    ! Close the fission file again
    close(fissionlib)

   end subroutine read_fission_rates_probability_format


   !> Read fission rates in half-life format
   !!
   !! This subroutine reads fission rates in half-life format.
   !! The subroutine can be called for spontaneous and beta-delayed fission
   !! and throws an error otherwise.
   !! An example of an entry can be:
   !! \file{
   !!  th198    4.760787e-09    kh20
   !!  th199    2.808520e-08    kh20
   !!  th200    2.841403e-08    kh20
   !! ...}
   !! where the first columns gives the nucleus, the second the halflife in seconds
   !! and the third the source. Note that in this format, no beta-delayed neutrons
   !! can be specified.
   !!
   !! @author M. Reichert
   !! @date 30.05.24
   subroutine read_fission_rates_halflife_format(fission_path,reac_type,start_idx)
    use benam_class,     only: benam
    use global_class,    only: isotope, ineu
    use error_msg_class, only: raise_exception
    use file_handling_class
    use format_class
    implicit none
    character(len=*), intent(in) :: fission_path  !< Path to fission rates file
    integer, intent(in)          :: reac_type     !< Reaction type
    integer, intent(inout)       :: start_idx     !< Start index for fission rates
    integer                      :: read_stat     !< Read status
    integer                      :: fissionlib    !< File id of fisionrates
    character(5), dimension(6)   :: parts         !< Participating nuclei names
    integer, dimension(6)        :: parts_index   !< Participating nuclei indices
    real(r_kind), dimension(9)   :: params        !< Reaclib fitting parameter
    character(4)                 :: src           !< Reaclib source label
    integer                      :: group_index   !< Storage for the current chapter
    integer                      :: fissindex     !< Index for current fission rate
    integer                      :: j             !< Loop variable



    !-- Open the fission library file
    fissionlib = open_infile(fission_path)

    ! make sure counting index for reaction rates continues for fission rates
    fissindex = start_idx
    fission_loop: do
       ! Set parts to empty strings
       parts(1:6) = '     '
       params(:)  = 0d0
       read(fissionlib, *, iostat = read_stat)  &
           parts(1), params(1), src
       if (read_stat /= 0) exit fission_loop

       params(1) = dlog(dlog(2d0)/params(1)) ! Convert half-life to reaclib a0
       parts_index = 0
       inner_fission: do j=1,6
           if (parts(j) .eq. '     ') exit inner_fission
           parts_index(j) = benam(parts(j))
           ! parts_index(i) = 0 if nuclide i is not in network
           if (parts_index(j) .eq. 0) cycle fission_loop
       end do inner_fission

       ! Keep track for debugging
       if (reac_type .eq. rrt_nf) n_nf=n_nf+1 ! count neutron induced fission
       if (reac_type .eq. rrt_sf) n_sf=n_sf+1 ! count spontaneous fission
       if (reac_type .eq. rrt_bf) n_bf=n_bf+1 ! count beta delayed fission

       ! Per default set number of released neutrons to zero (important for bdel fission)
       fissrate(fissindex)%released_n = 0
      if (reac_type .eq.rrt_sf) then  ! spontaneous fission
          fissrate(fissindex)%reac_type = reac_type
          fissrate(fissindex)%param = params
          ! Save fission nucleus index
          fissrate(fissindex)%fissnuc_index = parts_index(1)

       else if (reac_type .eq. rrt_bf) then  ! beta-delayed fission
          fissrate(fissindex)%reac_type = reac_type
          fissrate(fissindex)%param = params
          ! Save fission nucleus index
          fissrate(fissindex)%fissnuc_index = parts_index(1)
       else
          call raise_exception("Fission type not implemented for half life format!",&
                               "read_fission_rates_halflife_format",&
                               190011)
       end if

       ! Set some reaclib chapter, note that n induced fission should
       ! not have a halflife format
       group_index = 2

       ! Put the source to the rate
       fissrate(fissindex)%src        = src

       fissindex = fissindex + 1
    end do fission_loop

    ! Return the number of fission rates that have been added
    start_idx = fissindex

    ! Close the fission file again
    close(fissionlib)

   end subroutine read_fission_rates_halflife_format






   !> Determines fission fragment mass distribution as described in Panov et al. 2001.
   !!
   !! This routine is called for \ref parameter_class::fissflag = 1 and fills
   !! the array \ref fissrates with indices of the fragment nuclei.
   !!
   !! @see [Panov et al., Nuc. Phys. A688 2001](https://ui.adsabs.harvard.edu/abs/2001NuPhA.688..587P/abstract)
   !! @author C. Winteler
   subroutine fiss_dist(fissrate_inout)
      use global_class,    only: isotope, ineu
      use benam_class,     only: minmax, findaz
      use error_msg_class, only: raise_exception, int_to_str
      implicit none
      type(fissionrate_type),intent(inout)  :: fissrate_inout  !< Temporary fission rate
      ! Internal variables
      integer                       :: mass      !< mass of fissioning nucleus
      integer                       :: pnr       !< proton number of fissioning nucleus
      integer                       :: af,zf     !< mass and proton number of "compound" nucleus
      integer                       :: a1,z1     !< mass and proton number of fragment 1
      integer                       :: a2,z2     !< mass and proton number of fragment 2
      integer                       :: nemiss    !< number of fission neutrons
      integer                       :: fiss_mode !< specifies fission mode
      real(r_kind)                  :: qval      !< Qvalue of fission reaction

      INFO_ENTRY("fiss_dist")

      ! Counter for additional neutrons in case of
      ! fragment beyond dripline.
      nemiss = 0

      ! Get the mass and proton number
      mass = isotope(fissrate_inout%fissnuc_index)%mass
      pnr  = isotope(fissrate_inout%fissnuc_index)%p_nr

      !neutron induced fission
      if (fissrate_inout%reac_type .eq. rrt_nf) then
         fiss_mode = 1
         af = mass+1
         zf = pnr
      !spontaneous fission
      else if (fissrate_inout%reac_type .eq. rrt_sf) then
         fiss_mode = 2
         af = mass
         zf = pnr
        !beta delayed fission
      else if (fissrate_inout%reac_type .eq. rrt_bf) then
         fiss_mode = 3
         ! Take care of beta-delayed neutron emission fission
         ! The amount of neutrons is stored in the second parameter
         ! of the fission rate
         nemiss = fissrate_inout%released_n
         af = mass-nemiss
         zf = pnr + 1
      end if

      select case (af)
      case(255:265) ! -> symmetric fission
         a2 = nint(af/2.d0)
         z2 = nint(zf/2.d0)
         a1 = af - a2
         z1 = zf - z2
      case default ! -> asymmetric fission
         a2 = 130
         z2 = nint(52.01 - (zf - 80)/1.d1)
   !       a2 = max(a2,af-a2)     !modification to prevent very n-rich
   !       z2 = max(z2,zf-z2)     !fragments
         a1 = af - a2
         z1 = zf - z2
      end select

      if (a1.gt.minmax(z1,2)) then
         nemiss = a1 - minmax(z1,2)
         a1 = minmax(z1,2)
      else if (a1.lt.minmax(z1,1)) then
         if (VERBOSE_LEVEL .ge. 2) then
            print *, 'a1 error:', z1, a1, zf, af
            print *, 'adjusting...'
         end if
         a1 = minmax(z1,1)
         a2 = af - a1
      end if

      if (a2.gt.minmax(z2,2)) then
         nemiss = nemiss + a2 - minmax(z2,2)
         a2 = minmax(z2,2)

      else if (a2.lt.minmax(z2,1)) then
         if (VERBOSE_LEVEL .ge. 2) then
            print *, 'a2 error:', z2, a2, zf, af
            print *, 'adjusting...'
         end if
         a2 = minmax(z2,1)
         a1 = af - a2
      end if
      nufiss = nufiss + nemiss

      fissrate_inout%fissnuc_index  = findaz(mass,pnr)
      fissrate_inout%channels       = 1
      ! TODO THROW ERROR HERE!!
      allocate(fissrate_inout%channelprob(1))
      allocate(fissrate_inout%q_value(1))
      fissrate_inout%channelprob(1) = 1.d0
      fissrate_inout%mode           = fiss_mode

      if (fiss_mode.eq.1) then
        fissrate_inout%dimens = 4
         allocate(fissrate_inout%fissparts(fissrate_inout%dimens))
         allocate(fissrate_inout%ch_amount(fissrate_inout%dimens))
         fissrate_inout%fissparts(1) = ineu
         fissrate_inout%ch_amount(1) = float(nemiss) - 1.d0
         fissrate_inout%fissparts(2) = fissrate_inout%fissnuc_index
         fissrate_inout%ch_amount(2) = -1.d0
         fissrate_inout%fissparts(3) = findaz(a1,z1)
         fissrate_inout%ch_amount(3) = 1.d0
         fissrate_inout%fissparts(4) = findaz(a2,z2)
         fissrate_inout%ch_amount(4) = 1.d0
         fissrate_inout%q_value(1)   = (1-nemiss) * isotope(ineu)%mass_exc + &
                                      isotope(fissrate_inout%fissnuc_index)%mass_exc - &
                                      isotope(fissrate_inout%fissparts(3))%mass_exc - &
                                      isotope(fissrate_inout%fissparts(4))%mass_exc
      else                                  ! fiss_mode 2 and 3
         if (nemiss .eq. 0) then            ! no fission neutrons
            fissrate_inout%dimens = 3
            allocate(fissrate_inout%fissparts(fissrate_inout%dimens))
            allocate(fissrate_inout%ch_amount(fissrate_inout%dimens))
         else                               ! fission neutrons are produced
            fissrate_inout%dimens = 4
            allocate(fissrate_inout%fissparts(fissrate_inout%dimens))
            allocate(fissrate_inout%ch_amount(fissrate_inout%dimens))
            fissrate_inout%fissparts(4) = ineu
            fissrate_inout%ch_amount(4) = float(nemiss)
         end if
         fissrate_inout%fissparts(1) = fissrate_inout%fissnuc_index
         fissrate_inout%ch_amount(1) = -1.d0
         fissrate_inout%fissparts(2) = findaz(a1,z1)
         fissrate_inout%ch_amount(2) = 1.d0
         fissrate_inout%fissparts(3) = findaz(a2,z2)
         fissrate_inout%ch_amount(3) = 1.d0

         if (fissrate_inout%fissparts(2) .eq. -1) then
            call raise_exception("Fragment was not found in nucleus library, does your library have 'holes'? "//&
                                 "(A="//int_to_str(a1)//", Z="//int_to_str(z1)//")",&
                                 "fiss_dist",190013)
         end if
         if (fissrate_inout%fissparts(3) .eq. -1) then
            call raise_exception("Fragment was not found in nucleus library, does your library have 'holes'? "//&
                                 "(A="//int_to_str(a2)//", Z="//int_to_str(z2)//")",&
                                 "fiss_dist",190013)
         end if

         fissrate_inout%q_value(1)   = isotope(fissrate_inout%fissparts(1))%mass_exc - &
                                      isotope(fissrate_inout%fissparts(2))%mass_exc - &
                                      isotope(fissrate_inout%fissparts(3))%mass_exc - &
                                      nemiss * isotope(ineu)%mass_exc
      end if

      qval = fissrate_inout%q_value(1)
      fissrate_inout%averageQ = qval
      INFO_EXIT("fiss_dist")

      return

   end subroutine fiss_dist



   !>
   !! Kodama-Takahashi distribution
   !!
   !! Subroutine that calculates fission fragment distribution as
   !! described in Kodama & Takahashi, Nuc. Phys. Sec.A, Vol. 239,
   !! Issue 3, p. 489-510, 1975.
   !!
   !! @see [Kodama & Takahashi 1975](https://ui.adsabs.harvard.edu/abs/1975NuPhA.239..489K/abstract)
   !!
   !! \b Edited:
   !!      - 11.01.14
   !!      - 05.12.17 (ME)
   !!      - 30.05.2024 (MR): Fixed bug related to neutron emission if fragment was not in Sunet
   !! .
   !! @author C. Winteler
   subroutine kodtakdist(fissrate_inout)
      use parameter_class, only: unit
      use global_class,    only: isotope,ineu
      use benam_class,     only: minmax,findaz
      use error_msg_class, only: raise_exception, int_to_str
      implicit none
      type(fissionrate_type),intent(inout)  :: fissrate_inout  !< Temporary fission rate
      ! Internal variables
      integer                 :: mass             !< Mass number of parent nucleus
      integer                 :: pnr              !< Proton number of parent nucleus
      integer                 :: nf               !< number of fission channels
      integer                 :: fiss_mode        !< Fission mode (1: n-induced, 2: spontaneous, 3: beta-delayed fission)
      integer                 :: af,zf            !< mass and proton number of fissioning nucleus
      integer                 :: a,z              !< mass and proton number of fragment
      integer                 :: a1,a2,z1,z2
      real(r_kind)            :: paz,ptot
      real(r_kind)            :: za,al,ah,cz,ca
      real(r_kind),parameter  :: lim = 1.d-6         !< Limit for fragments to be taken into account
      integer                 :: nemiss              !< Neutrons missing in case of fragment beyond dripline for one fragment
      integer                 :: alloc_stat          !< Allocation status
      integer                 :: i                   !< Loop variable
      integer                 :: dimens              !< Dimension of the parts in fission rate (amount fragments + 1)
      integer                 :: neutronflag         !< Flag to check whether any fragment is beyond dripline
      real(r_kind)            :: nemiss_total        !< Total number of neutrons emitted (sum of nemiss)
      integer                 :: additional_neutrons !< Neutrons that may be emitted by beta decay
      real(r_kind)            :: qval                !< Q-value of the fission reaction

      INFO_ENTRY("kodtakdist")

      ! Set neutron flag to zero
      neutronflag = 0
      additional_neutrons = 0d0

      ! Get the mass and proton number
      mass = isotope(fissrate_inout%fissnuc_index)%mass
      pnr  = isotope(fissrate_inout%fissnuc_index)%p_nr

      !neutron induced fission
      if (fissrate_inout%reac_type .eq. rrt_nf) then
        fiss_mode = 1
        af = mass+1
        zf = pnr
      !spontaneous fission
      else if (fissrate_inout%reac_type .eq. rrt_sf) then
        fiss_mode = 2
        af = mass
        zf = pnr
      !beta delayed fission
      else if (fissrate_inout%reac_type .eq. rrt_bf) then
        fiss_mode = 3
        ! Take care of beta-delayed neutron emission fission
        ! The amount of neutrons is stored in the second parameter
        ! of the fission rate
        additional_neutrons = fissrate_inout%released_n
        af = mass-additional_neutrons
        zf = pnr + 1
      end if

      al = 0.85d0*af - 104.98d0
      ah = 0.15d0*af + 103.87d0
      cz = 0.8d0
      ca = 78.d0

      nf = 0
      ptot = 0.d0

      do a = 1,af
         za = zf*(a+0.6d0)/af
         do z=1,zf
            paz = 0.5d0*dexp(-((z-za)**2)/cz) *                         &
            ((dexp(-((a-al)**2)/ca))+(dexp(-((a-ah)**2)/ca)))/      &
            (unit%pi*dsqrt(cz*ca))
            if(paz.ge.lim) then
                nf = nf+1
                ptot = ptot + paz
                a1 = a
                z1 = z
                a2 = af - a1
                z2 = zf - z1

               ! check if any neutrons are emitted (needed to determine fissrate()%dimens)
               if (a1.gt.minmax(z1,2)) then
                  neutronflag = 1
               else if (a1.lt.minmax(z1,1)) then
                  if (VERBOSE_LEVEL .ge. 2) then
                    print *, 'Warning in kodtakdist, fragment was lighter than included nuclei'
                    print *, 'Fragment was: ', a1, z1
                  end if
               end if

               if (a2.gt.minmax(z2,2)) then
                  neutronflag = 1
               else if (a2.lt.minmax(z2,1)) then
                  if (VERBOSE_LEVEL .ge. 2) then
                    print *, 'Warning in kodtakdist, fragment was lighter than included nuclei'
                    print *, 'Fragment was: ', a2, z2
                  end if
               end if

            end if
         end do
      end do

      fissrate_inout%channels       = nf
      fissrate_inout%fissnuc_index  = findaz(mass,pnr)
      fissrate_inout%mode           = fiss_mode

      ! Allocate and complain in case of failure
      allocate(fissrate_inout%channelprob(nf),fissrate_inout%q_value(nf), &
               stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Allocation of fission rate prob. arrays failed.",&
                                                  "kodtakdist",190001)


      select case (fiss_mode)
      case(1)
         dimens = 2 * fissrate_inout%channels + 2
         fissrate_inout%dimens = dimens
         allocate(fissrate_inout%fissparts(dimens),fissrate_inout%ch_amount(dimens),&
                  stat=alloc_stat)
         if (alloc_stat .ne. 0) call raise_exception("Allocation of fission rate parts failed.",&
                                                     "kodtakdist",190001)

         fissrate_inout%fissparts(1) = ineu
         fissrate_inout%fissparts(2) = fissrate_inout%fissnuc_index
         fissrate_inout%ch_amount(2) = -1.d0
      case(2:)
         dimens = 2 * fissrate_inout%channels + 1
         if (neutronflag.eq.1) dimens = dimens+1
         fissrate_inout%dimens = dimens
         allocate(fissrate_inout%fissparts(dimens),fissrate_inout%ch_amount(dimens),&
                  stat=alloc_stat)
         if (alloc_stat .ne. 0) call raise_exception("Allocation of fission rate parts failed.",&
                                                     "kodtakdist",190001)
         fissrate_inout%ch_amount(1) = -1.d0
         fissrate_inout%fissparts(1) = fissrate_inout%fissnuc_index
      end select

      ! Keep track of total number of neutrons emitted
      nemiss_total = float(additional_neutrons)
      !---- start writing the fission rates to rrate
      i = 0
      do z = 1,zf
         massloop: do a = 1,af
            za = zf*(a+0.6d0)/af
            paz = 0.5d0*dexp(-((z-za)**2)/cz) *                         &
                 ((dexp(-((a-al)**2)/ca))+(dexp(-((a-ah)**2)/ca)))/      &
                 (unit%pi*dsqrt(cz*ca))
            if(paz.lt.lim) cycle massloop
            i = i+1
            paz = paz/ptot  ! normalise probabilities to 1
            a1 = a
            z1 = z
            a2 = af - a1
            z2 = zf - z1

            nemiss = 0

            if (a1.gt.minmax(z1,2)) then
               nemiss = a1 - minmax(z1,2)
               a1 = minmax(z1,2)

            else if (a1.lt.minmax(z1,1)) then
               if (VERBOSE_LEVEL .ge. 2) then
                  print *, 'Warning in kodtakdist, fragment was lighter than included nuclei'
                  print *, 'Fragment was: ', a1, z1
               end if
            end if

            if (a2.gt.minmax(z2,2)) then
               nemiss = nemiss + a2 - minmax(z2,2)
               a2 = minmax(z2,2)
            else if (a2.lt.minmax(z2,1)) then
                ! TODO: Maybe borrow neutrons here?
               if (VERBOSE_LEVEL .ge. 2) then
                print *, 'Warning in kodtakdist, fragment was lighter than included nuclei'
                print *, 'Fragment was: ', a2, z2
               end if
            end if

            ! Count the total amount of neutrons emitted
            nemiss_total = nemiss_total+ float(nemiss)*paz
            ! fill in fissrate()%fissparts(); same as rrate()%parts(1:6), but with individual array sizes
            select case(fiss_mode)
            case(1)
               fissrate_inout%channelprob(i) = paz
               fissrate_inout%fissparts(i+2) = findaz(a1,z1)
               ! Check that the nucleus was really found!
               if (fissrate_inout%fissparts(i+2) .eq. -1) then
                call raise_exception("Fragment was not found in nucleus library, does your library have 'holes'? "//&
                                     "(A="//int_to_str(a1)//", Z="//int_to_str(z1)//")",&
                                     "kodtakdist",190013)
               end if

               fissrate_inout%ch_amount(i+2) = fissrate_inout%channelprob(i)   ! rate at which fragment is produced per destroyed parent nucleus
               fissrate_inout%fissparts(i+2+fissrate_inout%channels) = findaz(a2,z2)
               ! Check that the nucleus was really found!
               if (fissrate_inout%fissparts(i+2+fissrate_inout%channels) .eq. -1) then
                call raise_exception("Fragment was not found in nucleus library, does your library have 'holes'? "//&
                                     "(A="//int_to_str(a2)//", Z="//int_to_str(z2)//")",&
                                     "kodtakdist",190013)
               end if
               fissrate_inout%ch_amount(i+2+fissrate_inout%channels) = fissrate_inout%channelprob(i)
               fissrate_inout%q_value(i) = isotope(fissrate_inout%fissnuc_index)%mass_exc - &
                       (nemiss - 1) * isotope(ineu)%mass_exc - &                   ! (nemiss - 1) to account for one neutron that is destroyed
                       isotope(fissrate_inout%fissparts(i+2))%mass_exc - &
                       isotope(fissrate_inout%fissparts(i+2+fissrate_inout%channels))%mass_exc
            case(2:)
               fissrate_inout%channelprob(i) = paz
               fissrate_inout%fissparts(i+1) = findaz(a1,z1)
               ! Check that the nucleus was really found!
               if (fissrate_inout%fissparts(i+1) .eq. -1) then
                call raise_exception("Fragment was not found in nucleus library, does your library have 'holes'?",&
                                     "kodtakdist",190013)
               end if

               fissrate_inout%ch_amount(i+1) = fissrate_inout%channelprob(i)    ! rate at which fragment is produced per destroyed parent nucleus
               fissrate_inout%fissparts(i+1+fissrate_inout%channels) = findaz(a2,z2)
               ! Check that the nucleus was really found!
               if (fissrate_inout%fissparts(i+1+fissrate_inout%channels) .eq. -1) then
                call raise_exception("Fragment was not found in nucleus library, does your library have 'holes'?",&
                                     "kodtakdist",190013)
               end if

               fissrate_inout%ch_amount(i+1+fissrate_inout%channels) = fissrate_inout%channelprob(i)
               if (neutronflag.eq.1) then
                  fissrate_inout%fissparts(dimens) = ineu
                  fissrate_inout%ch_amount(dimens) = real(nemiss)
               end if
               fissrate_inout%q_value(i) = isotope(fissrate_inout%fissnuc_index)%mass_exc - &
                       nemiss * isotope(ineu)%mass_exc - &
                       isotope(fissrate_inout%fissparts(i+1))%mass_exc - &
                       isotope(fissrate_inout%fissparts(i+1+fissrate_inout%channels))%mass_exc
            end select
         end do massloop
      end do

      ! Set the neutrons correctly
      select case(fiss_mode)
      case(1)
        fissrate_inout%ch_amount(1)   = nemiss_total - 1.d0
      case(2:)
        if (neutronflag.eq.1) then
            fissrate_inout%fissparts(dimens) = ineu
            fissrate_inout%ch_amount(dimens) = nemiss_total
        end if
      end select

      qval = 0.d0
      do i=1,fissrate_inout%dimens
          qval = qval - isotope(fissrate_inout%fissparts(i))%mass_exc*fissrate_inout%ch_amount(i) ! weighted average Q-value
      end do
      fissrate_inout%averageQ = qval

      nufiss = nufiss + int(nemiss_total)

      INFO_EXIT("kodtakdist")
      return

   end subroutine kodtakdist



   !> Read the fission distribution from a file
   !!
   !! The fission distribution is read from a file and stored in the
   !! fragment_array. The default file contains the fission distribution
   !! according to [Mumpower et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020PhRvC.101e4607M/abstract).
   !! However, the format has been modified compared to the original publication.
   !! An example of the file looks like:
   !! \file{
   !! 93  220  574
   !! 0 0 0
   !!  019  048  1.96607e-06
   !!  020  048  2.24289e-05
   !!  020  049  1.88041e-05
   !!  020  050  7.65944e-06
   !! }
   !! The first line contains the Z, A and number of the fissioning nucleus.
   !! The next line is a dummy line. Then followed by Z, A and the yield Y(Z,A).
   !! Note that \f$\sum Y(Z,A) \cdot A \f$ should be the mass number of the
   !! parent nucleus. Furthermore,  \f$\sum Y(Z,A) = 2\f$ in case 2 fragments are created.
   !! A large part of the subroutine deals with the problem what to do if the fragment
   !! is not included in the network. In case of a nucleus heavier than the ones contained in
   !! the network, the fragment is split into a lighter isotope plus neutrons.
   !! For lighter nuclei where neutrons would be necessary to create a nucleus that is in the network,
   !! the neutrons are "borrowed" from the heaviest fragment. This latter case should almost never happen.
   !! If the remapping fails for some reason, an error will be raised.
   !!
   !! @see [Mumpower et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020PhRvC.101e4607M/abstract).
   !!
   !! @author M. Reichert
   !! @date 14.02.23
   subroutine read_fiss_fragment_file(file,fragment_array)
    use file_handling_class
    use benam_class,     only: findaz, minmax
    use error_msg_class, only: raise_exception
    implicit none

    character(len=*), intent(in)                           :: file           !< File path to fission fragment file
    type(fragtype), dimension(:), allocatable, intent(out) :: fragment_array !< Array of fission fragments

    integer :: file_id       !< Integer id of the file
    integer :: frag_counter  !< Count the number of distributions
    integer :: zf            !< Z fission nucleus
    integer :: af            !< A fission nucleus
    integer :: nof           !< Number of fragments
    integer :: i,j,m         !< Loop variable
    integer :: read_stat     !< Status variable
    integer :: astat         !< Allocation status
    integer :: nem           !< Number of emitted neutrons
    integer :: bn            !< Number of borrowed neutrons
    real(r_kind) :: Ynem     !< Yield of emitted neutrons
    real(r_kind) :: helper   !< Yield of emitted neutrons
    integer :: mina,maxa     !< Minimum and maximum A


    INFO_ENTRY("read_fiss_fragment_file")

    ! Open the file
    file_id =open_infile(file)

    ! Count the number of relevant
    ! fragment distributions
    frag_counter = 0
    do
        read(file_id,*,iostat=read_stat) zf, af, nof
        if (read_stat .ne. 0) exit ! end of file
        read(file_id,*)

        do i=1,nof
            read(file_id,*)
        end do

        if (((findaz(af,zf) .ne. -1) .or. &
             (findaz(af-1,zf) .ne. -1) .or. &
             (findaz(af,zf-1) .ne. -1))) then
            frag_counter = frag_counter + 1
        end if
    end do

    ! Allocate the distribtions
    allocate(fragment_array(frag_counter),stat=astat)
    if (astat .ne. 0) then
        call raise_exception("Could not allocate 'fragment_array'.", &
                             "read_fiss_fragment_file",190001)
    end if

    ! Read the file again
    rewind(file_id)

    fragloop: do i = 1,frag_counter
        read(file_id,*,iostat=read_stat) zf, af, nof
        read(file_id,*)


        ! Nucleus not included nor the decay of
        ! a nucleus in the network, nor reachable
        ! via a neutron capture
        if (.not. ((findaz(af,zf)   .ne. -1) .or. &
                   (findaz(af-1,zf) .ne. -1) .or. &
                   (findaz(af,zf-1) .ne. -1))) then
            do j=1,nof
                read(file_id,*)
            end do
            cycle fragloop
        end if


        ! Allocate the fragments
        allocate(fragment_array(i)%net_idx(nof),&
                 fragment_array(i)%Z(nof),&
                 fragment_array(i)%A(nof),&
                 fragment_array(i)%Y(nof),&
                 stat=astat)
        if (astat .ne. 0) then
            call raise_exception("Could not allocate fragment properties.", &
                                 "read_fiss_fragment_file",190001)
        end if

        ! Store the properties
        fragment_array(i)%nr_frags = nof
        fragment_array(i)%Zp = zf
        fragment_array(i)%Ap = af

        ! Count the neutrons emitted
        nem = 0
        Ynem = 0

        ! Read the fragments
        do j=1,nof
            read(file_id,*) fragment_array(i)%Z(j),fragment_array(i)%A(j),fragment_array(i)%Y(j)
            fragment_array(i)%net_idx(j) =findaz(fragment_array(i)%A(j),fragment_array(i)%Z(j))

            ! What to do when the fragment is not in the network?
            ! Ideally make neutrons until the fragment is in the network
            if (fragment_array(i)%net_idx(j) .eq. -1) then
                ! Tell the world that this is done
                if (VERBOSE_LEVEL .ge. 2) then
                    write(*,*) "Fragment ",fragment_array(i)%Z(j),fragment_array(i)%A(j),&
                               " not in network. Splitting into lighter isotope and neutrons."
                end if

                mina = minmax(fragment_array(i)%Z(j),1)
                maxa = minmax(fragment_array(i)%Z(j),2)
                if (fragment_array(i)%A(j) .gt. maxa) then
                    ! Rest goes to neutrons
                    nem = nem + (fragment_array(i)%A(j)-maxa)
                    ! Yield of neutrons is the amount scaled with the
                    ! yield of the not included fragment
                    ! Y * A should be conserved
                    Ynem = Ynem + fragment_array(i)%Y(j)*(fragment_array(i)%A(j)-maxa)
                    fragment_array(i)%Y(j) = fragment_array(i)%Y(j) ! THis stays conserved
                    fragment_array(i)%A(j) = maxa
                    fragment_array(i)%net_idx(j) = findaz(maxa,fragment_array(i)%Z(j))
                end if
            end if
        end do

        ! Check how many neutrons are avaiable and give it to the missing
        ! fragments
        do j=1,nof
            if (fragment_array(i)%net_idx(j) .eq. -1) then

                ! Say something if verbose
                if (VERBOSE_LEVEL .ge. 2) then
                    write(*,*) "Fission fragment Z=",fragment_array(i)%Z(j), &
                               "A=",fragment_array(i)%A(j), &
                               "Y=",fragment_array(i)%Y(j), &
                               "not in network, trying to remap."
                end if

                ! Needed neutrons are:
                mina = minmax(fragment_array(i)%Z(j),1)
                bn = (mina-fragment_array(i)%A(j)) ! Amount of missing neutrons

                ! Borrow neutrons, take heaviest fragment that has enough
                ! abundance
                do m=nof,1,-1
                    if ((fragment_array(i)%net_idx(m) .ne. -1) .and. &
                        ((bn*fragment_array(i)%Y(j) / fragment_array(i)%A(m)) .lt. fragment_array(i)%Y(m))) then
                        exit
                    end if
                end do

                ! None found? Too bad, make an error
                if (m .eq. 0) then
                    call raise_exception("Could not remap fission fragments, "// &
                                         "consider to include more nuclides in the network. "//&
                                         NEW_LINE('A')//&
                                         "Nucleus Z="//trim(int_to_str(fragment_array(i)%Z(j)))// &
                                         ", A="//trim(int_to_str(fragment_array(i)%A(j)))// &
                                         " was not included. "//NEW_LINE('A')//&
                                         "Need "//int_to_str(bn)//" neutrons.",&
                                         "read_fiss_fragment_file",&
                                         190008)
                end if

                ! Otherwise borrow the neutrons from the heavier fragment
                fragment_array(i)%Y(m) = fragment_array(i)%Y(m) - (bn*fragment_array(i)%Y(j) / fragment_array(i)%A(m))
                fragment_array(i)%Y(j) = fragment_array(i)%Y(j)
                ! Give the fragment a new identity
                fragment_array(i)%A(j) = mina
                fragment_array(i)%net_idx(j) = findaz(mina,fragment_array(i)%Z(j))
            end if
        end do

        ! Proof the yields, this should be Af!
        helper = Ynem
        do j=1,nof
            helper = helper+fragment_array(i)%A(j)*fragment_array(i)%Y(j)
        end do
        if (abs(helper - float(fragment_array(i)%Ap)) .gt. 1e-2) then
            call raise_exception("Yield of fission fragments is not conserved. "// &
                                 "This happened for fragment distribution of "// &
                                 "Nucleus Z="//trim(int_to_str(fragment_array(i)%Zp))// &
                                 ", A="//trim(int_to_str(fragment_array(i)%Ap))// &
                                 ". Ensure that the sum is 2.","read_fiss_fragment_file",&
                                 190009)
        end if

        fragment_array(i)%neutrons = nem
        fragment_array(i)%Yn       = Ynem
        nufiss                     = nufiss + Ynem
    end do fragloop

    ! Close the file
    close(file_id)

    INFO_EXIT("read_fiss_fragment_file")

   end subroutine read_fiss_fragment_file


   !> Fill the rates with the correct fragments
   !!
   !! This subroutine fills the rates with the correct fragments
   !! from [Mumpower et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020PhRvC.101e4607M/abstract).
   !! This subroutine makes use of the fragtype struct.
   !! If no fragment is found, the fragment distribution is set to the
   !! one given by missing_frags.
   !! For beta-delayed fission the fissioning nucleus is the one of (Z+1, A)
   !! for neutron induced fission the nucleus with (Z, A+1), and for spontanous
   !! fission (not used within fissflag = 3) it is (Z,A).
   !!
   !! @author M. Reichert
   !! @date 14.02.2023
   subroutine file_fiss_frag(fissrate_inout,missing_frags)
    use global_class, only: isotope,ineu
    use benam_class, only: minmax,findaz
    use error_msg_class, only: raise_exception, int_to_str
    implicit none
       type(fissionrate_type),intent(inout)   :: fissrate_inout !< Temporary fission rate
       integer,optional, intent(in)           :: missing_frags!< Distribution to use if not in the file (1:panov, 2:kodama)
       ! Internal variables
       integer                                :: mass         !< Mass number of reacting nucleus
       integer                                :: pnr          !< Proton number of reacting nucleus
       real(r_kind)                           :: qval         !< Q-value of the reaction
       integer                                :: nfrac_distr !< Number of fragment distributions
       integer                                :: i           !< Loop variable
       integer                                :: astat       !< Allocation status variable
       integer                                :: ind_parent  !< Index of reacting nucleus
       integer                                :: afiss       !< A of fissioning nucleus
       integer                                :: zfiss       !< Z of fissioning nucleus
       type(fragtype)                         :: fissfrags_tmp!< Temporary fragment distribution
       logical                                :: found       !< Helper variable flag if the fragment
                                                             !  distribution was found
       integer                                :: additional_neutrons    !< Neutrons to be added for beta-delayed fission
       integer                                :: missing_frags_internal !< Internal variable for missing fragments

       INFO_ENTRY("file_fiss_frag")

       ! Set the default missing fragments
       if (.not. present(missing_frags)) then
          missing_frags_internal = 0
       else
          missing_frags_internal = missing_frags
       end if

      ! Get the mass and proton number
      mass = isotope(fissrate_inout%fissnuc_index)%mass
      pnr  = isotope(fissrate_inout%fissnuc_index)%p_nr

      ! Loop through the fission fragments
       if (fissrate_inout%reac_type .eq. rrt_bf) then
          nfrac_distr = size(fissfrags_beta_delayed)
       else if (fissrate_inout%reac_type .eq. rrt_nf) then
          nfrac_distr = size(fissfrags_n_induced)
       else if (fissrate_inout%reac_type .eq. rrt_sf) then
          nfrac_distr = size(fissfrags_spontaneous)
       end if

       ! Assume first that there are no additional neutrons
       additional_neutrons = fissrate_inout%released_n

       ! Check which nuclei is fissioning
       zfiss = pnr
       afiss = mass - additional_neutrons

       ! Z-1 for beta delayed fission
       if (fissrate_inout%reac_type .eq. rrt_bf) then
          zfiss = zfiss+1
       ! A+1 for neutron induced fission
       elseif (fissrate_inout%reac_type .eq. rrt_nf) then
          afiss = afiss+1
       end if
       ! Don't change anything for spontaneous fission

       ! Find relevant indices
       ind_parent = findaz(mass,pnr)
       fissrate_inout%fissnuc_index = ind_parent
       found = .false.
        do i=1,nfrac_distr

            if (fissrate_inout%reac_type .eq. rrt_bf) then
               fissfrags_tmp = fissfrags_beta_delayed(i)
            else if (fissrate_inout%reac_type .eq. rrt_nf) then
               fissfrags_tmp = fissfrags_n_induced(i)
            else if (fissrate_inout%reac_type .eq. rrt_sf) then
               fissfrags_tmp = fissfrags_spontaneous(i)
            end if

            if ((fissfrags_tmp%Zp .eq. zfiss) .and. (fissfrags_tmp%Ap .eq. afiss)) then
                found = .true.
                exit
            end if
        end do

        ! Not found, check what to do
        if (.not. found) then

            ! Raise an error when there is no alternative distribution specified
            if (missing_frags_internal .eq. 0) then
                call raise_exception("No fragment distribution found for Z="//int_to_str(zfiss)//", "//&
                                     "A="//int_to_str(afiss)//" in fragment file. Either add it to the fission fragment file or set "//&
                                     "'fission_frag_missing' to 1 (Panov et al. 2001) or 2 (Kodama & Takahashi 1975)",&
                                     "file_fiss_frag", 190016)
            end if

            ! Output this information
            if (VERBOSE_LEVEL .ge.2) then
                write(*,*) "No fragment distribution found for Z=",zfiss,", A=",afiss
                if (missing_frags_internal .eq. 1) then
                    write(*,*) "Using Panov 2001 instead."
                else if (missing_frags_internal .eq. 2) then
                    write(*,*) "Using Kodama & Takahashi 1975 instead."
                end if
            end if

            if (missing_frags_internal .eq. 1) then
                ! Get panov distribution
                call fiss_dist(fissrate_inout)
            else if (missing_frags_internal .eq. 2) then
                ! Get kodama distribution
                call kodtakdist(fissrate_inout)
            else
                call raise_exception("Unknown missing fragment distribution. "//&
                                     "Got "//int_to_str(missing_frags_internal)//".",&
                                     "file_fiss_frag",190015)
            end if
        else
            ! Found, do mumpower
            if (fissrate_inout%reac_type .eq. rrt_bf) then
                fissfrags_tmp = fissfrags_beta_delayed(i)
             else if (fissrate_inout%reac_type .eq. rrt_nf) then
                fissfrags_tmp = fissfrags_n_induced(i)
             else if (fissrate_inout%reac_type .eq. rrt_sf) then
                fissfrags_tmp = fissfrags_spontaneous(i)
             end if

            ! Allocate the parts, the first two are the fissioning nucleus and neutrons
            allocate(fissrate_inout%fissparts(fissfrags_tmp%nr_frags+2),&
                     fissrate_inout%ch_amount(fissfrags_tmp%nr_frags+2),&
                     stat=astat)
            ! Complain if not possible
            if (astat .ne. 0) then
                call raise_exception("Could not allocate memory for fission fragments.",&
                                     "file_fiss_frag",190001)
            end if

            ! Set up the number of participants,
            ! Number of fragments + parent nucleus + neutrons
            fissrate_inout%dimens = fissfrags_tmp%nr_frags+2

            ! Use the correct order of the reacting nuclei
            if (.not. (fissrate_inout%reac_type .eq. rrt_nf)) then
                fissrate_inout%fissparts(1) = ind_parent
                fissrate_inout%fissparts(2) = ineu
                fissrate_inout%ch_amount(1) = -1d0
                fissrate_inout%ch_amount(2) = fissfrags_tmp%Yn+float(additional_neutrons)
            else
                fissrate_inout%fissparts(1) = ineu
                fissrate_inout%fissparts(2) = ind_parent
                fissrate_inout%ch_amount(1) = fissfrags_tmp%Yn-1d0+float(additional_neutrons)
                fissrate_inout%ch_amount(2) = -1d0
            end if


            do i=1,fissfrags_tmp%nr_frags
                fissrate_inout%fissparts(i+2) = fissfrags_tmp%net_idx(i)
                fissrate_inout%fissparts(i+2) = fissfrags_tmp%net_idx(i)
                fissrate_inout%ch_amount(i+2) = fissfrags_tmp%Y(i)
            end do
        end if

        ! Calculate the Q-value
        qval = 0.d0
        qvalue: do i=1,fissrate_inout%dimens
           qval = qval - isotope(fissrate_inout%fissparts(i))%mass_exc*fissrate_inout%ch_amount(i) ! weighted average Q-value
        end do qvalue
        fissrate_inout%averageQ = qval

       INFO_EXIT("file_fiss_frag")

    end subroutine file_fiss_frag



   !> Calculates the (neutron-induced) fission fragment mass distribution according to
   !! the ABLA07 model: Kelic, Ricciardi, & Schmidt (arXiv:0906.4193)
   !!
   !! @see [Kelic et al. 2009](https://ui.adsabs.harvard.edu/abs/2009arXiv0906.4193K/abstract)
   subroutine abla_nfiss(pos,mass,pnr,neutfission,qval)
      use global_class, only: isotope,ineu
      use benam_class, only: minmax,findaz
      use error_msg_class, only: raise_exception

         integer,intent(in)                     :: mass,pnr
         integer, intent(in)                    :: pos
         real(r_kind),intent(out)               :: qval
         integer                                :: fmass,fpnr
         integer                                :: read_stat
         integer                                :: alloc_stat
         integer                                :: af,zf
         integer,dimension(:),allocatable       :: a1,z1
         integer,dimension(:),allocatable       :: a2,z2
         integer,intent(in)                     :: neutfission
         integer                                :: i,nof, nemiss
         integer                                :: dimens
         real(r_kind)                           :: nrate,npre,nafter
         real(r_kind),dimension(:),allocatable  :: paz
         real(r_kind)                           :: psum
         real(r_kind)                           :: n_mult
         character(2)                           :: dummy

         INFO_ENTRY("abla_nfiss")

         nemiss = 0
         do
           read(neutfission,"(I3,1X,I3,5X,I3)",iostat=read_stat) zf, af, nof
           if (read_stat /= 0) then            !< this condition can happen and the code should be allowed to continue (stop -> exit)
   !           print*, 'fission reaction not found in abla_nfiss: ', fmass, fpnr
              exit     ! parent nucleus not found in ABLA table -> continue with old Panov model (below, corresponds to fissflag=1 for this nucleus)
           endif

           read(neutfission,*,iostat=read_stat) nrate, npre, nafter
           if (read_stat /= 0) then
              call raise_exception("Could not read fission file.","abla_nfiss",190004)
           endif

           if ((af.eq.mass).and.(zf.eq.pnr).and.(nof.ne.0)) then

             n_mult = npre + nafter
             nemiss = int(dnint(n_mult))
             nufiss = nufiss + nemiss

             allocate(a1(nof),stat=alloc_stat)
             allocate(a2(nof),stat=alloc_stat)
             allocate(z1(nof),stat=alloc_stat)
             allocate(z2(nof),stat=alloc_stat)
             allocate(paz(nof),stat=alloc_stat)

             fissrate(pos)%channels       = nof
             fissrate(pos)%fissnuc_index  = findaz(af,zf)
             fissrate(pos)%mode           = 1

             allocate(fissrate(pos)%channelprob(nof))
             allocate(fissrate(pos)%q_value(nof))
             dimens = 2 * fissrate(pos)%channels + 2
             fissrate(pos)%dimens = dimens
             allocate(fissrate(pos)%fissparts(dimens))
             allocate(fissrate(pos)%ch_amount(dimens))
             fissrate(pos)%fissparts(1) = ineu
             fissrate(pos)%ch_amount(1) = float(nemiss) - 1.d0
             fissrate(pos)%fissparts(2) = fissrate(pos)%fissnuc_index
             fissrate(pos)%ch_amount(2) = -1.d0


             psum = 0.d0
             inner: do i=1,nof
               read(neutfission,*,iostat=read_stat) &
              &                   a1(i),z1(i),paz(i)
               if (read_stat /= 0) then
                  call raise_exception("Could not read fission file.","abla_nfiss",190004)
               endif
   ! check if fragment 1 is in the network
               if (a1(i).gt.minmax(z1(i),2)) then
                  nemiss = nemiss + a1(i) - minmax(z1(i),2)
                  a1(i) = minmax(z1(i),2)
               else if (a1(i).lt.minmax(z1(i),1)) then
                  print *, 'a1 error:'
               end if
               a2(i) = mass + 1 - a1(i) - nemiss
               z2(i) = pnr - z1(i)
               psum = psum + paz(i)

   ! check if fragment 2 is in the network
               if (a2(i).gt.minmax(z2(i),2)) then
                  nemiss = nemiss + a2(i) - minmax(z2(i),2)
                  a2(i) = minmax(z2(i),2)
               else if (a2(i).lt.minmax(z2(i),1)) then
                  print *, 'a2 error:', a2(i), z2(i), a1(i), z1(i), mass, pnr
               end if
             end do inner

             again: do i=1,nof
               paz(i) = paz(i)/psum        ! normalize probabilities
               fissrate(pos)%channelprob(i) = paz(i)

   !---------rate for neutron-induced fission
               fissrate(pos)%fissparts(i+2) = findaz(a1(i),z1(i))
               fissrate(pos)%ch_amount(i+2) = fissrate(pos)%channelprob(i)   ! rate at which fragment is produced per destroyed parent nucleus
               fissrate(pos)%fissparts(i+2+fissrate(pos)%channels) = findaz(a2(i),z2(i))
               fissrate(pos)%ch_amount(i+2+fissrate(pos)%channels) = fissrate(pos)%channelprob(i)
               fissrate(pos)%q_value(i) = (1-nemiss) * isotope(ineu)%mass_exc + &  ! one neutron is destroyed, a number equal to nemiss are produced
                       isotope(fissrate(pos)%fissnuc_index)%mass_exc - &
                       isotope(fissrate(pos)%fissparts(i+2))%mass_exc - &
                       isotope(fissrate(pos)%fissparts(i+2+fissrate(pos)%channels))%mass_exc

             end do again

             qval = 0.d0
             qvalue: do i=1,fissrate(pos)%channels
                qval = qval + fissrate(pos)%q_value(i)*fissrate(pos)%channelprob(i) ! weighted average Q-value
             end do qvalue
             fissrate(pos)%averageQ = qval


             deallocate(a1)
             deallocate(a2)
             deallocate(z1)
             deallocate(z2)
             deallocate(paz)

             rewind(neutfission,iostat=read_stat)
             if (read_stat /= 0) then
                call raise_exception("Could not rewind fission file.",&
                                     "abla_nfiss",190005)
             endif
             return

           else           ! read nof lines to advance to next parent nucleus in the table

             do i=1,nof
               read(neutfission,*,iostat=read_stat) dummy
               if (read_stat /= 0) then
                  call raise_exception("Could not read fission file.",&
                                       "abla_nfiss",190004)
               endif
             end do

           end if

         end do

   !--------old Panov model for nuclei not found in ABLA fission table---------------
         allocate(a1(1),stat=alloc_stat)
         allocate(a2(1),stat=alloc_stat)
         allocate(z1(1),stat=alloc_stat)
         allocate(z2(1),stat=alloc_stat)

         fmass = mass + 1
         fpnr  = pnr
         select case (fmass)
             case(:240) ! -> symmetric fission
               a2(1) = nint(fmass/2.d0)
               z2(1) = nint(fpnr/2.d0)
               a1(1) = fmass - a2(1)
               z1(1) = fpnr - z2(1)
             case(255:265) ! -> symmetric fission
               a2(1) = nint(fmass/2.d0)
               z2(1) = nint(fpnr/2.d0)
               a1(1) = fmass - a2(1)
               z1(1) = fpnr - z2(1)
             case default ! -> asymmetric fission
               a2(1) = 130
               z2(1) = nint(52.01 - (fpnr - 80)/1.d1)
   !            a2 = max(a2,af-a2)     !modification to prevent very n-rich
   !            z2 = max(z2,zf-z2)     !fragments
               a1(1) = fmass - a2(1)
               z1(1) = fpnr - z2(1)
         end select

         nemiss = 0

         if (a1(1).gt.minmax(z1(1),2)) then
             nemiss = a1(1) - minmax(z1(1),2)
             a1(1) = minmax(z1(1),2)
         else if (a1(1).lt.minmax(z1(1),1)) then
             print *, 'a1 error:'
         end if

         if (a2(1).gt.minmax(z2(1),2)) then
             nemiss = nemiss + a2(1) - minmax(z2(1),2)
             a2(1) = minmax(z2(1),2)
         else if (a2(1).lt.minmax(z2(1),1)) then
             print *, 'a2 error:'
         end if
         nufiss = nufiss + nemiss

         fissrate(pos)%dimens = 4
         allocate(fissrate(pos)%fissparts(fissrate(pos)%dimens))
         allocate(fissrate(pos)%ch_amount(fissrate(pos)%dimens))
         allocate(fissrate(pos)%channelprob(1))
         allocate(fissrate(pos)%q_value(1))
         fissrate(pos)%channels       = 1
         fissrate(pos)%fissnuc_index  = findaz(mass,pnr)
         fissrate(pos)%mode           = 1
         fissrate(pos)%channelprob(1) = 1.d0

         fissrate(pos)%fissparts(1) = ineu
         fissrate(pos)%ch_amount(1) = float(nemiss) - 1.d0
         fissrate(pos)%fissparts(2) = fissrate(pos)%fissnuc_index
         fissrate(pos)%ch_amount(2) = -1.d0

         fissrate(pos)%fissparts(3) = findaz(a1(1),z1(1))
         fissrate(pos)%ch_amount(3) = 1.d0
         fissrate(pos)%fissparts(4) = findaz(a2(1),z2(1))
         fissrate(pos)%ch_amount(4) = 1.d0

         fissrate(pos)%averageQ  = (1-nemiss) * isotope(ineu)%mass_exc + &
                                 isotope(fissrate(pos)%fissparts(2))%mass_exc - &
                                 isotope(fissrate(pos)%fissparts(3))%mass_exc - &
                                 isotope(fissrate(pos)%fissparts(4))%mass_exc

         deallocate(a1)
         deallocate(a2)
         deallocate(z1)
         deallocate(z2)


         rewind(neutfission,iostat=read_stat)
         if (read_stat /= 0) then
            call raise_exception("Could not rewind fission file.",&
                                 "abla_nfiss",190005)
         endif

         INFO_EXIT("abla_nfiss")

   end subroutine abla_nfiss


   !> Calculates the (beta-delayed and spontaneous) fission fragment mass distribution according to
   !! the ABLA07 model: Kelic, Ricciardi, & Schmidt (arXiv:0906.4193)
   !!
   !! @see [Kelic et al. 2009](https://ui.adsabs.harvard.edu/abs/2009arXiv0906.4193K/abstract)
   subroutine abla_betafiss(pos,mass,pnr,betafission,qval)
      use global_class,    only: isotope,ineu
      use benam_class,     only: minmax,findaz
      use error_msg_class, only: raise_exception
         integer,intent(in)                     :: mass,pnr
         integer, intent(inout)                 :: pos
         real(r_kind),intent(out)               :: qval
         integer                                :: nemiss,fmass,fpnr
         integer                                :: read_stat
         integer                                :: alloc_stat
         integer                                :: af,zf
         integer,dimension(:),allocatable       :: a1,z1
         integer,dimension(:),allocatable       :: a2,z2
         integer,intent(in)                     :: betafission
         integer                                :: i,nof
         integer                                :: dimens
         real(r_kind)                           :: nrate,npre,nafter
         real(r_kind),dimension(:),allocatable  :: paz
         real(r_kind)                           :: psum
         real(r_kind)                           :: n_mult
         character(2)                           :: dummy

         INFO_ENTRY("abla_betafiss")

         nemiss = 0
         do
            read(betafission,*, iostat = read_stat) zf, af, nof

           if (read_stat /= 0) exit
           read(betafission,*,iostat=read_stat) nrate, npre, nafter
           if (read_stat /= 0) then
              call raise_exception("Could not read fission file.",&
                                   "abla_betafiss",190006)
           endif
           if ((af.eq.mass).and.(zf.eq.pnr).and.(nof.ne.0)) then

             n_mult = npre + nafter
             nemiss = int(dnint(n_mult))
             nufiss = nufiss + nemiss

             allocate(a1(nof),stat=alloc_stat)
             allocate(a2(nof),stat=alloc_stat)
             allocate(z1(nof),stat=alloc_stat)
             allocate(z2(nof),stat=alloc_stat)
             allocate(paz(nof),stat=alloc_stat)

             fissrate(pos)%channels       = nof
             fissrate(pos)%fissnuc_index  = findaz(mass,pnr)
             if (fissrate(pos)%reac_type.eq. rrt_sf) then
                fissrate(pos)%mode = 2                ! spontaneous fission
             else if (fissrate(pos)%reac_type.eq. rrt_bf) then
                fissrate(pos)%mode = 3                ! beta-delayed fission
             end if

             allocate(fissrate(pos)%channelprob(nof))
             allocate(fissrate(pos)%q_value(nof))
             dimens = 2 * fissrate(pos)%channels + 1
             if (nemiss.gt.0) dimens = dimens + 1                            ! neutrons which would otherwise not be part of the reaction
             fissrate(pos)%dimens = dimens
             allocate(fissrate(pos)%fissparts(dimens))
             allocate(fissrate(pos)%ch_amount(dimens))
             fissrate(pos)%fissparts(1) = fissrate(pos)%fissnuc_index
             fissrate(pos)%ch_amount(1) = -1.d0

             if (nemiss.gt.0) then
                fissrate(pos)%fissparts(dimens) = ineu
                fissrate(pos)%ch_amount(dimens) = float(nemiss)
             end if

             psum = 0.d0
             inner: do i=1,nof
               read(betafission,"(I3,1X,I3,2X,F6.4)",iostat=read_stat) a1(i), z1(i), paz(i)
               if (read_stat /= 0) then
                  call raise_exception("Could not read fission file.",&
                                       "abla_betafiss",190006)
               endif

   ! check if fragment 1 is in the network
               if (a1(i).gt.minmax(z1(i),2)) then
                  nemiss = nemiss + a1(i) - minmax(z1(i),2)
                  a1(i) = minmax(z1(i),2)
               else if (a1(i).lt.minmax(z1(i),1)) then
                  print *, 'a1 error:'
               end if

               a2(i) = mass - a1(i) - nemiss
               z2(i) = pnr - z1(i)
               if (fissrate(pos)%mode.eq.3) z2(i) = z2(i) + 1      ! beta-delayed fission
               psum = psum + paz(i)

   ! check if fragment 2 is in the network
               if (a2(i).gt.minmax(z2(i),2)) then
                  nemiss = nemiss + a2(i) - minmax(z2(i),2)
                  a2(i) = minmax(z2(i),2)
               else if (a2(i).lt.minmax(z2(i),1)) then
                  print *, 'a2 error:',a2(i),minmax(z2(i),1)
               end if
             end do inner

             again: do i=1,nof
               paz(i) = paz(i)/psum            ! normalize probabilities

               fissrate(pos)%channelprob(i) = paz(i)

   !---------rate for spontaneous and beta-delayed fission
               fissrate(pos)%fissparts(i+1) = findaz(a1(i),z1(i))                         ! first fragment
               fissrate(pos)%ch_amount(i+1) = fissrate(pos)%channelprob(i)                ! rate at which fragment is produced per destroyed parent nucleus
               fissrate(pos)%fissparts(i+1+fissrate(pos)%channels) = findaz(a2(i),z2(i))  ! second fragment
               fissrate(pos)%ch_amount(i+1+fissrate(pos)%channels) = fissrate(pos)%channelprob(i)
               fissrate(pos)%q_value(i) = isotope(fissrate(pos)%fissnuc_index)%mass_exc - &
                       nemiss * isotope(ineu)%mass_exc - &
                       isotope(fissrate(pos)%fissparts(i+1))%mass_exc - &
                       isotope(fissrate(pos)%fissparts(i+1+fissrate(pos)%channels))%mass_exc
             end do again


             qval = 0.d0
             qvalue: do i=1,fissrate(pos)%channels
                qval = qval + fissrate(pos)%q_value(i)*fissrate(pos)%channelprob(i)  ! weighted average Q-value
             end do qvalue
             fissrate(pos)%averageQ = qval

             deallocate(a1)
             deallocate(a2)
             deallocate(z1)
             deallocate(z2)
             deallocate(paz)

             rewind(betafission,iostat=read_stat)
             if (read_stat /= 0) then
                call raise_exception("Could not rewind fission file.",&
                                     "abla_betafiss",190007)
             endif
             return

           else

             do i=1,nof
               read(betafission,*,iostat=read_stat) dummy
               if (read_stat /= 0) then
                  call raise_exception("Could not read fission file.",&
                                       "abla_betafiss",190007)
               endif
             end do

           end if

         end do

  !--------old Panov model for nuclei not found in ABLA fission file---------------
         allocate(a1(1),stat=alloc_stat)
         allocate(a2(1),stat=alloc_stat)
         allocate(z1(1),stat=alloc_stat)
         allocate(z2(1),stat=alloc_stat)

         fmass = mass
         if (fissrate(pos)%reac_type.eq. rrt_sf) then
            fpnr = pnr
            fissrate(pos)%mode = 2             ! spontaneous fission
         else if (fissrate(pos)%reac_type.eq. rrt_bf) then
            fpnr  = pnr + 1
            fissrate(pos)%mode = 3              ! beta-delayed fission
         end if

         select case (fmass)
             case(:240) ! -> symmetric fission
               a2(1) = nint(fmass/2.d0)
               z2(1) = nint(fpnr/2.d0)
               a1(1) = fmass - a2(1)
               z1(1) = fpnr - z2(1)
             case(255:265) ! -> symmetric fission
               a2(1) = nint(fmass/2.d0)
               z2(1) = nint(fpnr/2.d0)
               a1(1) = fmass - a2(1)
               z1(1) = fpnr - z2(1)
             case default ! -> asymmetric fission
               a2(1) = 130
               z2(1) = nint(52.01 - (fpnr - 80)/1.d1)
   !            a2 = max(a2,af-a2)     !modification to prevent very n-rich
   !            z2 = max(z2,zf-z2)     !fragments
               a1(1) = fmass - a2(1)
               z1(1) = fpnr - z2(1)
         end select

         nemiss = 0

         if (a1(1).gt.minmax(z1(1),2)) then
             nemiss = a1(1) - minmax(z1(1),2)
             a1(1) = minmax(z1(1),2)
         else if (a1(1).lt.minmax(z1(1),1)) then
             print *, 'a1 error:'
         end if

         if (a2(1).gt.minmax(z2(1),2)) then
             nemiss = nemiss + a2(1) - minmax(z2(1),2)
             a2(1) = minmax(z2(1),2)
         else if (a2(1).lt.minmax(z2(1),1)) then
             print *, 'a2 error:',a2(1),minmax(z2(1),1)
         end if
         nufiss = nufiss + nemiss

         if (nemiss .eq. 0) then
            fissrate(pos)%dimens = 3
         else
            fissrate(pos)%dimens = 4
         end if

         allocate(fissrate(pos)%fissparts(fissrate(pos)%dimens))
         allocate(fissrate(pos)%ch_amount(fissrate(pos)%dimens))
         allocate(fissrate(pos)%channelprob(1))
         allocate(fissrate(pos)%q_value(1))
         fissrate(pos)%channels       = 1
         fissrate(pos)%fissnuc_index  = findaz(mass,pnr)
         fissrate(pos)%channelprob(1) = 1.d0

         fissrate(pos)%fissparts(1) = fissrate(pos)%fissnuc_index
         fissrate(pos)%ch_amount(1) = -1.d0

         fissrate(pos)%fissparts(2) = findaz(a1(1),z1(1))
         fissrate(pos)%ch_amount(2) = 1.d0
         fissrate(pos)%fissparts(3) = findaz(a2(1),z2(1))
         fissrate(pos)%ch_amount(3) = 1.d0
         if (nemiss.gt.0) then
            fissrate(pos)%fissparts(4) = ineu
            fissrate(pos)%ch_amount(4) = float(nemiss)
         end if

         fissrate(pos)%averageQ = isotope(fissrate(pos)%fissparts(1))%mass_exc - &
                                 isotope(fissrate(pos)%fissparts(2))%mass_exc - &
                                 isotope(fissrate(pos)%fissparts(3))%mass_exc - &
                                 nemiss * isotope(ineu)%mass_exc

         deallocate(a1)
         deallocate(a2)
         deallocate(z1)
         deallocate(z2)


         rewind(betafission,iostat=read_stat)
         if (read_stat /= 0) then
            call raise_exception("Could not rewind fission file.",&
                                 "abla_betafiss",190007)
         endif

         INFO_EXIT("abla_betafiss")

   end subroutine abla_betafiss



end module fission_rate_module
