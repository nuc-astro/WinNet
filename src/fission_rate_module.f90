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
   use global_class, only: reactionrate_type
   implicit none


   !> fission rate type, designed to save fragment distribution at initialization
   type,public                     :: fissionrate_type
      integer                      :: fissnuc_index        !< index of parent nucleus
      integer                      :: channels             !< number of fragment pairs for each fission rate
      integer                      :: mode                 !< fission mode; 1:n-induced fission 2: spontaneous and beta-delayed fission
      integer                      :: dimens               !< dimension of fissparts array
      real(r_kind)                 :: cached               !< computed rate
      real(r_kind)                 :: averageQ             !< average Q value, weighted over all fragment channels
      integer                      :: reac_type            !< "rrt_sf": spontaneous fission, "rrt_bf": beta-delayed fission, "rrt_nf": neutron-induced fission
      integer,dimension(:),allocatable      :: fissparts   !< corresponds to rrate()%parts() array
      integer,dimension(:,:),allocatable    :: cscf_ind    !< cscf_ind (i,j) gives the position of the entry \f$\frac{d\dot{Y}_{j}}{dY_{i}}\f$ in the cscf data array
      real(r_kind),dimension(9)             :: param       !< rate parameters a0-a9, as given in reaclib file
      real(r_kind),dimension(:),allocatable :: q_value     !< Q value of each fission channel
      real(r_kind),dimension(:),allocatable :: channelprob !< probability for each fragment pair
      real(r_kind),dimension(:),allocatable :: ch_amount   !< corresponds to rrate(i)%ch_amount(j); for fragment i ch_amount(i) = +1 * channelprob(i)
   end type fissionrate_type

   type(fissionrate_type),dimension(:),allocatable,public   :: fissrate !< Array storing fission reactions
   type(reactionrate_type),dimension(:),allocatable,private :: rrate_fiss !< Array storing fission reactions in rrate array


   character(len=*), private, parameter                     :: fiss_binary_name='fiss_rates.windat' !< Name of the binary file containing the fission rates


    type,private                              :: fragtype
        integer                               :: Zp  !< Z of fissioning nucleus
        integer                               :: Ap  !< A of fissioning nucleus
        integer                               :: neutrons !< number of neutrons emitted
        integer                               :: nr_frags
        real(r_kind)                          :: Yn  !< Yield for neutrons
        integer,dimension(:),allocatable      :: net_idx !< index of fragment in network
        integer,dimension(:),allocatable      :: Z   !< Z of fragment
        integer,dimension(:),allocatable      :: A   !< A of fragment
        real(r_kind),dimension(:),allocatable :: Y   !< Yields of fragment
    end type fragtype
    type(fragtype),dimension(:),allocatable,private   :: fissfrags !< Array storing fragment distributions




   integer, public, parameter    :: fiss_neglect=5 !< Amount of fission fragments not to be neglected
                                                   !  in the jacobian in case of vanishing parent
   integer, public               :: nfiss        !< Amount of fission rates
   integer, private              :: nufiss       !< total number of fission neutrons

   integer, private  :: n_nf,n_bf,n_sf !< Amount of individual reactions
   !
   ! Public and private fields and methods of the module
   !
   public:: &
      init_fission_rates, merge_fission_rates, output_binary_fission_reaction_data
   private:: &
      count_fission_rates, read_fission_rates, fiss_dist, kodtakdist,&
      abla_nfiss, abla_betafiss, read_mumpower_fissfile, mumpower_fiss,&
      reorder_fragments
contains


   !> Initialize the fission reactions
   !!
   !! This subroutine counts and reads fission reactions.
   !! After calling it, the rate array \ref rrate_fiss and
   !! \ref fissrate will be filled as well as the integer
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
         allocate(rrate_fiss(nfiss),stat=alloc_stat)
         if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_fiss" failed.',&
                                                    "init_fission_rates",190001)
         allocate(fissrate(nfiss),stat=alloc_stat)
         if ( alloc_stat /= 0) call raise_exception('Allocation of "fissrate" failed.',&
                                                    "init_fission_rates",190001)

         !-- Read the fission rates into both arrays
         call read_fission_rates()

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

    close(file_id)


   end subroutine output_binary_fission_reaction_data





   !> Merge fission rates with larger array.
   !!
   !! This subroutine will merge the fission rates
   !! that are represented by a usual \ref global_class::reactionrate_type
   !! into the large rate array. As the reactions should not be
   !! already included in the large array, it will simply append the reactions
   !! behind the rrate_array array.
   !!
   !! @note The other fission rate array \ref fissrate will exist independently
   !!       and are used to calculate the correct equations for the
   !!       fission fragments in \ref jacobian_class::jacobi_init,
   !!       \ref jacobian_class::abchange.
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine merge_fission_rates(rrate_array,rrate_length)
      use error_msg_class, only: raise_exception
      use parameter_class, only: use_prepared_network
      implicit none
      type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate_array  !< Large rate array, containing all reactions
      integer,intent(inout)                                          :: rrate_length !< length of rrate_array
      type(reactionrate_type),dimension(:),allocatable               :: rrate_tmp    !< Temporary rate array
      integer                                                        :: alloc_stat   !< Allocation state
      integer                                                        :: new_length   !< New length of rrate_array

      !-- Only do something if there are fission rates
      if (allocated(rrate_fiss) .and. (.not. use_prepared_network)) then
         ! New length of the array
         new_length = rrate_length+nfiss
         if (nfiss .ne. 0) then
           if (.not. allocated(rrate_array)) then
              !-- Allocate the reaclib rate array
              allocate(rrate_array(nfiss),stat=alloc_stat)
              if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                         "merge_fission_rates",190001)
              rrate_array(1:nfiss) = rrate_fiss(1:nfiss)
           else
              !-- Allocate a temporary array to store the content of rrate_array
              allocate(rrate_tmp(rrate_length),stat=alloc_stat)
              if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_tmp" failed.',&
                                                         "merge_fission_rates",190001)
              rrate_tmp(1:rrate_length) = rrate_array(1:rrate_length)

              !-- Deallocate the input array
              deallocate(rrate_array)
              allocate(rrate_array(new_length),stat=alloc_stat)
              if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                         "merge_fission_rates",190001)
              rrate_array(1:rrate_length)             = rrate_tmp(1:rrate_length)
              rrate_array(rrate_length+1:new_length)  = rrate_fiss(1:nfiss)

              deallocate(rrate_tmp)
           end if
           !-- Deallocate the fission rate array
           deallocate(rrate_fiss)
           !-- Output the new length
           rrate_length = new_length
         end if
      end if
   end subroutine merge_fission_rates


   !> Count the amount of fission rates
   !!
   !! This subroutine counts the amount of fission rates and
   !! stores the result in \ref nfiss
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine count_fission_rates()
      use file_handling_class
      use parameter_class, only: fission_rates
      use benam_class,     only: benam
      use format_class
      implicit none
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

      ifiss = 0

      ! Open fission file
      fissionlib = open_infile(fission_rates)

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
 !----- parts_index(i) = 0 if nuclide i is not in network
            if (parts_index(j) .eq. 0) cycle fission
         end do indices
         ifiss = ifiss + 1
      end do fission

      ! Close the fission file again
      close(fissionlib)


      ! Store the amount of fission rates
      nfiss = ifiss

   end subroutine


   !> Reads fission rates
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
   !! @returns \ref fissrate and \ref rrate_fiss filled with fission rates
   !! @author M. Eichler
   !!
   !! \b Edited:
   !!   - 25.01.21 (MR): Moved from init_network to this new independent subroutine
   !! -
   subroutine read_fission_rates()
      use parameter_class, only: fission_rates, fissflag, &
                                 nfission_file, bfission_file
      use benam_class,     only: benam
      use global_class,    only: isotope
      use file_handling_class
      use format_class
      use error_msg_class
      implicit none
      integer                    :: neutfission=0 !< File IDs for abla fiss. distr. file
      integer                    :: betafission=0 !< File IDs for abla fiss. distr. file
      integer                    :: read_stat     !< Read status
      integer                    :: fissionlib    !< File id of fisionrates
      character(5), dimension(6) :: parts         !< Participating nuclei names
      integer, dimension(6)      :: parts_index   !< Participating nuclei indices
      real(r_kind), dimension(9) :: params        !< Reaclib fitting parameter
      character(4)               :: src           !< Reaclib source label
      character(1)               :: res           !< Reaclib weak flag label
      character(1)               :: rev           !< Reaclib reverse label
      real(r_kind)               :: q             !< Reaclib Q-Value
      integer                    :: grp           !< Reaclib chapter
      integer                    :: group_index   !< Storage for the current chapter
      integer                    :: fissindex     !< Index for current fission rate
      integer                    :: j             !< Loop variable
      integer                    :: astat         !< Allocation status



      !-- Open the fission library file
      fissionlib = open_infile(fission_rates)

      !-- Open the abla files, containing probabilities of the fiss. fragments
      if (fissflag .eq. 3) then
         call read_mumpower_fissfile(nfission_file)
      else if (fissflag .eq. 4) then
         neutfission= open_infile(nfission_file)
         betafission= open_infile(bfission_file)
      end if

      ! make sure counting index for reaction rates continues for fission rates
      fissindex = 1
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
         if (src .eq. 'fiss') then
            rrate_fiss(fissindex)%reac_type = rrt_nf
            n_nf=n_nf+1 !< count neutron induced fission
         else if ((src .eq. 'mp01') .or. &
                 (src .eq. 'ms99')) then
            rrate_fiss(fissindex)%reac_type = rrt_bf
            n_bf=n_bf+1 !< count beta delayed fission
         else if (src .eq. 'sfis') then
            rrate_fiss(fissindex)%reac_type = rrt_sf
            n_sf=n_sf+1 !< count spontaneous fission
         end if

         if (rrate_fiss(fissindex)%reac_type.eq.rrt_nf) then      ! neutron-induced fission
            rrate_fiss(fissindex)%is_const = .false.
            rrate_fiss(fissindex)%is_weak  = .false.
            fissrate(fissindex)%reac_type = rrate_fiss(fissindex)%reac_type
            fissrate(fissindex)%param = params
            select case (fissflag)
            case(1)
               call fiss_dist(fissindex,isotope(parts_index(2))%mass,    &
                    isotope(parts_index(2))%p_nr,src,q)
            case(2)
               call kodtakdist(fissindex,isotope(parts_index(2))%mass,  &
                    isotope(parts_index(2))%p_nr,src,q)
            case(3)
               call mumpower_fiss(fissindex,isotope(parts_index(2))%mass,  &
                    isotope(parts_index(2))%p_nr,src,q)
            case(4)
               call abla_nfiss(fissindex,isotope(parts_index(2))%mass,  &
                    isotope(parts_index(2))%p_nr,neutfission,q)
             case default
                call raise_exception("Fission flag ("//trim(adjustl(int_to_str(fissflag)))&
                                     //") not implemented yet. "//&
                                     "Set it to a supported value!",&
                                     "read_fission_rates",&
                                     190003)
            end select
         else if (rrate_fiss(fissindex)%reac_type.eq.rrt_sf) then  ! spontaneous fission
            rrate_fiss(fissindex)%is_const = .true.
            rrate_fiss(fissindex)%is_weak  = .true.
            fissrate(fissindex)%reac_type = rrate_fiss(fissindex)%reac_type
            fissrate(fissindex)%param = params
            select case (fissflag)
            case(1)
               call fiss_dist(fissindex,isotope(parts_index(1))%mass,    &
                    isotope(parts_index(1))%p_nr,src,q)
            case(2)
               call kodtakdist(fissindex,isotope(parts_index(1))%mass,  &
                    isotope(parts_index(1))%p_nr,src,q)
            case(3)
               call kodtakdist(fissindex,isotope(parts_index(1))%mass,  &
                    isotope(parts_index(1))%p_nr,src,q)
            case(4)
               call abla_betafiss(fissindex,isotope(parts_index(1))%mass,  &
                    isotope(parts_index(1))%p_nr,src,betafission,q)
            end select
         else if (rrate_fiss(fissindex)%reac_type.eq.rrt_bf) then  ! beta-delayed fission
            rrate_fiss(fissindex)%is_const = .true.
            rrate_fiss(fissindex)%is_weak  = .true.
            fissrate(fissindex)%reac_type = rrate_fiss(fissindex)%reac_type
            fissrate(fissindex)%param = params
            select case (fissflag)
            case(1)
               call fiss_dist(fissindex,isotope(parts_index(1))%mass,    &
                    isotope(parts_index(1))%p_nr,src,q)
            case(2)
               call kodtakdist(fissindex,isotope(parts_index(1))%mass,  &
                    isotope(parts_index(1))%p_nr,src,q)
            case(3)
               call mumpower_fiss(fissindex,isotope(parts_index(1))%mass,  &
                    isotope(parts_index(1))%p_nr,src,q)
            case(4)
               call abla_betafiss(fissindex,isotope(parts_index(1))%mass,  &
                    isotope(parts_index(1))%p_nr,src,betafission,q)
            end select
         end if
         rrate_fiss(fissindex)%group    = group_index
         rrate_fiss(fissindex)%parts    = parts_index
         rrate_fiss(fissindex)%source   = src
         rrate_fiss(fissindex)%reac_src = rrs_fiss
         rrate_fiss(fissindex)%q_value  = q    !  average Q-value weighted over all channels
         rrate_fiss(fissindex)%param    = params

         fissindex = fissindex + 1
      end do fission_loop

      ! Clean up
      ! Close the neutron and betafission files again
      if (neutfission.gt.0) call close_io_file(neutfission,nfission_file)
      if (betafission.gt.0) call close_io_file(betafission,bfission_file)
      ! Deallocate the fission fragment arrays in case of fissflag 3
      if (allocated(fissfrags)) then
         deallocate(fissfrags,stat=astat)
         if (astat.ne.0) call raise_exception("Could not deallocate fissfrags array.",&
                                              "read_fission_rates",190002)
      end if
      ! Close the fission file again
      close(fissionlib)

   end subroutine read_fission_rates




   !> Determines fission fragment mass distribution as described in Panov et al. 2001.
   !!
   !! This routine is called for \ref parameter_class::fissflag = 1 and fills
   !! the array \ref fissrates with indices of the fragment nuclei.
   !!
   !! @see [Panov et al., Nuc. Phys. A688 2001](https://ui.adsabs.harvard.edu/abs/2001NuPhA.688..587P/abstract)
   subroutine fiss_dist(pos,mass,pnr,src,qval)
      use global_class,only:isotope,ineu
      use benam_class, only: minmax,findaz

      integer,intent(in)            :: pos
      integer,intent(in)            :: mass, pnr !< mass and proton number of fissioning nucleus
      character(4),intent(in)       :: src
      real(r_kind),intent(out)      :: qval      !< Q-value of the fission reaction
      integer                       :: af,zf     !< mass and proton number of "compound" nucleus
      integer                       :: a1,z1     !< mass and proton number of fragment 1
      integer                       :: a2,z2     !< mass and proton number of fragment 2
      integer                       :: nemiss    !< number of fission neutrons
      integer                       :: fiss_mode !< specifies fission mode

      INFO_ENTRY("fiss_dist")

      !neutron induced fission
      if (fissrate(pos)%reac_type .eq. rrt_nf) then
         fiss_mode = 1
         af = mass+1
         zf = pnr
      !spontaneous fission
      else if (fissrate(pos)%reac_type .eq. rrt_sf) then
         fiss_mode = 2
         af = mass
         zf = pnr
        !beta delayed fission
      else if (fissrate(pos)%reac_type .eq. rrt_bf) then
         fiss_mode = 3
         af = mass
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

      nemiss = 0

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

      fissrate(pos)%fissnuc_index  = findaz(mass,pnr)
      fissrate(pos)%channels       = 1
      allocate(fissrate(pos)%channelprob(1))
      allocate(fissrate(pos)%q_value(1))
      fissrate(pos)%channelprob(1) = 1.d0
      fissrate(pos)%mode           = fiss_mode

      if (fiss_mode.eq.1) then
         fissrate(pos)%dimens = 4
         allocate(fissrate(pos)%fissparts(fissrate(pos)%dimens))
         allocate(fissrate(pos)%ch_amount(fissrate(pos)%dimens))
         fissrate(pos)%fissparts(1) = ineu
         fissrate(pos)%ch_amount(1) = float(nemiss) - 1.d0
         fissrate(pos)%fissparts(2) = fissrate(pos)%fissnuc_index
         fissrate(pos)%ch_amount(2) = -1.d0
         fissrate(pos)%fissparts(3) = findaz(a1,z1)
         fissrate(pos)%ch_amount(3) = 1.d0
         fissrate(pos)%fissparts(4) = findaz(a2,z2)
         fissrate(pos)%ch_amount(4) = 1.d0
         fissrate(pos)%q_value(1)   = (1-nemiss) * isotope(ineu)%mass_exc + &
                                      isotope(fissrate(pos)%fissnuc_index)%mass_exc - &
                                      isotope(fissrate(pos)%fissparts(3))%mass_exc - &
                                      isotope(fissrate(pos)%fissparts(4))%mass_exc
      else                                  ! fiss_mode 2 and 3
         if (nemiss .eq. 0) then            ! no fission neutrons
            fissrate(pos)%dimens = 3
            allocate(fissrate(pos)%fissparts(fissrate(pos)%dimens))
            allocate(fissrate(pos)%ch_amount(fissrate(pos)%dimens))
         else                               ! fission neutrons are produced
            fissrate(pos)%dimens = 4
            allocate(fissrate(pos)%fissparts(fissrate(pos)%dimens))
            allocate(fissrate(pos)%ch_amount(fissrate(pos)%dimens))
            fissrate(pos)%fissparts(4) = ineu
            fissrate(pos)%ch_amount(4) = float(nemiss)
         end if
         fissrate(pos)%fissparts(1) = fissrate(pos)%fissnuc_index
         fissrate(pos)%ch_amount(1) = -1.d0
         fissrate(pos)%fissparts(2) = findaz(a1,z1)
         fissrate(pos)%ch_amount(2) = 1.d0
         fissrate(pos)%fissparts(3) = findaz(a2,z2)
         fissrate(pos)%ch_amount(3) = 1.d0
         fissrate(pos)%q_value(1)   = isotope(fissrate(pos)%fissparts(1))%mass_exc - &
                                      isotope(fissrate(pos)%fissparts(2))%mass_exc - &
                                      isotope(fissrate(pos)%fissparts(3))%mass_exc - &
                                      nemiss * isotope(ineu)%mass_exc
      end if

      qval = fissrate(pos)%q_value(1)
      fissrate(pos)%averageQ = qval
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
   !! .
   subroutine kodtakdist(pos,mass,pnr,src,qval)
      use parameter_class, only:unit
      use global_class, only: isotope,ineu
      use benam_class, only: minmax,findaz
      use error_msg_class, only: raise_exception

      integer,intent(in)       :: mass, pnr
      character(4),intent(in)  :: src
      integer,intent(in)       :: pos
      real(r_kind),intent(out) :: qval     !< Q-value of the fission reaction (averaged and weighted over all channels)

      integer                 :: nf        !< number of fission channels
      integer                 :: fiss_mode
      integer                 :: af,zf
      integer                 :: a,z
      integer                 :: a1,a2,z1,z2
      real(r_kind)            :: paz,ptot
      real(r_kind)            :: za,al,ah,cz,ca
      real(r_kind),parameter  :: lim = 1.d-6
      integer                 :: nemiss
      integer                 :: i
      integer                 :: dimens
      integer                 :: neutronflag

      INFO_ENTRY("kodtakdist")

      !neutron induced fission
      if (fissrate(pos)%reac_type .eq. rrt_nf) then
        fiss_mode = 1
        af = mass+1
        zf = pnr
      !spontaneous fission
      else if (fissrate(pos)%reac_type .eq. rrt_sf) then
        fiss_mode = 2
        af = mass
        zf = pnr
      !beta delayed fission
      else if (fissrate(pos)%reac_type .eq. rrt_bf) then
        fiss_mode = 3
        af = mass
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

               neutronflag = 0

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

      fissrate(pos)%channels       = nf
      fissrate(pos)%fissnuc_index  = findaz(mass,pnr)
      fissrate(pos)%mode           = fiss_mode

      allocate(fissrate(pos)%channelprob(nf))
      allocate(fissrate(pos)%q_value(nf))

      select case (fiss_mode)
      case(1)
         dimens = 2 * fissrate(pos)%channels + 2
         fissrate(pos)%dimens = dimens
         allocate(fissrate(pos)%fissparts(dimens))
         allocate(fissrate(pos)%ch_amount(dimens))
         fissrate(pos)%fissparts(1) = ineu
         fissrate(pos)%fissparts(2) = fissrate(pos)%fissnuc_index
         fissrate(pos)%ch_amount(2) = -1.d0
      case(2:)
         dimens = 2 * fissrate(pos)%channels + 1
         if (neutronflag.eq.1) dimens = dimens+1
         fissrate(pos)%dimens = dimens
         allocate(fissrate(pos)%fissparts(dimens))
         allocate(fissrate(pos)%ch_amount(dimens))
         fissrate(pos)%ch_amount(1) = -1.d0
         fissrate(pos)%fissparts(1) = fissrate(pos)%fissnuc_index
      end select

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
               if (VERBOSE_LEVEL .ge. 2) then
                print *, 'Warning in kodtakdist, fragment was lighter than included nuclei'
                print *, 'Fragment was: ', a2, z2
               end if
            end if

   ! fill in fissrate()%fissparts(); same as rrate()%parts(1:6), but with individual array sizes
            select case(fiss_mode)
            case(1)
               fissrate(pos)%channelprob(i) = paz
               fissrate(pos)%ch_amount(1)   = float(nemiss) - 1.d0                  ! TODO: replace with average value (?)
               fissrate(pos)%fissparts(i+2) = findaz(a1,z1)
               fissrate(pos)%ch_amount(i+2) = fissrate(pos)%channelprob(i)   ! rate at which fragment is produced per destroyed parent nucleus
               fissrate(pos)%fissparts(i+2+fissrate(pos)%channels) = findaz(a2,z2)
               fissrate(pos)%ch_amount(i+2+fissrate(pos)%channels) = fissrate(pos)%channelprob(i)
               fissrate(pos)%q_value(i) = isotope(fissrate(pos)%fissnuc_index)%mass_exc - &
                       (nemiss - 1) * isotope(ineu)%mass_exc - &                   ! (nemiss - 1) to account for one neutron that is destroyed
                       isotope(fissrate(pos)%fissparts(i+2))%mass_exc - &
                       isotope(fissrate(pos)%fissparts(i+2+fissrate(pos)%channels))%mass_exc
            case(2:)
               fissrate(pos)%channelprob(i) = paz
               fissrate(pos)%fissparts(i+1) = findaz(a1,z1)
               fissrate(pos)%ch_amount(i+1) = fissrate(pos)%channelprob(i)    ! rate at which fragment is produced per destroyed parent nucleus
               fissrate(pos)%fissparts(i+1+fissrate(pos)%channels) = findaz(a2,z2)
               fissrate(pos)%ch_amount(i+1+fissrate(pos)%channels) = fissrate(pos)%channelprob(i)
               if (neutronflag.eq.1) then
                  fissrate(pos)%fissparts(dimens) = ineu
                  fissrate(pos)%ch_amount(dimens) = real(nemiss)
               end if
               fissrate(pos)%q_value(i) = isotope(fissrate(pos)%fissnuc_index)%mass_exc - &
                       nemiss * isotope(ineu)%mass_exc - &
                       isotope(fissrate(pos)%fissparts(i+1))%mass_exc - &
                       isotope(fissrate(pos)%fissparts(i+1+fissrate(pos)%channels))%mass_exc
            end select
         end do massloop
      end do

      qval = 0.d0
      qvalue: do i=1,fissrate(pos)%channels
         qval = qval + fissrate(pos)%q_value(i)*fissrate(pos)%channelprob(i)     ! weighted average Q-value
      end do qvalue
      fissrate(pos)%averageQ = qval

      nufiss = nufiss + int(paz*float(nemiss))

      INFO_EXIT("kodtakdist")
      return

   end subroutine kodtakdist



   !> Read the fission distribution from a file
   !!
   !! The fission distribution is read from a file and stored in the
   !! fissfrags array. The default file contains the fission distribution
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
   !! parent nucleus. Furthermore,  \f$\sum Y(Z,A) = 2\f$.
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
   subroutine read_mumpower_fissfile(file)
    use file_handling_class
    use benam_class, only: findaz, minmax
    implicit none

    character(len=*), intent(in) :: file !< File path to mumpower file

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


    INFO_ENTRY("read_mumpower_fissfile")

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
    allocate(fissfrags(frag_counter),stat=astat)
    if (astat .ne. 0) then
        call raise_exception("Could not allocate 'fissfrags'.", &
                             "read_mumpower_fissfile",190001)
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
        allocate(fissfrags(i)%net_idx(nof),&
                 fissfrags(i)%Z(nof),&
                 fissfrags(i)%A(nof),&
                 fissfrags(i)%Y(nof),&
                 stat=astat)
        if (astat .ne. 0) then
            call raise_exception("Could not allocate fragment properties.", &
                                 "read_mumpower_fissfile",190001)
        end if

        ! Store the properties
        fissfrags(i)%nr_frags = nof
        fissfrags(i)%Zp = zf
        fissfrags(i)%Ap = af

        ! Count the neutrons emitted
        nem = 0
        Ynem = 0

        ! Read the fragments
        do j=1,nof
            read(file_id,*) fissfrags(i)%Z(j),fissfrags(i)%A(j),fissfrags(i)%Y(j)
            fissfrags(i)%net_idx(j) =findaz(fissfrags(i)%A(j),fissfrags(i)%Z(j))

            ! What to do when the fragment is not in the network?
            ! Ideally make neutrons until the fragment is in the network
            if (fissfrags(i)%net_idx(j) .eq. -1) then
                ! Tell the world that this is done
                if (VERBOSE_LEVEL .ge. 2) then
                    write(*,*) "Fragment ",fissfrags(i)%Z(j),fissfrags(i)%A(j),&
                               " not in network. Splitting into lighter isotope and neutrons."
                end if

                mina = minmax(fissfrags(i)%Z(j),1)
                maxa = minmax(fissfrags(i)%Z(j),2)
                if (fissfrags(i)%A(j) .gt. maxa) then
                    ! Rest goes to neutrons
                    nem = nem + (fissfrags(i)%A(j)-maxa)
                    ! Yield of neutrons is the amount scaled with the
                    ! yield of the not included fragment
                    ! Y * A should be conserved
                    Ynem = Ynem + fissfrags(i)%Y(j)*(fissfrags(i)%A(j)-maxa)
                    fissfrags(i)%Y(j) = fissfrags(i)%Y(j) ! THis stays conserved
                    fissfrags(i)%A(j) = maxa
                    fissfrags(i)%net_idx(j) = findaz(maxa,fissfrags(i)%Z(j))
                end if
            end if
        end do

        ! Check how many neutrons are avaiable and give it to the missing
        ! fragments
        do j=1,nof
            if (fissfrags(i)%net_idx(j) .eq. -1) then

                ! Say something if verbose
                if (VERBOSE_LEVEL .ge. 2) then
                    write(*,*) "Fission fragment Z=",fissfrags(i)%Z(j), &
                               "A=",fissfrags(i)%A(j), &
                               "Y=",fissfrags(i)%Y(j), &
                               "not in network, trying to remap."
                end if

                ! Needed neutrons are:
                mina = minmax(fissfrags(i)%Z(j),1)
                bn = (mina-fissfrags(i)%A(j)) ! Amount of missing neutrons

                ! Borrow neutrons, take heaviest fragment that has enough
                ! abundance
                do m=nof,1,-1
                    if ((fissfrags(i)%net_idx(m) .ne. -1) .and. &
                        ((bn*fissfrags(i)%Y(j) / fissfrags(i)%A(m)) .lt. fissfrags(i)%Y(m))) then
                        exit
                    end if
                end do

                ! None found? Too bad, make an error
                if (m .eq. 0) then
                    call raise_exception("Could not remap fission fragments, "// &
                                         "consider to include more nuclides in the network. "//&
                                         NEW_LINE('A')//&
                                         "Nucleus Z="//trim(int_to_str(fissfrags(i)%Z(j)))// &
                                         ", A="//trim(int_to_str(fissfrags(i)%A(j)))// &
                                         " was not included. "//NEW_LINE('A')//&
                                         "Need "//int_to_str(bn)//" neutrons.",&
                                         "read_mumpower_fissfile",&
                                         190008)
                end if

                ! Otherwise borrow the neutrons from the heavier fragment
                fissfrags(i)%Y(m) = fissfrags(i)%Y(m) - (bn*fissfrags(i)%Y(j) / fissfrags(i)%A(m))
                fissfrags(i)%Y(j) = fissfrags(i)%Y(j)
                ! Give the fragment a new identity
                fissfrags(i)%A(j) = mina
                fissfrags(i)%net_idx(j) = findaz(mina,fissfrags(i)%Z(j))
            end if
        end do

        ! Proof the yields, this should be Af!
        helper = Ynem
        do j=1,nof
            helper = helper+fissfrags(i)%A(j)*fissfrags(i)%Y(j)
        end do
        if (abs(helper - float(fissfrags(i)%Ap)) .gt. 1e-2) then
            call raise_exception("Yield of fission fragments is not conserved. "// &
                                 "This happened for fragment distribution of "// &
                                 "Nucleus Z="//trim(int_to_str(fissfrags(i)%Zp))// &
                                 ", A="//trim(int_to_str(fissfrags(i)%Ap))// &
                                 ". Ensure that the sum is 2.","read_mumpower_fissfile",&
                                 190009)
        end if

        fissfrags(i)%neutrons = nem
        fissfrags(i)%Yn       = Ynem
        nufiss                = nufiss + nem
    end do fragloop

    ! Close the file
    close(file_id)

    INFO_EXIT("read_mumpower_fissfile")

   end subroutine read_mumpower_fissfile




   !> Fill the rates with the correct fragments
   !!
   !! This subroutine fills the rates with the correct fragments
   !! from [Mumpower et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020PhRvC.101e4607M/abstract).
   !! This subroutine makes use of the fragtype struct.
   !! If no fragment is found, the fragment distribution is set to the
   !! one of [Panov et al. (2001)](https://ui.adsabs.harvard.edu/abs/2001NuPhA.688..587P/abstract).
   !! For beta-delayed fission the fissioning nucleus is the one of (Z+1, A)
   !! for neutron induced fission the nucleus with (Z, A+1), and for spontanous
   !! fission (not used within fissflag = 3) it is (Z,A).
   !!
   !! @author M. Reichert
   !! @date 14.02.2023
   subroutine mumpower_fiss(pos,mass,pnr,src,qval)
    use global_class, only: isotope,ineu
    use benam_class, only: minmax,findaz
    use error_msg_class, only: raise_exception
       integer,intent(in)                     :: mass        !< Mass number of reacting nucleus
       integer,intent(in)                     :: pnr         !< Proton number of reacting nucleus
       integer, intent(in)                    :: pos         !< Position of the reaction
       real(r_kind),intent(out)               :: qval        !< Q-value of the reaction
       character(4), intent(in)               :: src         !< Source of the reaction
       ! Internal variables
       integer                                :: nfrac_distr !< Number of fragment distributions
       integer                                :: i           !< Loop variable
       integer                                :: astat       !< Allocation status variable
       integer                                :: ind_parent  !< Index of reacting nucleus
       integer                                :: afiss       !< A of fissioning nucleus
       integer                                :: zfiss       !< Z of fissioning nucleus
       type(fragtype)                         :: fissfrags_tmp!< Temporary fragment distribution
       logical                                :: found       !< Helper variable flag if the fragment
                                                             !  distribution was found

       INFO_ENTRY("mumpower_fiss")

      ! Loop through the fission fragments
       nfrac_distr = size(fissfrags)

       ! Check which nuclei is fissioning

       zfiss = pnr
       afiss = mass
       ! Z-1 for beta delayed fission
       if (fissrate(pos)%reac_type .eq. rrt_bf) zfiss = zfiss+1
       ! A+1 for neutron induced fission
       if (fissrate(pos)%reac_type .eq. rrt_nf) afiss = afiss+1
       ! Don't change anything for spontaneous fission

       ! Find relevant indices
       ind_parent = findaz(mass,pnr)
       fissrate(pos)%fissnuc_index = ind_parent
       found = .false.
        do i=1,nfrac_distr
            if ((fissfrags(i)%Zp .eq. zfiss) .and. (fissfrags(i)%Ap .eq. afiss)) then
                found = .true.
                exit
            end if
        end do

        ! Not found, do panof in this case
        if (.not. found) then
            ! Output this information
            if (VERBOSE_LEVEL .ge.2) then
                write(*,*) "No fragment distribution found for Z=",zfiss,", A=",afiss
                write(*,*) "Using Kodama & Takahashi 1975 instead."
            end if
            ! Get panov distribution
            call kodtakdist(pos,mass,pnr,src,qval)
        else
            ! Found, do mumpower
            fissfrags_tmp = fissfrags(i)

            ! Allocate the parts, the first two are the fissioning nucleus and neutrons
            allocate(fissrate(pos)%fissparts(fissfrags_tmp%nr_frags+2),&
                     fissrate(pos)%ch_amount(fissfrags_tmp%nr_frags+2),&
                     stat=astat)
            ! Complain if not possible
            if (astat .ne. 0) then
                call raise_exception("Could not allocate memory for fission fragments.",&
                                     "mumpower_fiss",190001)
            end if

            ! Set up the number of participants,
            ! Number of fragments + parent nucleus + neutrons
            fissrate(pos)%dimens = fissfrags_tmp%nr_frags+2

            ! Use the correct order of the reacting nuclei
            if (.not. (fissrate(pos)%reac_type .eq. rrt_nf)) then
                fissrate(pos)%fissparts(1) = ind_parent
                fissrate(pos)%fissparts(2) = ineu
                fissrate(pos)%ch_amount(1) = -1d0
                fissrate(pos)%ch_amount(2) = fissfrags_tmp%Yn
            else
                fissrate(pos)%fissparts(1) = ineu
                fissrate(pos)%fissparts(2) = ind_parent
                fissrate(pos)%ch_amount(1) = fissfrags_tmp%Yn-1
                fissrate(pos)%ch_amount(2) = -1d0
            end if


            do i=1,fissfrags_tmp%nr_frags
                fissrate(pos)%fissparts(i+2) = fissfrags_tmp%net_idx(i)
                fissrate(pos)%fissparts(i+2) = fissfrags_tmp%net_idx(i)
                fissrate(pos)%ch_amount(i+2) = fissfrags_tmp%Y(i)
            end do
        end if

        ! Calculate the Q-value
        qval = 0.d0
        qvalue: do i=1,fissrate(pos)%dimens
           qval = qval - isotope(fissrate(pos)%fissparts(i))%mass_exc*fissrate(pos)%ch_amount(i) ! weighted average Q-value
        end do qvalue
        fissrate(pos)%averageQ = qval

       INFO_EXIT("mumpower_fiss")

    end subroutine mumpower_fiss



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
   subroutine abla_betafiss(pos,mass,pnr,src,betafission,qval)
      use global_class,    only: isotope,ineu
      use benam_class,     only: minmax,findaz
      use error_msg_class, only: raise_exception
         integer,intent(in)                     :: mass,pnr
         character(4), intent(in)               :: src
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
             if (src.eq.'sfis') then
                fissrate(pos)%mode = 2                ! spontaneous fission
             else
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
         if (src.eq.'sfis') then      ! spontaneous fission
            fpnr = pnr
            fissrate(pos)%mode = 2
         else                         ! beta-delayed fission
            fpnr  = pnr + 1
            fissrate(pos)%mode = 3
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
