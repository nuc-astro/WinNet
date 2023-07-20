!> @file beta_decay_rate_module.f90
!!
!! The error file code for this file is ***W12***.
!! Contains the module \ref beta_decay_rate_module.


!>
!! @brief This module contains subroutines to include external beta decays
!!
!! The external beta decays are given in a different format compared
!! to reaclib. They tabulate the halflife of a nucleus together
!! with Pn probabilities (beta delayed neutron emission). The file
!! includes up to P10n (10 beta delayed neutrons), which can not be
!! calculated in the default reaclib format (There P2n is max. in chapter 3,
!! or P3n in the new reaclib format in chapter 11). An example entry looks like:
!!\file{...
!! ca73    6.200000e-04
!! 0.0000   0.0100   0.0200   0.0100   0.2800   0.0800   0.4000   0.0700   0.1100   0.0100   0.0100
!! ...}
!!
!! Some tables, however, tabulate until P10n, for example,
!! [Moeller et al. 2019](https://matthewmumpower.com/publications/paper/2019/moller/nuclear-properties-for-astrophysical-and-radioactive-ion-beam-applications-ii)
!!
!! @author Moritz Reichert
!! @date 24.01.21
#include "macros.h"
module beta_decay_rate_module
   use global_class, only: reactionrate_type
   implicit none

   !> Type for storing the beta decay data with Pn channels
   type,private                :: beta_pn_type
    integer                    :: idx_parent     !< Index of parent in abundance array
    integer,dimension(11)      :: idx_daughter   !< Index of daughter in abundance array
    integer                    :: Z_parent       !< Atomic number of parent
    integer                    :: N_parent       !< Neutron number of parent
    real(r_kind),dimension(11) :: Pn             !< Pn probabilities
    real(r_kind),dimension(11) :: Qval           !< Qvalue of the channels
    real(r_kind)               :: Av_Qtot        !< Average Q-Value
    real(r_kind)               :: Av_Qnu         !< Average Q-Value of neutrinos
    real(r_kind)               :: nu_frac        !< Energy fraction of neutrinos
    real(r_kind)               :: halflife       !< halflife of parent (s)
   end type beta_pn_type
   type(beta_pn_type),dimension(:),allocatable,private   :: beta_pn  !< Array storing the reaction rates
   integer,private                                       :: nbeta_pn !< Number of beta decays



   type(reactionrate_type),dimension(:),allocatable,public :: beta_decays !< array containing external beta decays
   integer,public                                          :: nbeta       !< total number of external beta decays
   logical,private                                         :: ext_decays  !< Flag if external beta decays are used
   ! Helper variables for ignoring specific sources
   character(len=4),allocatable,dimension(:),private :: src_ignore
   integer,private                                   :: src_ignore_length

   !
   ! Public and private fields and methods of the module
   !
   public:: &
       init_ext_beta_rates, merge_beta_decays
   private:: &
       count_reactions, remove_weak_rates, read_beta_decays

contains

   !> Initialize external beta decay rates.
   !!
   !! This subroutine counts and reads the rates into the array
   !! \ref beta_decays.
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine init_ext_beta_rates()
      use parameter_class, only: beta_decay_file,use_beta_decay_file,&
                                 beta_decay_src_ignore
      use nucstuff_class,  only: analyze_src_string
      use file_handling_class
      implicit none
      integer     :: beta_unit  !< File ID of external beta decay file
      integer     :: alloc_stat !< Allocation status

      ! External beta decay flag
      ext_decays = .False.

      if (use_beta_decay_file) then
         ! Flag that external beta decays are used
         ext_decays = .True.

         ! Check if some sources should be ignored
         call analyze_src_string(beta_decay_src_ignore,src_ignore,src_ignore_length)

         beta_unit= open_infile(beta_decay_file)
         call count_reactions(beta_unit)

         !----- allocate the array of beta decays
         allocate(beta_pn(nbeta_pn),stat=alloc_stat)
         if ( alloc_stat /= 0)  call raise_exception('Allocation of "beta_pn" failed',&
                                                     "init_ext_beta_rates",&
                                                     120001)


         ! Read the beta decays
         call read_beta_decays(beta_unit)
         ! Close the file again
         call close_io_file(beta_unit,beta_decay_file)

        !  ! Say how many there were
         if (VERBOSE_LEVEL .ge. 1) then
            call write_data_to_std_out("Amount beta-decay format rates",int_to_str(nbeta_pn*11))
         end if
      end if

   end subroutine init_ext_beta_rates


   !> Merge external beta decays into the larger rate array.
   !!
   !! The return value of this routine will be a larger rate array
   !! and a new length
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine merge_beta_decays(rrate_array,rrate_length)
      use error_msg_class,  only: raise_exception
      use mergesort_module, only: rrate_ms,rrate_sort
      implicit none
      type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate_array  !< Large rate array, containing all reactions
      integer,intent(inout)                                          :: rrate_length !< length of rrate_array
      integer                                                        :: alloc_stat   !< Allocation state



      if (ext_decays) then
        if (nbeta_pn .ne. 0) then
           if (.not. allocated(rrate_array)) then
              ! Create an array in the correct format
              call create_rrate_array(beta_decays,nbeta)

              rrate_length = nbeta
              !-- Allocate the reaclib rate array
              allocate(rrate_array(nbeta),stat=alloc_stat)
              if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                         "merge_beta_decays",&
                                                         120001)
              rrate_array(1:nbeta) = beta_decays(1:nbeta)
           else
              ! Check if we should not put some reactions into the array
              call ignore_reactions(rrate_array,rrate_length)
              !----- merge weak rates into rrate
              call remove_weak_rates(rrate_array,rrate_length)
           end if
           !-- Deallocate the reaclib rate array
           deallocate(beta_decays,stat=alloc_stat)
           if (alloc_stat .ne. 0) call raise_exception('Deallocation of "beta_decays" failed',&
                                                       "merge_beta_decays", 120002)
        end if
      end if

   end subroutine merge_beta_decays


   !> Subroutine to remove weak rates from the beta decay rate array
   !!
   !! This subroutine ensures that \ref parameter_class::beta_decay_src_ignore
   !! is working. It removes all reactions from the array \ref beta_pn
   !! that are listed in \ref parameter_class::beta_decay_src_ignore.
   !!
   !! @author M. Reichert
   subroutine ignore_reactions(rrate_array,rrate_length)
    use global_class,    only: isotope
    use error_msg_class, only: int_to_str, raise_exception
    implicit none
    integer,intent(in)                               :: rrate_length
    type(reactionrate_type),dimension(:),allocatable :: rrate_array
    type(beta_pn_type),dimension(:),allocatable      :: rrate_tmp

    integer :: i,j
    integer :: idx
    integer :: z_tmp,n_tmp
    integer :: replace_count
    logical,dimension(:),allocatable :: mask
    integer :: istat

    ! Check if you really want to add the rate
    if (src_ignore_length .gt. 0) then
        replace_count = 0
        allocate(mask(nbeta_pn),stat=istat)
        if (istat .ne. 0) call raise_exception('Allocation of "mask" failed',"ignore_reactions",&
                                              120001)
        mask(:) = .true.

        do i=1,nbeta_pn
            ! Get the parent again and compare,
            ! note that this could be vastly improved...
            do j=1,rrate_length
                if ((.not. rrate_array(j)%is_weak) .or. &
                    (.not. rrate_array(j)%reac_type .eq. rrt_betm) .or. &
                    (.not. rrate_array(j)%group .eq. 1)) cycle
                if (rrate_array(j)%parts(1) .ne. beta_pn(i)%idx_parent) cycle

                ! Found the chapter 1 rate, is it included in the src_ignore list?
                if (any(src_ignore .eq. adjustr(rrate_array(j)%source))) then
                    replace_count = replace_count+1
                    mask(i) = .false.
                end if
            end do
        end do

        if (replace_count .gt. 0) then

            ! Say something
            if (VERBOSE_LEVEL .ge. 2) then
                print*,"Ignoring "//int_to_str(replace_count)//" beta-decay rates by source criterium (beta_decay_src_ignore)."
            end if

            ! Now we have a mask, pack the array
            rrate_tmp = pack(beta_pn,mask)

            deallocate(beta_pn,stat=istat)
            if (istat .ne. 0) call raise_exception('Deallocation of "beta_pn" failed',"ignore_reactions",&
                                                120002)
            ! We removed rates so adjust the size
            nbeta_pn = nbeta_pn-replace_count

            ! Allocate beta_decays with new length
            allocate(beta_pn(nbeta_pn),stat=istat)
            if (istat .ne. 0) call raise_exception('Allocation of "beta_pn" failed',"ignore_reactions",&
                                                120001)
            ! Copy the temporary array back
            beta_pn(1:nbeta_pn) = rrate_tmp(1:nbeta_pn)

            ! Deallocate the temporary array again
            deallocate(rrate_tmp,stat=istat)
            if (istat .ne. 0) call raise_exception('Deallocation of "rrate_tmp" failed',"ignore_reactions",&
                                                120002)
        else
            ! Say something
            if (VERBOSE_LEVEL .ge. 2) then
                print*,"Warning: Specified 'beta_decay_src_ignore', but no rate ignored."
            end if
        end if

    end if


   end subroutine ignore_reactions



   !> Count the amount of external beta decays.
   !!
   !! After this routine has been called, nbeta will contain the amount of
   !! external beta decay rates.
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine count_reactions(sourcefile)
      use global_class, only: isotope
      use benam_class,  only: findaz, benam
      implicit none

      integer, intent(in)         :: sourcefile  !< file id to read beta decays

      real(r_kind),dimension(11)  :: Pn    !< P0 to P10 channel probabilities
      character(5)            :: parent    !< Name of the parent nucleus
      real(r_kind)            :: halflife  !< halflife of the parent nucleus
      integer                 :: count     !< Counter for the reaction rates
      integer                 :: read_stat !< Read status in the file
      integer                 :: index_tmp !< Temporary index of the nucleus in net_names
      integer                 :: z_tmp     !< Temporary storage for amount of protons
      integer                 :: n_tmp     !< Temporary storage for amount of neutrons
      integer                 :: a_tmp     !< Temporary storage for amount of nucleons
      integer                 :: j         !< Loop variable


      count = 0
      ! Count the number of reactions
      beta_loop: do
         read(sourcefile,*, iostat = read_stat) parent,halflife
         read(sourcefile,*, iostat = read_stat) Pn(1:11)
         index_tmp = benam(parent)
         if (read_stat /= 0) exit
         if (index_tmp .eq. 0) cycle

         ! Ensure that at least one daughter is
         ! also in the table
         ! Find out decay products
         z_tmp = isotope(index_tmp)%p_nr
         n_tmp = isotope(index_tmp)%n_nr
         a_tmp = z_tmp + n_tmp
         do j=1,11
            ! Daughter index
            index_tmp = findaz(a_tmp-(j-1),z_tmp+1)
            if (index_tmp .ne. -1) exit
         end do
         ! The nucleus was found
         if (j .ne. 12) count = count + 1

      end do beta_loop

      nbeta_pn = count

   end subroutine count_reactions


   !> Remove beta decays from the rate array
   !!
   !! This routine removes beta decays in case that they exist in
   !! another format ( see parameter_class::use_beta_decay_file ).
   !!
   !! @author Moritz Reichert
   subroutine remove_weak_rates(rrate_array,length_rate_array)
     use global_class,    only: reactionrate_type,ineu
     use error_msg_class, only: raise_exception
     implicit none

     logical,dimension(:),allocatable      :: mask,mask_bet
     logical :: delete
     integer  :: replace_count
     integer  :: istat
     integer  :: i,j,k,m
     integer  :: length_rate_array
     integer  :: nbeta_clean
     type(reactionrate_type),dimension(:),allocatable :: rrate_tmp
     type(reactionrate_type),dimension(:),allocatable :: rrate_array

        INFO_ENTRY("remove_weak_rates")
        replace_count = 0
        allocate(mask(length_rate_array),stat=istat)
        if (istat .ne. 0) call raise_exception('Allocation of "mask" failed',"remove_weak_rates",&
                                              120001)
        mask(:) = .true.

        do i = 1 , length_rate_array
           if ((.not. rrate_array(i)%is_weak) .or. &
               (.not. rrate_array(i)%reac_type .eq. rrt_betm)) cycle
           do j=1, nbeta_pn
              if (beta_pn(j)%idx_parent .eq. rrate_array(i)%parts(1)) then

                part_loop: do k=2,6
                  if (rrate_array(i)%parts(k) .eq. 0) exit part_loop
                  delete = .false.

                  parts_bet: do m=1,11
                    ! Daughter included?
                    if (beta_pn(j)%idx_daughter(m) .eq. -1) cycle parts_bet

                    if  ((beta_pn(j)%idx_daughter(m) .eq. rrate_array(i)%parts(k)) &
                        .or. (rrate_array(i)%parts(k) .eq. ineu)) then
                        delete = .true.
                        exit parts_bet
                     end if

                  end do parts_bet
                  ! Was something else than all possible daughters or neutrons?
                  if (delete .eqv. .false.) exit part_loop
                end do part_loop

                if (delete) then
                  mask(i) = .false.
                  replace_count = replace_count + 1
                end if
              end if
           end do
        end do

        if (VERBOSE_LEVEL.ge.2) then
           print*,'Replacing ',replace_count,' weak rates!'
        end if

        ! Remove the rates and implement other decay rates
        if ((replace_count .gt. 0) .or. (nbeta_pn .gt. 0))  then

           ! Get the beta decay rates in reactionrate_type format
           call create_rrate_array(beta_decays,nbeta)

           ! Deallocate the old array
           deallocate(beta_pn,stat=istat)
           if (istat .ne. 0) call raise_exception('Deallocation of "beta_pn" failed',"remove_weak_rates",&
                                                  120002)

           ! Now merge the two arrays
           allocate(rrate_tmp(length_rate_array-replace_count),stat=istat)
           if (istat .ne. 0) call raise_exception('Allocation of "rrate_tmp" failed',"remove_weak_rates",&
                                                  120001)

           rrate_tmp = pack(rrate_array,mask)

           deallocate(rrate_array,stat=istat)
           if (istat .ne. 0) call raise_exception('Deallocation of "rrate" failed',"remove_weak_rates",&
                                                  120002)
           length_rate_array = length_rate_array-replace_count

           allocate(rrate_array(length_rate_array+nbeta),stat=istat)
           if (istat .ne. 0) call raise_exception('Allocation of "rrate" failed',"remove_weak_rates",&
                                                  120001)

           rrate_array(1:length_rate_array) = rrate_tmp(1:length_rate_array)
           rrate_array(length_rate_array+1:length_rate_array+nbeta) = beta_decays(1:nbeta)
           length_rate_array = length_rate_array + nbeta

           deallocate(rrate_tmp,stat=istat)
           if (istat .ne. 0) call raise_exception('Deallocation of "rrate_tmp" failed',"remove_weak_rates",&
                                                  120002)
        end if

        INFO_EXIT("remove_weak_rates")
   end subroutine remove_weak_rates



   !> Convert \ref beta_pn_type to \ref reactionrate_type
   !!
   !! Creates a reactionrate_type array from the beta_pn_type array.
   !! Rates with Pn of 0 will not result in an extra rate.
   !!
   !! @author Moritz Reichert
   !! @date 28.02.2023
   subroutine create_rrate_array(rrate_beta,rrate_beta_length)
    use global_class,    only: reactionrate_type, ineu, Qnuloss
    use parameter_class, only: heating_frac, use_neutrino_loss_file
    use error_msg_class, only: raise_exception
    implicit none
    type(reactionrate_type), allocatable, intent(out) :: rrate_beta(:)     !< Reaction rates in reactionrate_type format
    integer, intent(out)                              :: rrate_beta_length !< Length of rrate_beta

    integer :: i,j   !< Loop variables
    integer :: istat !< Status variable
    integer :: count !< Counter for the reaction rates

    ! First count the amount of necessary beta decay rates
    count = 0
    do i=1,nbeta_pn
        pn_loop: do j=1,11
            if (beta_pn(i)%pn(j) .eq. 0) cycle pn_loop
            count = count + 1
        end do pn_loop
    end do

    rrate_beta_length = count

    ! Allocate the array
    allocate(rrate_beta(count),stat=istat)
    if (istat .ne. 0) call raise_exception('Allocation of "rrate_beta" failed',"create_rrate_array",&
                                           120001)

    ! Fill the array
    count = 0
    do i=1,nbeta_pn
        pn_loop2: do j=1,11
            if (beta_pn(i)%pn(j) .eq. 0) cycle pn_loop2
            count = count + 1
            rrate_beta(count)%parts(:)    = 0
            rrate_beta(count)%source      = "wext" ! weak-extern
            rrate_beta(count)%is_reverse  = .false.
            rrate_beta(count)%cached      = -1
            rrate_beta(count)%is_resonant = .false.
            rrate_beta(count)%is_weak     = .true.
            rrate_beta(count)%is_const    = .true.
            rrate_beta(count)%q_value     = beta_pn(i)%Qval(j)
            rrate_beta(count)%reac_src    = rrs_wext
            rrate_beta(count)%reac_type   = rrt_betm
            rrate_beta(count)%param(:)    = 0.0d0
            rrate_beta(count)%one_over_n_fac = 1.0d0

            ! Get the amount of energy radiated away by neutrinos
            if ((beta_pn(i)%Av_Qtot .eq. -1) .or. (beta_pn(i)%Av_Qnu .eq. -1)) then
                if (use_neutrino_loss_file) then
                    if (Qnuloss(beta_pn(i)%idx_parent) .ne. -1) then
                        rrate_beta(count)%nu_frac = Qnuloss(beta_pn(i)%idx_parent)/beta_pn(i)%Av_Qtot
                    else
                        rrate_beta(count)%nu_frac = heating_frac
                    end if
                else
                    rrate_beta(count)%nu_frac = heating_frac
                end if
            else
                rrate_beta(count)%nu_frac = beta_pn(i)%Av_Qnu/beta_pn(i)%Av_Qtot
            end if

            rrate_beta(count)%group        = min(j,2)   ! Put all remaining rates into chapter 2
            rrate_beta(count)%parts(1)     = beta_pn(i)%idx_parent ! parent nucleus
            rrate_beta(count)%ch_amount(1) = -1 ! Parent is destroyed

            if (j-1 .eq. 0) then
                rrate_beta(count)%parts(2) = beta_pn(i)%idx_daughter(j) ! decay product
                rrate_beta(count)%ch_amount(2) = 1                      ! product is created
            else
                rrate_beta(count)%parts(2) = ineu                       ! neutrons
                rrate_beta(count)%parts(3) = beta_pn(i)%idx_daughter(j) ! decay product
                rrate_beta(count)%ch_amount(2) = j-1                    ! j-1 nuclei are created
                rrate_beta(count)%ch_amount(3) = 1                      ! product is created
            end if

            ! Put the correct reaclib parameter
            rrate_beta(count)%param(1) =  dlog(beta_pn(i)%Pn(j)*(dLOG(2.0d0)/beta_pn(i)%halflife))

        end do pn_loop2
    end do

   end subroutine create_rrate_array





   !>
   !! Reads beta decays from a separate file.
   !!
   !! This routine is thought to read beta decays in a different format than reaclib.
   !! An advantage of this is that we can implement more than the P2n channel.
   !! The new format allows beta delayed neutrons up to 10 neutrons that become more
   !! important close to stability.
   !! An example input file could look like:
   !!\file{
   !! ca73    6.200000e-04
   !! 0.0000   0.0100   0.0200   0.0100   0.2800   0.0800   0.4000   0.0700   0.1100   0.0100   0.0100 }
   !! which will replace the decay of ca73 with a half life of 6.2e-04s and the listet Pn probabilities.
   !! Additionally, the average Q-Value and the average neutrino Q-Value can be given optionally in the first line.
   !! In this case an entry could look like:
   !! ca73    1.268195e-03  2.443690e+01  9.034672e+00
   !! 0.0030   0.0000   0.0010   0.0000   0.0020   0.9950   0.0000   0.0000   0.0000   0.0000   0.0000    }
   !! where the first line gives the name of the parent nucleus, the half life, the average Q-Value and the
   !! average neutrino Q-Value.
   !!
   !! @author Moritz Reichert
   !!
   !! \b Edited:
   !!         - 20.05.19
   !! .
   subroutine read_beta_decays(sourcefile)
      use global_class,    only: isotope,ineu
      use benam_class,     only: benam, findaz
      use error_msg_class, only: raise_exception, int_to_str
      use parameter_class, only: unit
      implicit none
      integer, intent(in)  :: sourcefile  !< file id to read beta decays
      !
      real(r_kind),dimension(11)  :: Pn     !< P0 to P10 channel probabilities
      character(5)            :: parent     !< Name of the parent nucleus
      real(r_kind)            :: halflife   !< halflife of the parent nucleus
      real(r_kind)            :: Qtot       !< Average total Q-Value
      real(r_kind)            :: Qnu        !< Average neutrino Q-Value
      real(r_kind)            :: probnot    !< Probability of not included daughters
      type(reactionrate_type) :: rr_tmp     !< temporary reaction rate
      integer                 :: count      !< Counter for the reaction rates
      integer                 :: read_stat  !< Read status in the file
      integer                 :: index_tmp  !< Temporary index of the nucleus in net_names
      integer                 :: z_tmp      !< Temporary storage for amount of protons
      integer                 :: n_tmp      !< Temporary storage for amount of neutrons
      integer                 :: a_tmp      !< Temporary storage for amount of nucleons
      character(200)          :: helper     !< helper variable to get the amount of columns in the file
      integer                 :: par_ind_tmp!< Index of parent nucleus in net_names
      integer                 :: col_count  !< Amount of columns in the file
      integer                 :: fformat    !< Format of the file
      integer                 :: j          !< Loop variable
      integer                 :: i          !< Loop variable

      INFO_ENTRY("read_beta_decays")

      rewind(sourcefile)

      ! Find out the number of columns in the first line of the file
      ! This defines the format
      read(sourcefile,"(A)", iostat = read_stat) helper
      col_count = 0
      do i = 2,200
         if ((helper(i-1:i-1) .ne. ' ') .and. (helper(i:i) .eq. ' ')) then
            ! Count the columns
            col_count = col_count +1
          end if
      end do

      ! Check if the file is compatible
      if (col_count .eq. 2) then
         fformat = 1
      else if (col_count .eq. 4) then
         fformat = 2
      else
        call raise_exception('First line of the beta decay file has an incompatible amount of columns. '//&
                             'It should be either 2 or 4, but got '//int_to_str(col_count)//".",&
                             "read_beta_decays")!! TODO Give error code!
      end if

      ! Rewind the file again to start reading
      rewind(sourcefile)

      i = 0
      read_loop: do
         if (fformat .eq. 1) then
            read(sourcefile,*, iostat = read_stat) parent,halflife
            Qnu  = -1
            Qtot = -1
         elseif (fformat .eq. 2) then
            read(sourcefile,*, iostat = read_stat) parent,halflife,Qtot,Qnu
         end if
         read(sourcefile,*, iostat = read_stat) Pn(1:11)

         if (read_stat /= 0) exit read_loop
         par_ind_tmp = benam(parent)
         if (par_ind_tmp .eq. 0) cycle read_loop

         ! Take care of nuclei that are not included.
         ! Find out decay products
         z_tmp = isotope(par_ind_tmp)%p_nr
         n_tmp = isotope(par_ind_tmp)%n_nr
         a_tmp = z_tmp + n_tmp

         do j=1,11
            ! Daughter index
            index_tmp = findaz(a_tmp-(j-1),z_tmp+1)
            if (index_tmp .ne. -1) exit
         end do

         ! There was no daughter nucleus included at all.
         ! Don't include this rate
         if (j .eq. 12) cycle read_loop

         ! Next rate
         i = i + 1

         ! Fill the struct with decay information
         beta_pn(i)%halflife   = halflife
         beta_pn(i)%idx_parent = par_ind_tmp
         beta_pn(i)%Z_parent   = z_tmp
         beta_pn(i)%N_parent   = n_tmp
         beta_pn(i)%Av_Qnu     = Qnu
         beta_pn(i)%Av_Qtot    = Qtot

         do j=1,11
            index_tmp = findaz(a_tmp-(j-1),z_tmp+1)
            beta_pn(i)%idx_daughter(j) = index_tmp
            beta_pn(i)%Pn(j)           = Pn(j)
            ! Fill the Q-value
            beta_pn(i)%Qval(j)         = 0
            if (beta_pn(i)%idx_daughter(j) .ne. -1) then
                beta_pn(i)%Qval(j)  = isotope(par_ind_tmp)%mass_exc - &
                                     (isotope(ineu)%mass_exc*float(j-1)+isotope(beta_pn(i)%idx_daughter(j))%mass_exc)
            end if
         end do
         probnot = 0
         do j=11,1,-1
            if (beta_pn(i)%idx_daughter(j) .eq. -1) then
                probnot = probnot + beta_pn(i)%Pn(j)
                beta_pn(i)%Pn(j) = 0
            end if

            if ((beta_pn(i)%idx_daughter(j) .ne. -1) .and. (probnot .ne. 0)) then
                beta_pn(i)%Pn(j) = beta_pn(i)%Pn(j) + probnot
                probnot = 0
            end if
         end do
         ! Take also care if the very left nuclei where not included
         ! Loop through it again but from left to right
         if (probnot .ne. 0) then
            do j=1,11,1
                if (beta_pn(i)%idx_daughter(j) .eq. -1) then
                    probnot = probnot + beta_pn(i)%Pn(j)
                    beta_pn(i)%Pn(j) = 0
                end if

                if ((beta_pn(i)%idx_daughter(j) .ne. -1) .and. (probnot .ne. 0)) then
                    beta_pn(i)%Pn(j) = beta_pn(i)%Pn(j) + probnot
                    probnot = 0
                end if
             end do
         end if

         if (beta_pn(i)%Av_Qtot .eq. -1) then
            do j=1,11,1
                beta_pn(i)%Av_Qtot = beta_pn(i)%Av_Qtot + beta_pn(i)%Pn(j)*beta_pn(i)%Qval(j)
            end do
         end if

         ! Lets see, this should always work. If not there is something conceptional wrong. Raise exception...
        if (probnot .ne. 0) then
            call raise_exception('No daughter nucleus of beta decay file seemed to be implemented.',&
                                 "read_beta_decays")!! TODO Give error code!
        end if
      end do read_loop

      INFO_EXIT("read_beta_decays")
   end subroutine read_beta_decays






end module beta_decay_rate_module
