!> @file alpha_decay_rate_module.f90
!!
!! The error file code for this file is ***W46***.
!! Contains the module \ref alpha_decay_rate_module.


!>
!! @brief This module contains subroutines to add alpha decays
!!
!! The alpha decays given by an external file. The default file in WinNet uses a file that has been
!! calculated using the Viola-Seaborg formula. Hereby, a parameterization of the alpha-decay half life
!! of [Dong & Ren 2005](https://ui.adsabs.harvard.edu/abs/2005EPJA...26...69D/abstract) was used.
!! Via the Viola-Seaborg formula it is possible to calculate alpha-decay half lifes by knowing only the Qvalue of the
!! alpha decay. Within this model, reaclib rates can get replaced or additional alpha-decay rates can be added.
!!
!! @author Moritz Reichert
!! @date 10.03.23
#include "macros.h"
module alpha_decay_rate_module
   use global_class, only: reactionrate_type
   implicit none

   type(reactionrate_type),dimension(:),allocatable,private   :: alpha_dec_rate       !< Array storing the reaction rates
   integer,private                                            :: nalpha_dec           !< Number of alpha decays
   logical,private                                            :: include_alpha_decays !< Flag whether alpha decays should be included or not

   ! Helper variables for ignoring specific sources
   character(len=4),allocatable,dimension(:),private :: src_ignore
   integer,private                                   :: src_ignore_length

   character(15), private, parameter :: fileformat="(A5,3X,1pE11.5)"

   !
   ! Public and private fields and methods of the module
   !
   public:: &
     init_alpha_decay_rates, merge_alpha_decays

   private:: &
     count_reactions, unify_rate_array, reaction_string, read_reactions, create_verbose_output_file


contains


   !> Initialize alpha decay rates.
   !!
   !! This subroutine counts and reads the rates into the array
   !! \ref beta_decays.
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine init_alpha_decay_rates()
        use parameter_class,     only: alpha_decay_src_ignore, use_alpha_decay_file
        use nucstuff_class,      only: analyze_src_string
        use file_handling_class
        implicit none
        integer     :: beta_unit  !< File ID of external beta decay file
        integer     :: alloc_stat !< Allocation status


        if (use_alpha_decay_file) then
            ! Alpha decays will be used
            include_alpha_decays = .True.

            ! Check if some reactions should not get replaced
            call analyze_src_string(alpha_decay_src_ignore,src_ignore,src_ignore_length)

            ! Count the amount of possible alpha decays to be added
            call count_reactions(nalpha_dec)

            ! Create an array with alpha decays
            call read_reactions(alpha_dec_rate,nalpha_dec)

            ! Say how many there were
            if (VERBOSE_LEVEL .ge. 1) then
                call write_data_to_std_out("Amount alpha-decay format rates",int_to_str(nalpha_dec))
            end if

        end if

   end subroutine init_alpha_decay_rates


   !> Creates a file with the decays
   !!
   !! This file is created for verbose_levels greater than 2. An example
   !! could look like
   !! \file{
   !!  Alpha decay rates that additionally have been added
   !! po181 =>   he4 + pb177  T_1/2: 1.49E-09 s;  a0: 2.00E+01  Qa: 1.02E+01
   !! po182 =>   he4 + pb178  T_1/2: 3.89E-09 s;  a0: 1.90E+01  Qa: 9.56E+00
   !! ..
   !! }
   !! @author M. Reichert
   !! @date 10.03.23
   subroutine create_verbose_output_file(alpha_decay_rates,nalpha)
    use file_handling_class
    implicit none
    type(reactionrate_type),dimension(:),intent(in) :: alpha_decay_rates !< Array storing the reaction rates
    integer,intent(in)                              :: nalpha            !< Number of alpha decays
    ! Internal variables
    integer :: i !< Loop variable
    integer :: alpha_unit !< File ID of external alpha decay file

    ! Open the file
    alpha_unit= open_outfile("debug_alpha_decay_rates.txt")

    ! Write the header
    write(alpha_unit,*) "Alpha decay rates that additionally have been added"
    write(alpha_unit,*)
    ! Write the header
    do i=1,nalpha
        write(alpha_unit,"(A)") reaction_string(alpha_decay_rates(i))
    end do

   end subroutine create_verbose_output_file


   !>
   !! @brief Return a string to represent a given reaction
   !!
   !! This routine is useful for error messages in case
   !! a rate is not working as intented and is used for
   !! verbose output only.
   !!
   !! ### Example
   !! For rrate(i) being the alpha-decay of po181
   !! ~~~~~~~~~~~.f90
   !! a = reaction_string(rrate(i))
   !! ~~~~~~~~~~~
   !! will return a="po181 =>   he4 + pb177  T_1/2: 1.49E-09 s;  a0: 2.00E+01  Qa: 1.02E+01".
   !!
   !! @author Moritz Reichert
   function reaction_string(reac)
    use global_class,    only: net_names, reactionrate_type
    use error_msg_class, only: num_to_str
    use benam_class,     only: get_net_name
    implicit none
    type(reactionrate_type),intent(in) :: reac
    character(150) :: reaction_string
    integer       :: i
    real(r_kind)  :: a0
    logical       :: s_prod

    reaction_string = trim(adjustl(net_names(reac%parts(1))))

    do i =2, 6
       s_prod = .false.
       ! Make a cut between educts and products
       select case(reac%group)
          case(1:3,11)
             if (i .eq. 2) then
                reaction_string = trim(adjustl(reaction_string))//" =>"
                s_prod = .true.
             end if
          case(4:7)
             if (i .eq. 3) then
                reaction_string = trim(adjustl(reaction_string))//" =>"
                s_prod = .true.
             end if
          case(8:9)
             if (i .eq. 4) then
                reaction_string = trim(adjustl(reaction_string))//" =>"
                s_prod = .true.
             end if
          case(10)
             if (i .eq. 5) then
                reaction_string = trim(adjustl(reaction_string))//" =>"
                s_prod = .true.
             end if
       end select

       ! Write the names of the isotopes
       if (reac%parts(i) .ne. 0) then
          if (.not. s_prod) then
             reaction_string = trim(adjustl(reaction_string))//" + "//&
                               get_net_name(reac%parts(i))
          else
             reaction_string = trim(adjustl(reaction_string))//" "//&
                               get_net_name(reac%parts(i))
          end if
       end if

    end do
    ! Add half life and reaclib parameter
    a0 =  reac%param(1)
    reaction_string = trim(adjustl(reaction_string))//"  T_1/2: "//num_to_str(dLOG(2.0d0)/dexp(a0))//" s;  a0: "//num_to_str(a0)//"  Qa: "//num_to_str(reac%q_value)

    return
 end function reaction_string


   !> Merge alpha decays into the larger rate array.
   !!
   !! The return value of this routine will be a larger rate array
   !! and a new length
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine merge_alpha_decays(rrate_array,rrate_length)
    use error_msg_class,  only: raise_exception
    use mergesort_module, only: rrate_ms,rrate_sort
    implicit none
    type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate_array  !< Large rate array, containing all reactions
    integer,intent(inout)                                          :: rrate_length !< length of rrate_array
    integer                                                        :: alloc_stat   !< Allocation state
    type(reactionrate_type),dimension(:),allocatable               :: helper       !< Large rate array, containing all reactions
    integer                                                        :: new_length   !< new length of the array

    if (include_alpha_decays) then
        ! Only merge if there are rates to be added
        if (nalpha_dec .ne. 0) then
            ! No reaclib rates are present, create the rate array only with alpha decays
            if (.not. allocated(rrate_array)) then

               rrate_length = nalpha_dec
               !-- Allocate the reaclib rate array
               allocate(rrate_array(nalpha_dec),stat=alloc_stat)
               if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                          "merge_alpha_decays",&
                                                          460001)
               rrate_array(1:nalpha_dec) = alpha_dec_rate(1:nalpha_dec)
            else
                ! Merge the rates as you desire
                call unify_rate_array(rrate_array,alpha_dec_rate,helper,rrate_length,nalpha_dec,new_length)
                ! Deallocate the array and allocate with new size
                deallocate(rrate_array,stat=alloc_stat)
                if ( alloc_stat /= 0) call raise_exception('Deallocation of "rrate_array" failed.',&
                                                           "merge_alpha_decays",&
                                                           460002)

                allocate(rrate_array(new_length),stat=alloc_stat)
                if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                           "merge_alpha_decays",&
                                                           460001)

                rrate_array(1:new_length) = helper(1:new_length)
                rrate_length = new_length

                ! Now helper array can be deallocated
                deallocate(helper,stat=alloc_stat)
                if ( alloc_stat /= 0) call raise_exception('Deallocation of "helper" failed.',&
                                                           "merge_alpha_decays",&
                                                           460002)
            end if

            ! Deallocate the alpha decay array to clean memory
            deallocate(alpha_dec_rate,stat=alloc_stat)
            if ( alloc_stat /= 0) call raise_exception('Deallocation of "alpha_dec_rate" failed.',&
                                                       "merge_alpha_decays",&
                                                       460002)
        end if

    end if

   end subroutine merge_alpha_decays


   !> Merge reaclib rates and alpha decays into one array.
   !!
   !! Hereby the rates can be taken as addition or also replace the reaclib rates.
   !!
   !! @author M. Reichert
   !! @date 25.01.21
   subroutine unify_rate_array(rrate_array,alpha_dec_rate_array,merged_array,nrrate,nalpha,ntot)
      use global_class,    only: ihe4,isotope
      use benam_class,     only: findaz
      use parameter_class, only: alpha_decay_ignore_all
      use error_msg_class, only: write_data_to_std_out,raise_exception,int_to_str
      implicit none
      ! Declare the input
      type(reactionrate_type),dimension(nrrate),intent(inout)      :: rrate_array          !< Reaclib rates
      type(reactionrate_type),dimension(nalpha),intent(inout)      :: alpha_dec_rate_array !< Alpha decay rates
      type(reactionrate_type),dimension(:),allocatable,intent(out) :: merged_array         !< Combined, merged array
      integer,intent(in)                                           :: nrrate               !< Length of rrate_array
      integer,intent(in)                                           :: nalpha               !< Length of alpha_dec_rate_array
      integer,intent(out)                                          :: ntot                 !< Length of merged_array
      ! Internal variables
      integer                    :: i, j            !< Loop variable
      logical, dimension(nrrate) :: mask_rrate      !< Mask for rrate_array
      logical, dimension(nalpha) :: mask_alpha      !< Mask for alpha_dec_rate_array
      integer                    :: replace_reaclib !< Number of replaced reaclib rates
      integer                    :: replace_alpha   !< Number of removed alpha decays
      integer                    :: alloc_stat      !< Allocation status variable
      integer                    :: ztmp            !< proton number
      integer                    :: atmp            !< mass number
      integer                    :: daughter        !< index of daughter nucleus

      ! Create masks for both arrays
      mask_rrate(:) = .True.
      mask_alpha(:) = .True.

      ! Keep track how many rates got replaced
      replace_reaclib = 0
      replace_alpha   = 0

      do i=1,nrrate
          ! Alpha decays are chapter 2, weak, and have an alpha in the products
          if (rrate_array(i)%group .ne. 2) cycle
          if (.not. rrate_array(i)%is_weak) cycle
          if (.not. ((rrate_array(i)%parts(2) .eq. ihe4) .or. (rrate_array(i)%parts(3) .eq. ihe4))) cycle

          ! Also check if it is not an beta delayed alpha-decay
          if (rrate_array(i)%parts(2) .eq. ihe4) then
              ! Daughter is the second particle
              daughter = rrate_array(i)%parts(3)
          else
              ! Daughter is the first particle
              daughter = rrate_array(i)%parts(2)
          end if

          if (isotope(daughter)%p_nr .ne. isotope(rrate_array(i)%parts(1))%p_nr-2) then
              ! It is a beta delayed alpha decay
              cycle
          end if

          ! Check which rate has to be replaced for this reaction
          if ((any(src_ignore .eq. adjustr(rrate_array(i)%source))) .or. (alpha_decay_ignore_all)) then
              inner_loop: do j = 1,nalpha
                  if (rrate_array(i)%parts(1) .eq. alpha_dec_rate_array(j)%parts(1)) then
                      ! The same reaction is present in both arrays, remove the alpha rate here
                      mask_alpha(j) = .False.
                      replace_alpha = replace_alpha + 1
                      exit inner_loop
                  end if
              end do inner_loop
          else
              ! Reaclib alpha decay rate is replaced by the alpha decay rate here
              mask_rrate(i) = .False.
              replace_reaclib = replace_reaclib + 1
          end if
      end do

      ! Total amount of reactions
      ntot = nrrate + nalpha - replace_alpha - replace_reaclib
      allocate(merged_array(ntot),stat=alloc_stat)
      if (alloc_stat /= 0) call raise_exception("Allocation of 'merged_array' failed.",&
                                                "unify_rate_array",&
                                                460001)

      ! Write things to the merged array
      merged_array(1:nrrate-replace_reaclib)      = pack(rrate_array(1:nrrate),mask_rrate)
      merged_array(nrrate-replace_reaclib+1:ntot) = pack(alpha_dec_rate_array(1:nalpha),mask_alpha)

      ! Say something
      if (VERBOSE_LEVEL .ge. 2) then
          call write_data_to_std_out("Amount alpha decays replaced in Reaclib",int_to_str(replace_reaclib))
          call write_data_to_std_out("Amount alpha decays replaced in V.S. eq.",int_to_str(replace_alpha))
      end if
      ! Write a file with the rates used
      if (VERBOSE_LEVEL .ge. 3) then
          call create_verbose_output_file(pack(alpha_dec_rate_array(1:nalpha),mask_alpha),nalpha - replace_alpha)
      end if

   end subroutine unify_rate_array


   !> Count the amount of possible alpha decays
   !!
   !! Check which isotopes between \ref alpha_decay_zmin and \ref alpha_decay_zmax should
   !! have an alpha decay.
   !!
   !! @author M. Reichert
   !! @date 10.03.23
   subroutine count_reactions(n)
    use global_class,    only: isotope
    use parameter_class, only: alpha_decay_zmin, alpha_decay_zmax,&
                               alpha_decay_file
    use benam_class,     only: findaz, benam
    use file_handling_class
    implicit none
    integer,intent(out) :: n !< Number of possible alpha decays
    integer             :: i             !< Loop variable
    integer             :: idx_daughter  !< Index of daughter
    integer             :: idx_parent    !< Index of parent
    integer             :: N_tmp         !< Number of neutrons
    integer             :: Z_tmp         !< Number of protons
    integer             :: file_id       !< File id
    logical             :: include_rate  !< Flag if the decay should be included
    character(5)        :: nametmp       !< Name of the isotope
    real(r_kind)        :: thalf         !< Half life of the isotope
    integer             :: read_stat     !< Read status variable

    ! Reset the counter
    n = 0
    ! Open the file
    file_id = open_infile(alpha_decay_file)

    ! Read the file
    loop: do
        read(file_id,fileformat, iostat = read_stat) nametmp, thalf
        if (read_stat /= 0) exit loop

        ! Skip nuclei that are not included
        idx_parent = benam(nametmp)
        if (idx_parent .eq. 0) cycle loop
        ! Skip nuclei that are not in desired Z range
        if ((isotope(idx_parent)%p_nr .lt. alpha_decay_zmin) .and. (isotope(idx_parent)%p_nr .gt. alpha_decay_zmax)) cycle loop
        idx_daughter = findaz(isotope(idx_parent)%mass-4,isotope(idx_parent)%p_nr-2)
        ! Skip nuclei if the daughter is not included
        if (idx_daughter .le. 0) cycle loop
        ! Increase amount of reactions
        n = n+1

    end do loop

   ! Close the fileagain
    close(file_id)

   end subroutine count_reactions



   !> Read the alpha decays
   !!
   !! All isotopes between \ref alpha_decay_zmin and \ref alpha_decay_zmax should
   !! have an alpha decay. This subroutine reads the file \ref alpha_decay_file
   !! and creates the alpha decay rate array.
   !! An example of the file could look like:
   !! \file{
   !! sb103   7.84058e+10
   !! sb104   4.29655e+05
   !! sb105   1.37511e+10
   !! sb106   1.03780e+17
   !! te103   8.30806e+02
   !! ...
   !! }
   !!
   !! @author M. Reichert
   !! @date 14.03.23
   subroutine read_reactions(alpha_rate_array,n)
    use global_class,    only: isotope,net_size,ihe4
    use error_msg_class, only: raise_exception
    use benam_class,     only: findaz,getcoefficients,benam
    use parameter_class, only: alpha_decay_zmin, alpha_decay_zmax,&
                               alpha_decay_file
    use file_handling_class
    implicit none
    type(reactionrate_type),dimension(:),allocatable, intent(out) :: alpha_rate_array  ! Alpha decay rate array
    integer,intent(in)        :: n            !< Length of the rates
    integer                   :: counter      !< Counter
    integer                   :: idx_daughter !< Index of daughter nucleus
    integer                   :: N_tmp        !< Number of neutrons
    integer                   :: Z_tmp        !< Number of protons
    real(r_kind)              :: mexc_alpha   !< Mass excess of alpha particle
    real(r_kind)              :: mexc_daughter!< Mass excess of daughter
    real(r_kind)              :: mexc_parent  !< Mass excess of parent
    real(r_kind)              :: Qa           !< Q-value of alpha decay [MeV]
    real(r_kind)              :: a0_tmp       !< Reaclib a0 parameter
    real(r_kind)              :: thalf        !< Half life of the isotope
    integer                   :: file_id      !< File id
    integer                   :: read_stat    !< Read status variable
    integer                   :: alloc_stat   !< Allocation status variable
    character(5)              :: nametmp      !< Name of the isotope
    integer                   :: idx_parent   !< Index of parent
    logical                   :: include_rate !< Flag if the decay should be included or was problematic
    type(reactionrate_type)   :: rrate_tmp    !< Temporary reaction rate


    ! Allocate the array
    allocate(alpha_rate_array(n),stat=alloc_stat)
    if ( alloc_stat /= 0)  call raise_exception('Allocation of "alpha_rate_array" failed',&
                                                "read_reactions",&
                                                460001)
    ! Initialize the counter
    counter = 1
    ! Get mass excess of alpha particle
    mexc_alpha = isotope(ihe4)%mass_exc

   ! Open the file
    file_id = open_infile(alpha_decay_file)

    ! Read the file
    loop: do
        read(file_id,fileformat, iostat = read_stat) nametmp, thalf
        if (read_stat /= 0) exit loop

        idx_parent = benam(nametmp)
        ! Skip nuclei that are not included
        if (idx_parent .eq. 0) cycle loop
        ! Skip nuclei that are not in desired Z range
        if ((isotope(idx_parent)%p_nr .lt. alpha_decay_zmin) .and. (isotope(idx_parent)%p_nr .gt. alpha_decay_zmax)) cycle loop
        idx_daughter = findaz(isotope(idx_parent)%mass-4,isotope(idx_parent)%p_nr-2)
        ! Skip nuclei if the daughter is not included
        if (idx_daughter .le. 0) cycle loop

        ! Calculate the Q-value
        mexc_daughter = isotope(idx_daughter)%mass_exc
        mexc_parent   = isotope(idx_parent)%mass_exc
        Qa = mexc_parent-mexc_daughter-mexc_alpha
        ! Also save the numbers of neutrons and protons
        N_tmp = isotope(idx_parent)%n_nr
        Z_tmp = isotope(idx_parent)%p_nr

        ! Calculate a0 parameter
        a0_tmp = dlog((dLOG(2.0d0)/(thalf)))

        ! Get the rate
        rrate_tmp%parts(:)    = 0            ! Initialize first with zero
        rrate_tmp%source      = "aext"       ! Added alpha-decay
        rrate_tmp%is_reverse  = .false.
        rrate_tmp%cached      = -1
        rrate_tmp%is_resonant = .false.
        rrate_tmp%is_weak     = .true.
        rrate_tmp%is_const    = .true.
        rrate_tmp%q_value     = Qa           ! Set the Qvalue
        rrate_tmp%reac_src    = rrs_aext     ! Remember that this is an added alpha-decay
        rrate_tmp%reac_type   = rrt_alpd     ! It's an alpha decay
        rrate_tmp%one_over_n_fac = 1.0d0
        rrate_tmp%group       =  2           ! Alpha decays are in chapter 2
        rrate_tmp%parts(1)    = idx_parent   ! parent nucleus index
        rrate_tmp%parts(2)    = ihe4         ! alpha index
        rrate_tmp%parts(3)    = idx_daughter ! daughter nucleus index
        rrate_tmp%param(:)    = 0            ! Reset parameter
        rrate_tmp%param(1)    = a0_tmp       ! a0 parameter as calculated
        rrate_tmp%nu_frac     = 0            ! Fraction radiated away by neutrinos
        ! Save the rate and increase the counter for the next rate
        alpha_rate_array(counter) = rrate_tmp
        counter = counter + 1
    end do loop

    ! Make coefficients right
    call getcoefficients(alpha_rate_array,n)

   end subroutine read_reactions


end module alpha_decay_rate_module