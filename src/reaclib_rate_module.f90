!> @file reaclib_rate_module.f90
!!
!! The error file code for this file is ***W38***.
!!
!! Contains the module \ref reaclib_rate_module


!> @brief Module that deals with reaclib reaction rates.
!!
!! This module initializes reaclib reaction rates and stores it into a
!! rate array.
#include "macros.h"
module reaclib_rate_module
   use global_class, only: reactionrate_type
   implicit none

  integer,private     :: nrea                                      !< Amount of reaclib rates
  integer,private     :: n_w,n_ng,n_gn,n_ap,n_pa,n_o,n_pg,n_gp     !< Amount of individual reaction types
  integer,private     :: n_np,n_pn,n_ag,n_ga,n_an,n_na             !< Amount of individual reaction types
  integer,private     :: n_bm,n_bp,n_ad,n_pe,n_ne,n_ec             !< Amount of individual reaction types
  type(reactionrate_type), dimension(:), allocatable, private :: rrate_reaclib !< Reaclib reaction rates

  real(r_kind),private :: infty
  !
  ! Public and private fields and methods of the module
  !
  public:: &
      init_reaclib_rates, merge_reaclib_rates, set_reaction_type
  private:: &
      read_reaclib, count_reaclib_rates, write_reac_verbose_out

contains


!> Count and read reaclib reactions
!!
!! This subroutine will fill \ref nrea with the amount
!! of reaclib reactions and the \ref rrate_reaclib array
!! with reaction rates from a reaclib file.
!!
!! @author M. Reichert
!! @date 25.01.21
subroutine init_reaclib_rates()
   use error_msg_class, only: raise_exception
   implicit none
   integer    :: alloc_stat  !< Status of allocation

   ! Amount of individual reaction types
   n_ne=0; n_ng=0; n_gn=0; n_ap=0; n_pa=0; n_o=0
   n_np=0; n_pn=0; n_ag=0; n_ga=0; n_an=0; n_na=0
   n_pg=0; n_gp=0; n_bm=0; n_bp=0; n_ad=0; n_pe=0
   n_ec=0

   ! Get a variable for comparison for overflows
   infty = huge(infty)

   !-- Count the amount of rates
   call count_reaclib_rates()

   !-- Allocate the reaclib rate array
   allocate(rrate_reaclib(nrea),stat=alloc_stat)
   if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_reaclib" failed.',&
                                              "init_reaclib_rates",380001)

   !-- Read the reaclib rates into rrate_reaclib
   call read_reaclib()

   call write_reac_verbose_out()

end subroutine init_reaclib_rates



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
      call write_data_to_std_out("Amount reaclib rates",int_to_str(nrea))
   end if
   if (VERBOSE_LEVEL .ge. 2) then
      if (nrea .gt. 0) write(*,"(A)") ""
      if (nrea .gt. 0) write(*,"(A)") "    Reaclib rates:  "
      if (nrea .gt. 0) write(*,"(A)") "   |----------------|"
      tmp = int_to_str(nrea)
      if (nrea .gt. 0) write(*,"(A)") "   | Total :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_ng)
      if (n_ng .gt. 0) write(*,"(A)") "   | (n,g) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_gn)
      if (n_gn .gt. 0) write(*,"(A)") "   | (g,n) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_pg)
      if (n_pg .gt. 0) write(*,"(A)") "   | (p,g) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_gp)
      if (n_gp .gt. 0) write(*,"(A)") "   | (g,p) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_ag)
      if (n_ag .gt. 0) write(*,"(A)") "   | (a,g) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_ga)
      if (n_ga .gt. 0) write(*,"(A)") "   | (g,a) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_an)
      if (n_an .gt. 0) write(*,"(A)") "   | (a,n) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_na)
      if (n_na .gt. 0) write(*,"(A)") "   | (n,a) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_ap)
      if (n_an .gt. 0) write(*,"(A)") "   | (a,p) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_pa)
      if (n_na .gt. 0) write(*,"(A)") "   | (p,a) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_pn)
      if (n_pn .gt. 0) write(*,"(A)") "   | (p,n) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_np)
      if (n_np .gt. 0) write(*,"(A)") "   | (n,p) :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_bm)
      if (n_bm .gt. 0) write(*,"(A)") "   | beta- :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_bp)
      if (n_bp .gt. 0) write(*,"(A)") "   | beta+ :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_ec)
      if (n_ec .gt. 0) write(*,"(A)") "   | EC    :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_ad)
      if (n_ad .gt. 0) write(*,"(A)") "   | a-dec.:"//adjustr(tmp)//" |"
      tmp = int_to_str(n_ne)
      if (n_ne .gt. 0) write(*,"(A)") "   | n emis:"//adjustr(tmp)//" |"
      tmp = int_to_str(n_pe)
      if (n_pe .gt. 0) write(*,"(A)") "   | p emis:"//adjustr(tmp)//" |"
      tmp = int_to_str(n_o)
      if (n_o  .gt. 0) write(*,"(A)") "   | other :"//adjustr(tmp)//" |"
      if (nrea .gt. 0) write(*,"(A)") "   |----------------|"
      if (nrea .gt. 0) write(*,"(A)") ""
   end if

end subroutine write_reac_verbose_out





!> Calculate a reaclib rate.
!!
!! The rate is calculated according to the reaclib fit function:
!! \f[ \lambda = \exp\left[ a_0 + \sum_{i=1}^5 a_i T_9^{(2i-5)/3} +a_6 \log T_9 \right]  \f]
!!
!! @see nucstuff_class::calc_t9_pow
!!
!! \b Edited:
!!      - 26.07.22, MR: created this subroutine from code parts in jacobian_class.
!! .
subroutine calculate_reacl_rate(rrate,rat_calc)
  use global_class,   only: reactionrate_type
  use nucstuff_class, only: t9_pow
  use benam_class,    only: reaction_string
  implicit none
  ! Declare the pass
  type(reactionrate_type),intent(in)  :: rrate    !< rate instance
  real(r_kind),intent(out)            :: rat_calc !< rate value
  ! Internal variables
  integer :: j !< Loop variable

  rat_calc = 0.d0
  do j=1,9
     rat_calc = rat_calc +rrate%param(j)*t9_pow(j)
  end do

  if (rat_calc .gt. dlog(infty)) then
    ! Complain about overflow
    if (VERBOSE_LEVEL .ge. 2) then
        print*,"Warning, rate overflow! Rate was "
        print*,reaction_string(rrate)
    end if
    rat_calc = dlog(infty)
  end if

  rat_calc = dexp(rat_calc)

end subroutine calculate_reacl_rate




!> Merge the reaclib rates into a larger array
!!
!! The subroutine combines a larger rate array with
!! \ref rrate_reaclib. Afterwards it deallocates \ref rrate_reaclib.
!!
!! @returns Large rate array with new size, containing the reaclib reactions
!!
!! @author M. Reichert
!! @date 25.01.21
subroutine merge_reaclib_rates(rrate_array,rrate_length)
   use error_msg_class, only: raise_exception
   implicit none
   type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate_array  !< Large rate array, containing all reactions
   integer,intent(inout)                                          :: rrate_length !< length of rrate_array
   type(reactionrate_type),dimension(:),allocatable               :: rrate_tmp    !< Temporary rate array
   integer                                                        :: alloc_stat   !< Allocation state
   integer                                                        :: new_length   !< New length of rrate_array

   new_length = rrate_length+nrea
   if (nrea .ne. 0) then
     if (.not. allocated(rrate_array)) then
        !-- Allocate the reaclib rate array
        allocate(rrate_array(nrea),stat=alloc_stat)
        if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                   "merge_reaclib_rates",380001)
        rrate_array(1:nrea) = rrate_reaclib(1:nrea)
     else
        !-- Allocate a temporary array to store the content of rrate_array
        allocate(rrate_tmp(rrate_length),stat=alloc_stat)
        if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_tmp" failed.',&
                                                   "merge_reaclib_rates",380001)
        rrate_tmp(1:rrate_length) = rrate_array(1:rrate_length)

        !-- Deallocate the input array
        deallocate(rrate_array)
        allocate(rrate_array(new_length),stat=alloc_stat)
        if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                   "merge_reaclib_rates",380001)
        rrate_array(1:rrate_length)             = rrate_tmp(1:rrate_length)
        rrate_array(rrate_length+1:new_length)  = rrate_reaclib(1:nrea)

        deallocate(rrate_tmp)
     end if
     !-- Deallocate the reaclib rate array
     deallocate(rrate_reaclib)
   end if
   !-- Output the new length
   rrate_length = new_length

end subroutine merge_reaclib_rates




!> Read reaclib reaction rates
!!
!! This routine reads the reaclib file and stores
!! the reactions into \ref rrate_reaclib. An example file
!! looks like:
!! \file{...
!!       li9  be9                            wc12w     1.36060e+01
!! 6.501820e-01 0.000000e+00 0.000000e+00 0.000000e+00
!! 0.000000e+00 0.000000e+00 0.000000e+00
!!      li11 be11                            wc12w     2.05510e+01
!! 2.322150e+00 0.000000e+00 0.000000e+00 0.000000e+00
!! 0.000000e+00 0.000000e+00 0.000000e+00
!! ...}
!!
!! @see [Reaclib file format](https://reaclib.jinaweb.org/help.php?topic=reaclib_format&intCurrentNum=0)
!!
!! \b Edited:
!!     - 25.01.21, MR - Moved it from network_init to this subroutine
!!     - 28.07.22, MR - Introduced new chapters
!! .
subroutine read_reaclib()
   use parameter_class, only: reaclib_file, heating_frac,&
                              use_neutrino_loss_file
   use benam_class,     only: benam,getcoefficients
   use nucstuff_class,  only: get_nr_reactants
   use global_class,    only: Qnuloss,ipro,ineu
   use file_handling_class
   use format_class
   implicit none
   integer                    :: i           !< Count of the rates
   integer                    :: j           !< Loop variable
   integer                    :: read_stat   !< Read status
   integer                    :: reaclib     !< File id of reaclib
   character(5), dimension(6) :: parts       !< Participating nuclei names
   integer, dimension(6)      :: parts_index !< Participating nuclei indices
   real(r_kind), dimension(9) :: params      !< Reaclib fitting parameter
   character(4)               :: src         !< Reaclib source label
   character(1)               :: res         !< Reaclib weak flag label
   character(1)               :: rev         !< Reaclib reverse label
   real(r_kind)               :: q           !< Reaclib Q-Value
   integer                    :: grp         !< Reaclib chapter
   integer                    :: group_index !< Reaclib chapter storage

   ! Open the reaclib
   reaclib= open_infile(reaclib_file)

   ! Set flags false by default
   rrate_reaclib%is_weak     = .false.
   rrate_reaclib%is_resonant = .false.
   rrate_reaclib%is_reverse  = .false.
   rrate_reaclib%cached      = -1
   rrate_reaclib%reac_type   = rrt_o     ! Initialize as "other" reaction
   rrate_reaclib%reac_src    = rrs_reacl ! Initialize as reaclib reaction
   rrate_reaclib%nu_frac     = 0         ! No neutrinos, only for beta decays later on


   i = 1
   outer_loop: do
      read(reaclib,my_format(1), iostat = read_stat)  &
           grp, parts(1:6), src, res, rev, q
      if (read_stat /= 0) exit outer_loop
      !!read(reaclib,my_format(2)) params(1:7) ! *** my_format(2) reads only 4 numbers but somehow this works...
      read(reaclib,"(4E13.6)") params(1:4)
      read(reaclib,"(5E13.6)") params(5:9)

 !--- fission rates are read in separately. For a transition period, the fission rates are probably still present in the reaclib file
      if ((src.eq.'ms99').or.(src.eq.'mp01').or.(src.eq.'sfis').or.(src.eq.'fiss')) cycle outer_loop
 !--- grp =/ 0 only for group separators in reaclib
 !    explanation: reaclib is structured to have group separators with same number of chars
 !    of blocks containing data for each reaction; variable grp gets overwritten with a 0
 !    in blocks containing data and is not 0 only for separators
      select case (grp)
      case(1:11)
         group_index = grp
         cycle
      case default
         parts_index = 0
         inner_loop: do j=1,6
            if (parts(j) .eq. '     ') exit inner_loop
            parts_index(j) = benam(parts(j))
 !----- parts_index(i) = 0 if nuclide i is not in network
            if (parts_index(j) .eq. 0) cycle outer_loop
         end do inner_loop
      end select

      rrate_reaclib(i)%group = group_index
      rrate_reaclib(i)%parts = parts_index
      rrate_reaclib(i)%source = src
      rrate_reaclib(i)%q_value = q
      ! Give errors if the reaclib introduces something new...
      if (res == "r") then
         rrate_reaclib(i)%is_resonant = .true.
      elseif (res =="w") then
        rrate_reaclib(i)%is_weak = .true.
      elseif (res =="s") then
        rrate_reaclib(i)%is_weak = .true.
      elseif (res =="p") then
        rrate_reaclib(i)%is_weak = .true.
      elseif ((res ==" ") .or. (res =="n")) then
        continue
      else
        call raise_exception('Flag within reaclib not known! (res = "'//res//'"). '//&
                             'Source flag was '//trim(adjustl(src))//'.',&
                             "read_reaclib",380003)
      end if
      ! Same for the rev flag
      if (rev == "v") then
        rrate_reaclib(i)%is_reverse  = .true.
      elseif (rev == " ") then
        continue
      else
        call raise_exception('Flag within reaclib not known! (rev = "'//rev//'"). '//&
                             'Source flag was '//trim(adjustl(src))//'.',&
                             "read_reaclib",380004)
      end if

      rrate_reaclib(i)%param = params

      !-- Set the reaction type (e.g., (Alpha,n) )
      call set_reaction_type(rrate_reaclib(i))
      i=i+1
   end do outer_loop

   ! Close the file again
   close(reaclib)

   ! Set the fraction of energy radiated away by neutrinos
   call set_heating_frac(rrate_reaclib,nrea)

   ! get the correct coefficients to prevent double counting
   call getcoefficients(rrate_reaclib,nrea)


end subroutine read_reaclib


!> Set the fraction of energy radiated away by neutrinos
!!
!! This subroutine sets the fraction of energy radiated away by neutrinos
!! for each reaction. This is done either by a user defined file,
!! or by a default value.
!!
!! @author M. Reichert
!! @date 31.03.2023
subroutine set_heating_frac(reac_rate_array,length)
    use global_class,    only: reactionrate_type, net_size, Qnuloss
    use parameter_class, only: heating_frac, use_neutrino_loss_file

    implicit none
    type(reactionrate_type),dimension(length), intent(inout)  :: reac_rate_array
    integer, intent(in)                                       :: length
    ! internal variables
    integer                          :: i,j
    real(r_kind),dimension(net_size) :: av_Q
    real(r_kind),dimension(net_size) :: heat_f
    real(r_kind),dimension(net_size) :: rate
    real(r_kind)                     :: rate_tmp

    if (use_neutrino_loss_file) then
        av_Q(:) = 0
        rate(:) = 0

        do i=1,length
            if (reac_rate_array(i)%is_weak) then
                ! Identify if it is a beta decay
                if ((reac_rate_array(i)%reac_type == rrt_betm) .or. (reac_rate_array(i)%reac_type == rrt_betp)) then
                    rate_tmp = dexp(reac_rate_array(i)%param(1))
                    ! Identify the nuclide that is being heated
                    av_Q(reac_rate_array(i)%parts(1)) = av_Q(reac_rate_array(i)%parts(1)) + reac_rate_array(i)%q_value*rate_tmp
                    rate(reac_rate_array(i)%parts(1)) = rate(reac_rate_array(i)%parts(1)) + rate_tmp
                end if
            end if
        end do

        ! Calculate average Q value
        do i=1,net_size
            if (rate(i) .gt. 0) then
                av_Q(i) = av_Q(i)/rate(i)
            else
                av_Q(i) = 0
            end if
        end do

        ! Calculate the heating fraction
        do i=1,net_size
            if ((av_Q(i) .gt. 0) .and. (Qnuloss(i) .ne. -1)) then
                heat_f(i) = Qnuloss(i)/av_Q(i)
            else
                heat_f(i) = heating_frac
            end if
        end do
    else
        heat_f(:) = heating_frac
    end if

    do i=1,length
        if (reac_rate_array(i)%is_weak) then
            ! Identify if it is a beta decay
            if ((reac_rate_array(i)%reac_type == rrt_betm) .or. (reac_rate_array(i)%reac_type == rrt_betp) .or. &
                (reac_rate_array(i)%reac_type == rrt_ec)) then
                reac_rate_array(i)%nu_frac = heat_f(reac_rate_array(i)%parts(1))
            else
                reac_rate_array(i)%nu_frac = 0
            end if
        end if
    end do


end subroutine set_heating_frac



!> Set a flag for the reaction type
!!
!! Possible reaction types are
!! <table>
!! <caption id="multi_row">Reaction types</caption>
!! <tr><th> Reaction type    <th> Meaning
!! <tr><td> rrt_betm         <td> Beta minus
!! <tr><td> rrt_betp         <td> Beta plus
!! <tr><td> rrt_ec           <td> Electron capture
!! <tr><td> rrt_alpd         <td> Alpha-decay
!! <tr><td> rrt_nemi         <td> neutron emission
!! <tr><td> rrt_pemi         <td> proton emission
!! <tr><td> rrt_ng           <td> \f$(n,\gamma)\f$
!! <tr><td> rrt_gn           <td> \f$(\gamma,n)\f$
!! <tr><td> rrt_ag           <td> \f$(\alpha,\gamma)\f$
!! <tr><td> rrt_ga           <td> \f$(\gamma,\alpha)\f$
!! <tr><td> rrt_pg           <td> \f$(p,\gamma)\f$
!! <tr><td> rrt_gp           <td> \f$(\gamma,p)\f$
!! <tr><td> rrt_na           <td> \f$(n,\alpha)\f$
!! <tr><td> rrt_an           <td> \f$(\alpha,n)\f$
!! <tr><td> rrt_np           <td> \f$(n,p)\f$
!! <tr><td> rrt_pn           <td> \f$(p,n)\f$
!! <tr><td> rrt_pa           <td> \f$(p,\alpha)\f$
!! <tr><td> rrt_ap           <td> \f$(\alpha,p)\f$
!! <tr><td> rrt_nu           <td> Neutrino
!! <tr><td> rrt_nf           <td> Neutron induced fission
!! <tr><td> rrt_bf           <td> \f$\beta\f$-delayed fission
!! <tr><td> rrt_sf           <td> Spontaneous fission
!! <tr><td> rrt_o            <td> Other reaction
!! </table>
!! Not all flags are set in this routine (for example, neutrino and fission types).
!! The rrt_* variables are defined as preprocessor integer variables in
!! macros.h.
!!
!! @author M. Reichert
!!
!! \b Edited:
!!     - 25.01.21, MR - Moved it to separate subroutine
!! .
subroutine set_reaction_type(rr_tmp)
   use benam_class,    only: get_net_name
   use nucstuff_class, only: get_nr_products, get_nr_reactants
   use global_class,   only: ineu,ipro,ihe4,isotope
   implicit none
   type(reactionrate_type),intent(inout)      :: rr_tmp      !< Reaction rate to classify
   integer                                    :: group_index !< Reaclib chapter
   integer, dimension(6)                      :: parts_index !< Participating nuclei indices
   logical                                    :: trimmed     !< Flag to return the name trimmed
   integer                                    :: Zsum_r      !< Count for amount of neutrons
   integer                                    :: Zsum_p      !< Count for amount of neutrons
   integer                                    :: Nsum_r      !< Count for amount of neutrons
   integer                                    :: Nsum_p      !< Count for amount of neutrons
   integer                                    :: neutron_c_r !< Count for amount of neutrons
   integer                                    :: neutron_c_p !< Count for amount of neutrons
   integer                                    :: proton_c_r  !< Count for amount of protons
   integer                                    :: proton_c_p  !< Count for amount of protons
   integer                                    :: alpha_c_r   !< Count for amount of alphas
   integer                                    :: alpha_c_p   !< Count for amount of alphas
   integer                                    :: i           !< Loop variable

   trimmed = .True. ! Set some value

   parts_index = rr_tmp%parts ! Shorthand variable for the part indices
   group_index = rr_tmp%group ! Shorthand variable for the chapter

    ! Classify weak reactions (Remove hardcoded flags in future)
         if (rr_tmp%is_weak) then
            neutron_c_r = 0 ; neutron_c_r = 0 ; proton_c_r  = 0
            alpha_c_r = 0
            Zsum_r=0 ; Nsum_r=0
            ! Count amount of neutrons in reactants
            do i=1,get_nr_reactants(group_index)
               if (parts_index(i) .le. 0) cycle
               if (parts_index(i) .eq. ineu) then
                  neutron_c_r = neutron_c_r +1
               else if (parts_index(i) .eq. ipro) then
                    proton_c_r = proton_c_r +1
               else if (parts_index(i) .eq. ihe4) then
                    alpha_c_r = alpha_c_r +1
               end if
               Zsum_r = Zsum_r + isotope(parts_index(i))%p_nr
               Nsum_r = Nsum_r + isotope(parts_index(i))%n_nr
            end do

            ! Same for products
            neutron_c_p = 0 ; neutron_c_p = 0 ; proton_c_p  = 0
            alpha_c_p = 0;
            Zsum_p=0 ; Nsum_p=0

            ! Count amount of neutrons in products
            do i=get_nr_reactants(group_index)+1,get_nr_products(group_index)+get_nr_reactants(group_index)+1
               if (parts_index(i) .le. 0) cycle
               if (parts_index(i) .eq. ineu) then
                  neutron_c_p = neutron_c_p +1
               else if (parts_index(i) .eq. ipro) then
                    proton_c_p = proton_c_p +1
               else if (parts_index(i) .eq. ihe4) then
                    alpha_c_p = alpha_c_p +1
               end if
               Zsum_p = Zsum_p + isotope(parts_index(i))%p_nr
               Nsum_p = Nsum_p + isotope(parts_index(i))%n_nr
            end do

            ! Electron capture
            if (rr_tmp%source .eq. '  ec') then
               rr_tmp%reac_type = rrt_ec
               n_ec = n_ec + 1
            ! Neutron emission
            elseif ((Zsum_p .eq. Zsum_r) .and. (neutron_c_p .gt. 0)) then
                rr_tmp%reac_type = rrt_nemi
                n_ne = n_ne + 1
            ! Proton emission
            elseif ((Zsum_p .eq. Zsum_r) .and. (proton_c_p .gt. 0)) then
                rr_tmp%reac_type = rrt_pemi
                n_pe = n_pe + 1
            ! Alpha decay
            elseif ((Zsum_p .eq. Zsum_r) .and. (alpha_c_p .gt. 0)) then
                rr_tmp%reac_type = rrt_alpd
                n_ad = n_ad + 1
            ! Beta minus
            elseif ((Zsum_p-1 .eq. Zsum_r) .and. (Nsum_p+1 .eq. Nsum_r)) then
                rr_tmp%reac_type = rrt_betm
                n_bm = n_bm + 1
            ! Beta plus
            elseif ((Zsum_p+1 .eq. Zsum_r) .and. (Nsum_p-1 .eq. Nsum_r)) then
                rr_tmp%reac_type = rrt_betp
                n_bp = n_bp + 1
            else
                rr_tmp%reac_type = rrt_o
                n_o = n_o + 1
            end if

    ! Classify neutron capture reactions
         else if ((group_index.eq.4) .and. &
                  ((get_net_name(parts_index(1),trimmed).eq.'n').or.&
                   (get_net_name(parts_index(2),trimmed).eq.'n'))) then
            rr_tmp%reac_type = rrt_ng
            n_ng = n_ng +1 ! Count reaction type for statistics
    ! Classify gamma-neutron
         else if (((group_index.eq.2)).and.&
                  ((get_net_name(parts_index(2),trimmed).eq.'n').or.&
                   (get_net_name(parts_index(3),trimmed).eq.'n')))then
            rr_tmp%reac_type = rrt_gn
            n_gn = n_gn +1 ! Count reaction type for statistics
    ! Classify alpha capture reactions
         else if ((group_index.eq.4).and. &
                  ((get_net_name(parts_index(1),trimmed).eq.'he4').or. &
                   (get_net_name(parts_index(2),trimmed).eq.'he4'))) then
            rr_tmp%reac_type = rrt_ag
            n_ag = n_ag +1 ! Count reaction type for statistics
    ! Classify gamma-alpha
         else if (((group_index.eq.2)).and.&
                  ((get_net_name(parts_index(2),trimmed).eq.'he4').or.&
                   (get_net_name(parts_index(3),trimmed).eq.'he4')))then
            rr_tmp%reac_type = rrt_ga
            n_ga = n_ga +1 ! Count reaction type for statistics
    ! Classify proton capture reactions
         else if ((group_index.eq.4).and.&
                  ((get_net_name(parts_index(1),trimmed).eq.'p').or. &
                   (get_net_name(parts_index(2),trimmed).eq.'p'))) then
            rr_tmp%reac_type = rrt_pg
            n_pg = n_pg +1 ! Count reaction type for statistics
    ! Classify gamma-p
         else if (((group_index.eq.2)).and.&
                  ((get_net_name(parts_index(2),trimmed).eq.'p').or.&
                   (get_net_name(parts_index(3),trimmed).eq.'p'))) then
            rr_tmp%reac_type = rrt_gp
            n_gp = n_gp +1 ! Count reaction type for statistics
    ! Classify n-a
         else if ( (group_index.eq.5).and.&
                   ((get_net_name(parts_index(1),trimmed).eq.'n')  .or.&
                    (get_net_name(parts_index(2),trimmed).eq.'n')).and.&
                   ((get_net_name(parts_index(3),trimmed).eq.'he4').or.&
                    (get_net_name(parts_index(4),trimmed).eq.'he4'))) then
            rr_tmp%reac_type = rrt_na
            n_na = n_na +1 ! Count reaction type for statistics
    ! Classify a-n
         else if ( (group_index.eq.5).and.&
                   ((get_net_name(parts_index(1),trimmed).eq.'he4').or.&
                    (get_net_name(parts_index(2),trimmed).eq.'he4')).and.&
                   ((get_net_name(parts_index(3),trimmed).eq.'n')  .or.&
                    (get_net_name(parts_index(4),trimmed).eq.'n'))) then
            rr_tmp%reac_type = rrt_an
            n_an = n_an +1 ! Count reaction type for statistics
    ! Classify n-p
         else if ( (group_index.eq.5).and.&
                   ((get_net_name(parts_index(1),trimmed).eq.'n').or.&
                    (get_net_name(parts_index(2),trimmed).eq.'n')).and.&
                   ((get_net_name(parts_index(3),trimmed).eq.'p').or.&
                    (get_net_name(parts_index(4),trimmed).eq.'p'))) then
            rr_tmp%reac_type = rrt_np
            n_np = n_np +1 ! Count reaction type for statistics
    ! Classify p-n
         else if ( (group_index.eq.5).and.&
                   ((get_net_name(parts_index(1),trimmed).eq.'p').or.&
                    (get_net_name(parts_index(2),trimmed).eq.'p')).and.&
                   ((get_net_name(parts_index(3),trimmed).eq.'n').or.&
                    (get_net_name(parts_index(4),trimmed).eq.'n'))) then
            rr_tmp%reac_type = rrt_pn
            n_pn = n_pn +1 ! Count reaction type for statistics
    ! Classify p-a
         else if ( (group_index.eq.5).and.&
                   ((get_net_name(parts_index(1),trimmed).eq.'p')  .or.&
                    (get_net_name(parts_index(2),trimmed).eq.'p')).and.&
                   ((get_net_name(parts_index(3),trimmed).eq.'he4').or.&
                    (get_net_name(parts_index(4),trimmed).eq.'he4'))) then
           rr_tmp%reac_type = rrt_pa
           n_pa = n_pa +1 ! Count reaction type for statistics
    ! Classify a-p
         else if ( (group_index.eq.5).and.&
                   ((get_net_name(parts_index(1),trimmed).eq.'he4').or.&
                    (get_net_name(parts_index(2),trimmed).eq.'he4')).and.&
                   ((get_net_name(parts_index(3),trimmed).eq.'p')  .or.&
                    (get_net_name(parts_index(4),trimmed).eq.'p'))) then
           rr_tmp%reac_type = rrt_ap
           n_ap = n_ap +1 ! Count reaction type for statistics
         else
           rr_tmp%reac_type = rrt_o
           n_o = n_o +1 ! Count reaction type for statistics
         end if

end subroutine set_reaction_type





!> Count the amount of reaclib reactions
!!
!! After the subroutine has been called
!! \ref nrea will store the amount of reaclib reaction rates.
!!
!! \b Edited:
!!    - 25.01.21, MR - Moved this routine from network_init_module to this new subroutine
!! .
subroutine count_reaclib_rates()
   use file_handling_class
   use parameter_class, only: reaclib_file
   use benam_class,     only: benam
   use format_class
   implicit none
   integer                    :: i           !< Count of the rates
   integer                    :: j           !< Loop variable
   integer                    :: read_stat   !< Read status
   integer                    :: reaclib     !< File id of reaclib
   character(5), dimension(6) :: parts       !< Participating nuclei names
   integer, dimension(6)      :: parts_index !< Participating nuclei indices
   real(r_kind), dimension(9) :: params      !< Reaclib fitting parameter
   character(4)               :: src         !< Reaclib source label
   character(1)               :: res         !< Reaclib weak flag label
   character(1)               :: rev         !< Reaclib reverse label
   real(r_kind)               :: q           !< Reaclib Q-Value
   integer                    :: grp         !< Reaclib chapter

   i=0
   reaclib= open_infile(reaclib_file)

   outer: do
      read(reaclib,my_format(1), iostat = read_stat)  &
           grp, parts(1:6), src, res, rev, q

      if (read_stat /= 0) exit
      read(reaclib,"(4E13.6)") params(1:4)
      read(reaclib,"(5E13.6)") params(5:9)
 !--- fission rates are read in separately. For a transition period, the fission rates are probably still present in the reaclib file
      if ((src.eq.'ms99').or.(src.eq.'mp01').or.(src.eq.'sfis').or.(src.eq.'fiss')) cycle outer
 !--- grp =/ 0 only for group separators in reaclib
 !    explaination: reaclib is structured to have group separators with same number of chars
 !    of blocks containing data for each reaction; variable grp gets overwritten with a 0
 !    in blocks containing data and is not 0 only for separators
      if (grp.ne.0) cycle outer
      parts_index = 0
      inner: do j=1,6
         if (parts(j) .eq. '     ') exit inner
         parts_index(j) = benam(parts(j))
 !----- parts_index(i) = 0 if nuclide i is not in network
         if (parts_index(j) .eq. 0) cycle outer
      end do inner

      i = i+1

   end do outer

   ! Close the fileagain
   close(reaclib)

   ! Store the amount of rates
   nrea = i

end subroutine count_reaclib_rates







end module reaclib_rate_module
