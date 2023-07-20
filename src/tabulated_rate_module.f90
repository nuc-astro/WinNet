!> @file tabulated_rate_module.f90
!!
!! The error file code for this file is ***W32***.
!!
!! Contains the module \ref tabulated_rate_module


!>
!! @brief This module contains everything for the tabulated rates that
!!        can replace reaclib rates.
!!
!! The tabulated rates were implemented by D. Martin. They replace all occurences
!! of the corresponding reaclib rate (including resonances). It is useful
!! to use rates that are, for example, output of the TALYS code.
!!
!! @note The temperature grid for the tabulated rates can be changed here.
!!
!! @see [TALYS](https://tendl.web.psi.ch/tendl_2019/talys.html)
!!
!! \b Edited:
!!    - 17.07.22, M.R., made the maximum temperature grid points a variable (nt_tab)
!! .
!!
!! @author Moritz Reichert
!! @date 24.01.21
#include "macros.h"
module tabulated_rate_module
   use global_class, only: reactionrate_type
   implicit none


   integer, private                :: ntab              !< number of tabulated rates (e.g. calculated with TALYS)
   logical, public                 :: tabulated         !< switch for tabulated rates
   integer,dimension(2), private   :: tab_index         !< Multi-index for the tabulated rates
   type(reactionrate_type),dimension(:),allocatable,public :: tabulated_rate !< array containing all tabulated reaction rates


   integer, parameter                       :: nt_tab = 30  !< Number of temperature grid points,
                                                            !< this has to be still changed manually
                                                            !< in the global clase
   real(r_kind), dimension(nt_tab), private :: temp_grid_tab =                   &
    (/1.0d-4,5.0d-4,1.0d-3,5.0d-3,1.0d-2,5.0d-2,1.0d-1,1.5d-1,2.0d-1,2.5d-1, &
      3.0d-1,4.0d-1,5.0d-1,6.0d-1,7.0d-1,8.0d-1,9.0d-1,1.0d+0,1.5d+0,2.0d+0, &
      2.5d+0,3.0d+0,3.5d+0,4.0d+0,5.0d+0,6.0d+0,7.0d+0,8.0d+0,9.0d+0,1.0d+1 /) !< Temperature grid of tabulated reaction rates [GK]
                                                                               !< @warning This grid is hard-coded and has to be changed here
                                                                               !< in case of a different grid when using tabulated rates
                                                                               !< @see parameter_class::use_tabulated_rates,
                                                                               !< parameter_class::tabulated_rates_file
   !
   ! Public and private fields and methods of the module
   !
   public:: &
      init_tabulated_rates, merge_tabulated_rates, tabulated_index,&
      calculate_tab_rate
   private:: &
      readtabulated,tabulated_inter
contains

   !> Initialize tabulated rates
   !!
   !! This subroutine creates the shorthand flag "tabulated"
   !! to indicate that tabulated rates are used,
   !! reads and counts the tabulated rates
   !!
   !! @author Moritz Reichert
   !! @date 24.01.21
   subroutine init_tabulated_rates()
      use parameter_class, only: use_tabulated_rates, &
                                 tabulated_rates_file
      use file_handling_class
      implicit none
      integer :: tab_unit !< File unit id

      ! Give the flag a different name
      tabulated=use_tabulated_rates

      ! Read and count tabulated rates
      ntab = 0
      if (tabulated) then
         tab_unit= open_infile(tabulated_rates_file)

    !----- readtabulated returns number of tabulated rates
         call readtabulated(tab_unit,ntab)
    !----- Give a verbose output
         call write_reac_verbose_out()
      endif
   end subroutine init_tabulated_rates


   !> Write the verbose output of the reaction rates
   !!
   !! The rates are always counted, for a certain verbose level they
   !! are also printed to the OUT file
   !!
   !! @author M. Reichert
   !! @date 27.01.21
   subroutine write_reac_verbose_out()
      use error_msg_class, only: int_to_str,write_data_to_std_out
      implicit none

      if (VERBOSE_LEVEL .ge. 1) then
         call write_data_to_std_out("Amount tabulated rates",int_to_str(ntab))
      end if
   end subroutine write_reac_verbose_out


   !> Merge tabulated rates into larger rate array.
   !!
   !! This subroutine merges and replaces rates from the input array with
   !! tabulated rates from \ref tabulated_rate.
   !!
   !! @author D. Martin
   !!
   !! \b Edited:
   !!  - M. Reichert (25.01.21): Made it a separate subroutine
   !!  - M. Reichert (05.08.22): Temporary save rates in other array
   !! .
   subroutine merge_tabulated_rates(rrate_array,rrate_length)
      use error_msg_class,  only: raise_exception
      use mergesort_module, only: rrate_qs_replace
      implicit none
      type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate_array  !< Large rate array, containing all reactions
      integer,intent(inout)                                          :: rrate_length !< length of rrate_array
      type(reactionrate_type),dimension(:),allocatable               :: rrate_tmp    !< Temporary rate array
      integer                                                        :: alloc_stat   !< Allocation state
      integer                                                        :: new_length   !< New length of rrate_array

      if (tabulated) then
        if (ntab .ne. 0) then
           if (.not. allocated(rrate_array)) then
              rrate_length = ntab
              !-- Allocate the reaclib rate array
              allocate(rrate_array(ntab),stat=alloc_stat)
              if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                         "merge_tabulated_rates",420001)
              rrate_array(1:ntab) = tabulated_rate(1:ntab)
           else
              !----- merge tabulated rates into rrate
              call rrate_qs_replace(rrate_array(1:rrate_length),rrate_length,tabulated_rate(1:ntab),ntab,1,rrate_tmp,new_length)
              rrate_length = new_length
              deallocate(rrate_array)
              allocate(rrate_array(rrate_length))
              ! TODO throw exception
              rrate_array(:) = rrate_tmp(:)
           end if
           !-- Deallocate the reaclib rate array
           deallocate(tabulated_rate)
           !-- Output the new length
        end if
      end if

   end subroutine merge_tabulated_rates



   !> Calculates the tabulated rate.
   !!
   !! This subroutine serves as an interface between the tabulated_rate_module
   !! and the jacobian class. It interpolates the rate on the grid and returns
   !! the rate value at a given temperature
   !!
   !! \b Edited:
   !!  - 26.07.22, MR: Created this subroutine to be in line with the other
   !!                  reaction rate types.
   subroutine calculate_tab_rate(rrate, temp, rat_calc)
     use global_class, only: reactionrate_type
     implicit none
     ! Declare the pass
     type(reactionrate_type),intent(in)  :: rrate    !< rate instance
     real(r_kind),intent(in)             :: temp     !< Temperature [GK]
     real(r_kind),intent(out)            :: rat_calc !< rate value

     ! Interpolate the rate
     rat_calc = tabulated_inter(rrate%tabulated,temp)
   end subroutine calculate_tab_rate


   !>
   !! Reads tabulated reaction rates.
   !!
   !! The tabulated reaction rate file is given in the same
   !! chapter style as a usual reaclib file, but has different entries
   !! instead of the fit values. The temperature grid of the tabulated reactions
   !! if given with nucstuff_class::temp_grid_tab . This function is only called
   !! if parameter_class::use_tabulated_rates is set to true. An example file
   !! could look like:
   !! \l_file{...
   !!       be7    n  be6                       tabln    -1.06774e+01
   !! 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 1.323e-94 1.437e-76 1.129e-63 5.421e-54 1.868e-46 2.028e-40 2.890e-22 3.813e-13 1.201e-07 5.800e-04 2.562e-01 2.534e+01 1.669e+04 1.341e+06 3.232e+07 3.659e+08 2.502e+09 1.201e+10 }
   !!
   !! @see parameter_class::use_tabulated_rates,
   !!      parameter_class::tabulated_rates_file
   !!
   !! @author D. Martin
   !!
   !! \b Edited:
   !!      - 14.04.15
   !!      - 23.01.21, MR - set the source to always "tabl"
   !!      - 17.07.22, MR - introduced custom reading format,
   !!                       depending on length of the temperature grid
   !! .
   subroutine readtabulated(sourcefile,cntTab)
     use error_msg_class, only: int_to_str
     use format_class
     use benam_class
     use reaclib_rate_module, only: set_reaction_type
     implicit none

     integer, intent(in)  :: sourcefile  !< file to read tabulated rates from
     integer, intent(out) :: cntTab      !< total count of tabulated rates
     !
     integer                         :: grp
     integer                         :: group_index
     character(5), dimension(6)      :: parts
     integer, dimension(6)           :: parts_index
     character(4)                    :: src
     character(1)                    :: res
     character(1)                    :: rev
     real(r_kind),dimension(nt_tab)  :: curTable
     real(r_kind)                    :: qvalue
     integer                         :: j,k,read_stat,alloc_stat
     character(20)                   :: fmt_dyn

     INFO_ENTRY("readtabulated")

     ! Create a custom format, depending on the length
     fmt_dyn = "("//int_to_str(nt_tab)//"f10.3)"

     k=0
   !----- first read the input file in order to determine the number of tabulated reactions
     talloc_loop: do
   !----- read tabulated rate from file
        read(sourcefile,my_format(1), iostat = read_stat)  &
             grp, parts(1:6), src, res, rev, qvalue
        if (read_stat /= 0) exit
        read(sourcefile,fmt_dyn) curTable
        if (grp.ne.0) cycle talloc_loop
        k=k+1
     end do talloc_loop
     ntab = k
   !----- allocate the array of tabulated reactions
     allocate(tabulated_rate(ntab),stat=alloc_stat)
     if ( alloc_stat /= 0)  call raise_exception('Allocation of "tabulated_rate" failed',&
                                                 "readtabulated",420001)
     rewind(sourcefile)

     k=1
   !----- read the input file again and fill the array of tabulated reactions
     touter_loop: do
   !----- read names of participating nuclides and Q-value
        read(sourcefile,my_format(1), iostat = read_stat)  &
             grp, parts(1:6), src, res, rev, qvalue
        if (read_stat /= 0) exit
        read(sourcefile,fmt_dyn) curTable
        select case (grp)
        case (1:11)
           group_index = grp
           cycle
        case default
           parts_index = 0
           tinner_loop: do j=1,6
              if (parts(j) .eq. '     ') exit tinner_loop
              parts_index(j) = benam(parts(j))
   !----- parts_index(j)==0 means that nuclide j is not part of the network
              if (parts_index(j) .eq. 0) cycle touter_loop
           end do tinner_loop
   !----- if both participants are part of the network, write the rate into
   !----- tabulated_rate
        end select
        tabulated_rate(k)%group       = group_index
        tabulated_rate(k)%parts       = parts_index
        tabulated_rate(k)%source      = src
        tabulated_rate(k)%q_value     = qvalue
        tabulated_rate(k)%is_resonant = (res == "r")
        tabulated_rate(k)%is_weak     = (res == "w")
        tabulated_rate(k)%is_reverse  = (rev == "v")
        tabulated_rate(k)%param       = 0.0
        tabulated_rate(k)%tabulated   = curTable
        tabulated_rate(k)%reac_src    = rrs_tabl
        tabulated_rate(k)%cached      = -1

        cntTab=k
        ! Set the reaction type
        call set_reaction_type(tabulated_rate(k))
        ! Next rate
        k=k+1
     end do touter_loop

   ! get the correct coefficients to prevent double counting
   call getcoefficients(tabulated_rate,ntab)

   INFO_EXIT("readtabulated")
   end subroutine readtabulated


   !> Interpolate tabulated rates from the table
   !!
   !! This function uses a lin-log interpolation to calculate
   !! a given reaction rate (contained in \ref parameter_class::tabulated_rates_file)
   !! An example entry of a tabulated rate looks like:
   !!
   !! \l_file{...
   !!       be7    n  be6                       tabln    -1.06774e+01
   !! 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 1.323e-94 1.437e-76 1.129e-63 5.421e-54 1.868e-46 2.028e-40 2.890e-22 3.813e-13 1.201e-07 5.800e-04 2.562e-01 2.534e+01 1.669e+04 1.341e+06 3.232e+07 3.659e+08 2.502e+09 1.201e+10 }
   !!
   !! @note The function contains also a log-log interpolation that is commented out
   !!       and that can be replaced with the current lin-log interpolation.
   !!
   !! @warning This function should be changed in case of a different temperature grid
   !!          for tabulated rates
   !!
   !! @see tabulated_inter, temp_grid_tab, parameter_class::use_tabulated_rates,
   !!      parameter_class::tabulated_rates_file
   !!
   !! @author D. Martin
   function tabulated_inter(rate,temp) result (tabr)
      implicit none

      real(r_kind),dimension(nt_tab) :: rate  !< Rate entries
      real(r_kind)                   :: temp  !< Temperature [GK]
      real(r_kind)                   :: tabr  !< Interpolated rate at temperature

      tabr = 0.0d0
      if(tab_index(1) .eq. tab_index(2)) then
        tabr = rate(tab_index(1))
      else
        if((rate(tab_index(1)).lt.1.d-100) .or. (rate(tab_index(2)).lt.1.d-100)) then
          tabr = 0.0d0
        else
          !!! lin-log
          tabr = rate(tab_index(1))*(rate(tab_index(2))/rate(tab_index(1)))**((temp-temp_grid_tab(tab_index(1)))&
                 /(temp_grid_tab(tab_index(2))-temp_grid_tab(tab_index(1))))
          !!! log-log
          !tabr = rate(tab_index(1))*(rate(tab_index(2))/rate(tab_index(1)))**(dlog10(temp/temp_grid_tab(tab_index(1)))/dlog10(temp_grid_tab(tab_index(2))/temp_grid_tab(tab_index(1))))
        endif
      endif

      return

   end function tabulated_inter


   !>
   !! @brief Set \ref tab_index for a given temperature
   !!
   !! Sets the indices for the tabulated rate interpolation
   !! in \ref tabulated_inter. Temperature values below or above the
   !! temperature grid will result in two equal indices
   !! with either 1 or 30, respectively.
   !!
   !! ### Example
   !!~~~~~~~~~~~~~~.f90
   !! temp = 7e-4
   !! call tabulated_index(temp)
   !! write(*,*) tab_index(1), tab_index(2)
   !!~~~~~~~~~~~~~~
   !! will output 2 and 3.
   !!
   !! @returns \ref tab_index, array with indices for interpolation
   !!
   !! @see tabulated_inter, temp_grid_tab, parameter_class::use_tabulated_rates,
   !!      parameter_class::tabulated_rates_file, readtabulated
   !!
   !! @author D. Martin
   subroutine tabulated_index (temp)
      implicit none
      real(r_kind), intent(in) :: temp !< Temperature [GK]
      integer                  :: i

      if (temp .gt. temp_grid_tab(nt_tab)) then
        tab_index = nt_tab
      elseif (temp .lt. temp_grid_tab(1)) then
        tab_index = 1
      else
        do i=1,nt_tab
          if (temp .gt. temp_grid_tab(i)) cycle
          tab_index(1) = i-1
          tab_index(2) = i
          exit
        enddo
      endif

   end subroutine tabulated_index


end module tabulated_rate_module
