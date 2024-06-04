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
!!    - 04.10.23  M. Jacobi : Tabulated rates of variable lenghts
!! .
!!
!! @author Moritz Reichert
!! @date 24.01.21
#include "macros.h"
module tabulated_rate_module
   use global_class, only: reactionrate_type
   implicit none

   type,private :: tabulated_rate_type !< type for tabulated reaction rates
                                       !< so far contains only the tabulated reaction rates but can be extended
      real(r_kind),dimension(:),allocatable:: tabulated !< tabulated reaction rates
   end type tabulated_rate_type

   integer, private                :: ntab              !< number of tabulated rates (e.g. calculated with TALYS)
   integer                         :: nt_tab            !< number of temperature grid points,
   logical, public                 :: tabulated         !< switch for tabulated rates
   integer,dimension(2), private   :: tab_index         !< Multi-index for the tabulated rates

   character(len=*), private, parameter                        :: tabulated_binary_name='tabulated_rates.windat' !< Filename of binary file to save weak rates

   real(r_kind), dimension(:), allocatable, private :: temp_grid_tab
   real(r_kind), dimension(30), private             :: temp_grid_tab_default = &
    (/1.0d-4,5.0d-4,1.0d-3,5.0d-3,1.0d-2,5.0d-2,1.0d-1,1.5d-1,2.0d-1,2.5d-1, &
      3.0d-1,4.0d-1,5.0d-1,6.0d-1,7.0d-1,8.0d-1,9.0d-1,1.0d+0,1.5d+0,2.0d+0, &
      2.5d+0,3.0d+0,3.5d+0,4.0d+0,5.0d+0,6.0d+0,7.0d+0,8.0d+0,9.0d+0,1.0d+1 /) !< default Temperature grid of tabulated reaction rates [GK]

   type(reactionrate_type),dimension(:),allocatable,public       :: rrates_tabulated !< array containing all tabulated reaction rates in rrate format
   type(tabulated_rate_type), dimension(:), allocatable,public   :: tabulated_rate   !< array containing all tabulated reaction rates

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
   !!
   !! \b Edited:
   !!   - 04.10.23, M. Jacobi - support for flexible tabulated temperature grids
   !! .
   subroutine init_tabulated_rates()
      use parameter_class, only: use_tabulated_rates, &
                                 tabulated_rates_file, &
                                 use_prepared_network, &
                                 prepared_network_path

      use file_handling_class
      implicit none
      integer :: tab_unit !< File unit id

      ! Give the flag a different name
      tabulated=use_tabulated_rates

      ! Read and count tabulated rates
      ntab = 0
      nt_tab = 0
      if (tabulated) then
        if (use_prepared_network) then
            call read_binary_tabulated_reaction_data(prepared_network_path)
        else
           tab_unit= open_infile(tabulated_rates_file)
           ! readtabulatedtemps reads in the tabulated rate
           ! temperature grid or sets default if no specific
           ! grid is given
           call readtabulatedtemps()
           ! readtabulated returns number of tabulated rates
           call readtabulated(tab_unit,ntab)
         end if
         ! Give a verbose output
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
         call write_data_to_std_out("Size tabulated rate temperature grid",int_to_str(nt_tab))
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
      use parameter_class,  only: use_prepared_network
      implicit none
      type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate_array  !< Large rate array, containing all reactions
      integer,intent(inout)                                          :: rrate_length !< length of rrate_array
      type(reactionrate_type),dimension(:),allocatable               :: rrate_tmp    !< Temporary rate array
      integer                                                        :: alloc_stat   !< Allocation state
      integer                                                        :: new_length   !< New length of rrate_array

      if (tabulated .and. (.not. use_prepared_network)) then
        if (ntab .ne. 0) then
           if (.not. allocated(rrate_array)) then
              rrate_length = ntab
              !-- Allocate the reaclib rate array
              allocate(rrate_array(ntab),stat=alloc_stat)
              if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                         "merge_tabulated_rates",420001)
              rrate_array(1:ntab) = rrates_tabulated(1:ntab)
           else
              !----- merge tabulated rates into rrate
              call rrate_qs_replace(rrate_array(1:rrate_length),rrate_length,rrates_tabulated(1:ntab),ntab,1,rrate_tmp,new_length)
              rrate_length = new_length
              deallocate(rrate_array)
              allocate(rrate_array(rrate_length))
              ! TODO throw exception
              rrate_array(:) = rrate_tmp(:)
           end if
           !-- Deallocate the reaclib rate array
           deallocate(rrates_tabulated)
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
   !!  - 04.10.23, MJ: Added support for flexible tabulated temperature grids
   !! .
   subroutine calculate_tab_rate(rrate, temp, rat_calc)
     use global_class, only: reactionrate_type
     implicit none
     ! Declare the pass
     type(reactionrate_type),intent(in)  :: rrate    !< rate instance
     real(r_kind),intent(in)             :: temp     !< Temperature [GK]
     real(r_kind),intent(out)            :: rat_calc !< rate value

     integer :: tab_rate_index !< index of the tabulated rate

     ! Interpolate the rate
     tab_rate_index = int(rrate%param(1))
     rat_calc = tabulated_inter(tabulated_rate(tab_rate_index)%tabulated,temp)
   end subroutine calculate_tab_rate


   !> Reads tabulated reaction rate temperature grid.
   !!
   !! If tabulated_temperature_file is not given, a default
   !! temperature grid is used.
   !!
   !! @see parameter_class::use_tabulated_rates,
   !!      parameter_class::tabulated_temperature_file
   !!
   !! \b Edited:
   !!     - 04.06.24, MR: - Added check for monotonicity
   !! .
   !!
   !! @author M. Jacobi
   subroutine readtabulatedtemps()
     use error_msg_class, only: int_to_str
     use parameter_class, only: tabulated_temperature_file, max_fname_len
     use file_handling_class

     implicit none

     integer :: i
     integer :: read_stat,alloc_stat
     integer :: tabtemp_unit !< File unit id
     character(max_fname_len) :: help_reader  !< Helper variable

     INFO_ENTRY("readtabulated")

     if (tabulated_temperature_file .eq. "default") then
        nt_tab = 30
        allocate(temp_grid_tab(30),stat=alloc_stat)
        temp_grid_tab = temp_grid_tab_default
        return
     end if

     ! Count how many lines to skip and how many are there
     tabtemp_unit= open_infile(tabulated_temperature_file)

     do ! determine the number of records
        read(tabtemp_unit,"(A)",iostat=read_stat)help_reader
        ! Go out, file ended
        if (read_stat .ne. 0)exit
        help_reader = trim(adjustl(help_reader))
        ! Check if the line is to skip or not
        if ((len_trim(help_reader) .eq. 0) .or. (help_reader(1:1) .eq. "#")) then
          cycle
        else
          exit
        end if
     end do
     close(tabtemp_unit)

     ! Count the number of space separated entries
     nt_tab = 0
     do i = 1, len_trim(help_reader)-1
       if ((help_reader(i:i) == ' ') .and. (help_reader(i+1:i+1) .ne. ' ')) then
         nt_tab = nt_tab + 1
       end if
     end do
     nt_tab = nt_tab + 1  ! Add 1 to account for the last entry

     ! Allocate the temperature grid
     allocate(temp_grid_tab(nt_tab),stat=alloc_stat)

     ! Read the temperature grid
     read(help_reader,*) temp_grid_tab

     ! Check that it is monotonically increasing
     do i = 1, nt_tab-1
        if (temp_grid_tab(i) >= temp_grid_tab(i+1)) then
            call raise_exception("Temperature grid is not monotonically increasing",&
                                "readtabulatedtemps",420004)

        end if
     end do


   end subroutine readtabulatedtemps

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
   !!      - 04.10.23, MJ - support for flexible tabulated temperature grids
   !!      - 31.05.24, MR - Fixed bug related to reading rates.
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
     real(r_kind),dimension(:), allocatable  :: curTable
     real(r_kind)                    :: qvalue
     integer                         :: j,k,read_stat,alloc_stat
     character(20)                   :: fmt_dyn

     INFO_ENTRY("readtabulated")

     allocate(curTable(nt_tab),stat=alloc_stat)

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

        ! Check if reaction is in the network
        select case (grp)
        case (1:11)
           group_index = grp
           cycle
        case default
        parts_index = 0
        tinner_loop_cnt: do j=1,6
           if (parts(j) .eq. '     ') exit tinner_loop_cnt
           parts_index(j) = benam(parts(j))
           !----- parts_index(j)==0 means that nuclide j is not part of the network
           if (parts_index(j) .eq. 0) cycle talloc_loop
        end do tinner_loop_cnt
        !----- if both participants are part of the network, write the rate into
        !----- tabulated_rate
        end select

        k=k+1
     end do talloc_loop

     ntab = k

   !----- allocate the array of tabulated reactions
     allocate(rrates_tabulated(ntab),stat=alloc_stat)
   !----- allocate the array of tabulated rates
     allocate(tabulated_rate(ntab),stat=alloc_stat)
     do k=1, ntab
        allocate(tabulated_rate(k)%tabulated(nt_tab),stat=alloc_stat)
     end do

     if ( alloc_stat /= 0)  call raise_exception('Allocation of "tabulated_rate" failed',&
                                                 "readtabulated",420001)
     rewind(sourcefile)

     k=1
     cntTab=0
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
        rrates_tabulated(k)%group       = group_index
        rrates_tabulated(k)%parts       = parts_index
        rrates_tabulated(k)%source      = src
        rrates_tabulated(k)%q_value     = qvalue
        rrates_tabulated(k)%is_resonant = (res == "r")
        rrates_tabulated(k)%is_weak     = (res == "w")
        rrates_tabulated(k)%is_reverse  = (rev == "v")
        rrates_tabulated(k)%param       = 0.0
        rrates_tabulated(k)%param(1)    = dble(k)
        rrates_tabulated(k)%reac_src    = rrs_tabl
        rrates_tabulated(k)%cached      = -1

        tabulated_rate(k)%tabulated   = curTable

        cntTab=k
        ! Set the reaction type
        call set_reaction_type(rrates_tabulated(k))
        ! Next rate
        k=k+1
     end do touter_loop

     ! Make the reading bullet proof
     if (ntab .ne. cntTab) then
        call raise_exception('Number of tabulated rates does not match while reading!',&
                             "readtabulated",420003)
     end if

   ! get the correct coefficients to prevent double counting
   call getcoefficients(rrates_tabulated,ntab)

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

      real(r_kind),dimension(:) :: rate  !< Rate entries
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
   !! with either 1 or the last index, respectively.
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


  !> Read the tabulated rates from a unformatted binary file
  !!
  !! In case the binary file is read, no other data has to be read.
  !!
  !! @author M. Jacobi
  !! @date 04.10.23
  subroutine read_binary_tabulated_reaction_data(path)
    use file_handling_class, only: open_unformatted_infile
    use error_msg_class,     only: raise_exception
    implicit none
    character(len=*), intent(in)        :: path            !< Path to binary file
    integer                             :: file_id         !< File identifier
    integer                             :: i               !< Loop variable
    integer                             :: status          !< Status of allocation


    file_id = open_unformatted_infile(trim(adjustl(path))//trim(adjustl(tabulated_binary_name)))

    read(file_id) ntab
    read(file_id) nt_tab

   !----- allocate the array containing the temperature grid
     allocate(temp_grid_tab(nt_tab), stat=status)
   !----- allocate the array of tabulated reactions
     allocate(rrates_tabulated(ntab),stat=status)
   !----- allocate the array of tabulated rates
     allocate(tabulated_rate(ntab),stat=status)
     do i=1, ntab
        allocate(tabulated_rate(i)%tabulated(nt_tab),stat=status)
     end do

     if ( status /= 0)  call raise_exception('Allocation of "tabulated_rate" failed',&
                                                 "readtabulated",420001)

     read(file_id) temp_grid_tab

     do i=1, ntab
         read(file_id) tabulated_rate(i)%tabulated
     end do


     close(file_id)

  end subroutine read_binary_tabulated_reaction_data




  !> Save the theoretical tabulated rates to a unformatted binary file
  !!
  !! @author M. Jacobi
  !! @date 04.10.23
  subroutine output_binary_tabulated_reaction_data(path)
    use file_handling_class, only: open_unformatted_outfile

    implicit none
    character(len=*), intent(in) :: path
    integer                             :: i
    integer                             :: file_id

    if (.not. tabulated) return

    ! Open an unformatted file
    file_id = open_unformatted_outfile(trim(adjustl(path))//trim(adjustl(tabulated_binary_name)))
    ! Save the data
    write(file_id) ntab
    write(file_id) nt_tab
    write(file_id) temp_grid_tab

    do i=1,ntab
      write(file_id) tabulated_rate(i)%tabulated
    end do

    close(file_id)

   end subroutine output_binary_tabulated_reaction_data

end module tabulated_rate_module
