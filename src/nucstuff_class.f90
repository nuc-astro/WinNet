!> @file nucstuff_class.f90
!!
!! The error file code for this file is ***W32***.
!!
!! @brief Module \ref nucstuff_class with nuclear physics helpers
!!

!> Module with helper functions such as rate tabulation, partition functions,
!! and theoretical weak rates.
#include "macros.h"
module nucstuff_class
  use global_class, only: t9_data, isotope, net_size

  implicit none
  real(r_kind),dimension(9), public    :: t9_pow       !< Powers of T, used for the rates
                                                       !< @see calc_t9_pow
  integer,public                               :: ntgp !< (24/72) Number of temp grid points for the partition functions
  real(r_kind),dimension(:),allocatable,public :: pf   !< partition functions

     !
     ! Public and private fields and methods of the module.
     !
     public:: &
           calc_t9_pow, inter_partf, masscalc, el_ab, analyze_src_string, is_stable

     ! private:: &
     !     nothing here

contains



!> Initialize nucstuff class
!!
!! Allocate the partition function array.
!!
!! @author M. Reichert
!! @date 26.07.22
subroutine init_nucstuff()
  use global_class,    only: net_size
  use error_msg_class, only: raise_exception
  implicit none
  integer :: istat !< Allocation variable

  ! Allocate the partition function array
  allocate(pf(0:net_size),stat = istat)

  ! Throw exception
  if (istat .ne. 0) then
    call raise_exception('Allocation of "pf" failed.',"init_nucstuff",320001)
  end if

end subroutine init_nucstuff

!> Function to decide whether a given isotope is stable or not
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = is_stable(2,2)
!!~~~~~~~~~~~~~~
!! b will be true.
!!
!! @author M. Reichert
function is_stable(Z,N)
    implicit none
    integer, intent(in) :: Z         !< proton number
    integer, intent(in) :: N         !< neutron number
    logical             :: is_stable !< is the isotope stable?
    integer             :: i         !< loop variable
    integer,parameter   :: numstable=253 !< number of stable isotopes
    integer, dimension(numstable) ::  Nstable = &
        (/ 0,   1,   1,   2,   3,   4,   5,   5,   6,   6,   7,   7,   8,   8,   9,  10,  10,  10, &
          11,  12,  12,  12,  13,  14,  14,  14,  15,  16,  16,  16,  17,  18,  20,  18,  20,  18, &
          20,  22,  20,  22,  20,  22,  23,  24,  26,  24,  24,  25,  26,  27,  28,  28,  26,  28, &
          29,  30,  30,  28,  30,  31,  32,  32,  30,  32,  33,  34,  36,  34,  36,  34,  36,  37, &
          38,  40,  38,  40,  38,  40,  41,  42,  42,  40,  42,  43,  44,  46,  44,  46,  42,  44, &
          46,  47,  48,  50,  48,  46,  48,  49,  50,  50,  50,  51,  52,  54,  52,  50,  52,  53, &
          54,  55,  56,  52,  54,  55,  56,  57,  58,  60,  58,  56,  58,  59,  60,  62,  64,  60, &
          62,  58,  60,  62,  63,  64,  66,  64,  62,  64,  65,  66,  67,  68,  69,  70,  72,  74, &
          70,  72,  68,  70,  71,  72,  73,  74,  74,  70,  72,  74,  75,  76,  77,  78,  80,  78, &
          74,  76,  78,  79,  80,  81,  82,  82,  78,  80,  82,  84,  82,  82,  83,  85,  86,  88, &
          82,  87,  88,  90,  92,  90,  90,  91,  92,  93,  94,  96,  94,  90,  92,  94,  95,  96, &
          97,  98,  98,  94,  96,  98,  99, 100, 102, 100,  98, 100, 101, 102, 103, 104, 106, 104, &
         104, 105, 106, 107, 108, 108, 108, 109, 110, 112, 110, 111, 112, 113, 114, 116, 114, 116, &
         114, 116, 117, 118, 120, 118, 116, 118, 119, 120, 121, 122, 124, 122, 124, 122, 124, 125, &
         126/)
    integer, dimension(numstable) ::  Zstable = &
       (/ 1,  1,  2,  2,  3,  3,  4,  5,  5,  6,  6,  7,  7,  8,  8,  8,  9, 10, 10, 10, 11, 12, 12, 12, &
         13, 14, 14, 14, 15, 16, 16, 16, 16, 17, 17, 18, 18, 18, 19, 19, 20, 20, 20, 20, 20, 21, 22, 22, &
         22, 22, 22, 23, 24, 24, 24, 24, 25, 26, 26, 26, 26, 27, 28, 28, 28, 28, 28, 29, 29, 30, 30, 30, &
         30, 30, 31, 31, 32, 32, 32, 32, 33, 34, 34, 34, 34, 34, 35, 35, 36, 36, 36, 36, 36, 36, 37, 38, &
         38, 38, 38, 39, 40, 40, 40, 40, 41, 42, 42, 42, 42, 42, 42, 44, 44, 44, 44, 44, 44, 44, 45, 46, &
         46, 46, 46, 46, 46, 47, 47, 48, 48, 48, 48, 48, 48, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, &
         51, 51, 52, 52, 52, 52, 52, 52, 53, 54, 54, 54, 54, 54, 54, 54, 54, 55, 56, 56, 56, 56, 56, 56, &
         56, 57, 58, 58, 58, 58, 59, 60, 60, 60, 60, 60, 62, 62, 62, 62, 62, 63, 64, 64, 64, 64, 64, 64, &
         65, 66, 66, 66, 66, 66, 66, 66, 67, 68, 68, 68, 68, 68, 68, 69, 70, 70, 70, 70, 70, 70, 70, 71, &
         72, 72, 72, 72, 72, 73, 74, 74, 74, 74, 75, 76, 76, 76, 76, 76, 77, 77, 78, 78, 78, 78, 78, 79, &
         80, 80, 80, 80, 80, 80, 80, 81, 81, 82, 82, 82, 82 /)

    is_stable = .false.
    if ((N .gt. maxval(Nstable)) .or. (N .lt. minval(Nstable))) then
        is_stable = .false.
    elseif ((Z .gt. maxval(Zstable)) .or. (Z .lt. minval(Zstable))) then
        is_stable = .false.
    else
        do i=1,numstable
            if ((N .eq. Nstable(i)) .and. (Z .eq. Zstable(i))) then
                is_stable = .true.
                exit
            end if
        end do
    end if

end function is_stable




!>
!! @brief Calculates partition function
!!
!! The partition functions are interpolated from
!! the temperature grid (stored in \ref t9_data, read from the winvn in
!! \ref benam_class::get_nuclear_properties (index 1-24) and the htpf file
!! in \ref benam_class::read_htpf (index 25-72)) and the tabulated
!! partition functions \ref global_class::isotope\%part_func.
!!
!! \b Edited:
!!       - MR : 19.01.2021 - set ntgp to 72
!! .
subroutine inter_partf (temp, interpol)
   implicit none

   real(r_kind),intent(in)     :: temp                         !< [in] temperature [GK]
   real(r_kind),dimension(0:net_size),intent(out) :: interpol  !< [out] partition function
   !
   integer                     :: i,k          !< Loop variables
   real(r_kind)                :: grad         !< Gradient
   real(r_kind)                :: lpfk, lpfkl  !< log values of the partition functions at different indices

   interpol(0)=1.d0

   if (temp .gt. t9_data(ntgp)) then
      k=ntgp+1
   else
      do k=1,ntgp
         if (temp .le. t9_data(k)) exit
      end do
   end if

   if (k .eq. 1) then ! The input temperature is lower than the grid
      do i=1,net_size
         interpol(i) = isotope(i)%part_func(1)
      end do
   elseif (k .eq. ntgp+1) then ! The input temperature is hotter than the grid
      do i=1,net_size
         interpol(i) = isotope(i)%part_func(ntgp)
      end do
   else ! The input temperature lies in the grid
      do i=1,net_size
         lpfk  = dlog(isotope(i)%part_func(k))
         lpfkl = dlog(isotope(i)%part_func(k-1))
         grad=(lpfk-lpfkl)/(t9_data(k)-t9_data(k-1))
         interpol(i) = dexp(lpfkl + (temp-t9_data(k-1))*grad)
      end do
   end if

   return
end subroutine inter_partf


!>
!! A helper to compute powers of temperature used in the reaction rates
!!
!! The reaclib reaction rates are calculated with the following fit formula:
!! \f[
!! \lambda = \mathrm{exp} \left[ a_0 + \sum \limits _{i=1} ^5 a_i T_9 ^{\frac{2i-5}{3}} a_6 \, \mathrm{ln} T_9 \right].
!! \f]
!! The helper array \ref t9_pow stores the powers of the temperature used for this calculation
!! (starting with index 2 for i=1, index 1 stores unity). The \ref t9_pow array is after the call given
!! by \f[ \mathrm{t9\_pow} = [ T_9 ^{0}, T_9 ^{-1}, T_9 ^{-1/3}, T_9 ^{1/3}, T_9 ^{1},
!! T_9 ^{5/3}, \log T_9, T_9 ^{7/3}, T_9 ^{9/3} ] \f]
!!
!! @returns Array (\ref t9_pow) that contains the temperature with different exponents
!!
!! \b Edited:
!!       - MR : 19.01.2021 - renamed subroutine and removed rho as input
!! .
subroutine calc_t9_pow(t9)
   implicit none

   real(r_kind),intent(in)    :: t9   !< Temperature [GK]

   t9_pow(1)   = 1.d0               ! t9^(0)
   t9_pow(2)   = 1.d0/t9            ! t9^-1
   t9_pow(4)   = t9**(1.d0/3.d0)    ! t9^(1/3)
   t9_pow(3)   = 1.d0/t9_pow(4)     ! t9^(-1/3)
   t9_pow(5)   = t9                 ! t9^(1)
   t9_pow(6)   = t9_pow(4)**5.d0    ! t9^(5/3)
   t9_pow(7)   = dlog(t9)           ! log(t9)
   t9_pow(8)   = t9_pow(4)**7.d0    ! t9^(7/3)
   t9_pow(9)   = t9_pow(4)**9.d0    ! t9^(9/3)

end subroutine calc_t9_pow


!>
!! Computes the electron fraction
!!
!! The electron fraction is defined as:
!! \f[ Y_e = \sum Z(i) \cdot Y(i) \f],
!! with the number of protons Z and the abundance Y
!! of each isotope i.
!!
!! \b Edited:
!!       - MR : 19.01.2021
!! .
subroutine el_ab (Y, Ye)
   implicit none

   real(r_kind),dimension(:),intent(in) :: Y       !< Abundances
   real(r_kind),intent(out)             :: Ye      !< Electron fraction

   Ye= sum(isotope(1:net_size)%p_nr*Y(1:net_size))

end subroutine el_ab




!> Get the reaclib chapter based on nr. of prods. and educts.
!!
!! This function translates the number of educts and products
!! into the corresponding reaclib chapter.
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = get_chapter(2,2)
!!~~~~~~~~~~~~~~
!! b will be 5 afterwards.
!!
!! @Author M. Reichert
!! @date 01.08.22
function get_chapter(nr_react,nr_prod)
  use error_msg_class, only: raise_exception, int_to_str
  implicit none
  integer,intent(in) :: nr_react    !< Nr. reactants
  integer,intent(in) :: nr_prod     !< Nr. products
  integer            :: get_chapter !< Reaclib chapter


  if     ((nr_react .eq. 1) .and. (nr_prod .eq. 1)) then
    get_chapter = 1
  elseif ((nr_react .eq. 1) .and. (nr_prod .eq. 2)) then
    get_chapter = 2
  elseif ((nr_react .eq. 1) .and. (nr_prod .eq. 3)) then
    get_chapter = 3
  elseif ((nr_react .eq. 2) .and. (nr_prod .eq. 1)) then
    get_chapter = 4
  elseif ((nr_react .eq. 2) .and. (nr_prod .eq. 2)) then
    get_chapter = 5
  elseif ((nr_react .eq. 2) .and. (nr_prod .eq. 3)) then
    get_chapter = 6
  elseif ((nr_react .eq. 2) .and. (nr_prod .eq. 4)) then
    get_chapter = 7
  elseif ((nr_react .eq. 3) .and. (nr_prod .eq. 1)) then
    get_chapter = 8
  elseif ((nr_react .eq. 3) .and. (nr_prod .eq. 2)) then
    get_chapter = 9
  elseif ((nr_react .eq. 4) .and. (nr_prod .eq. 2)) then
    get_chapter = 10
  elseif ((nr_react .eq. 1) .and. (nr_prod .eq. 4)) then
    get_chapter = 11
  else
    call raise_exception("Corresponding Reaclib chapter not known ("//&
                         "nr. educts: "//int_to_str(nr_react)//&
                         ", nr. products: "//int_to_str(nr_prod)//").",&
                         "get_nr_educts",320004)
  end if

  return
end function get_chapter



!> Get number of reactants of the reaction based on the reaclib chapter
!!
!! Returns the number of reactants in the reaction.
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = get_nr_reactants(2)
!!~~~~~~~~~~~~~~
!! b will be 1 afterwards.
!!
!! @see [Reaclib help](https://reaclib.jinaweb.org/help.php?topic=reaclib_format&intCurrentNum=0)
!!
!! @Author M. Reichert
!! @date 01.08.22
function get_nr_reactants(group)
  use error_msg_class, only: raise_exception, int_to_str
  implicit none
  integer,intent(in)  :: group         !< Reaclib chapter
  integer             :: get_nr_reactants !< Amount of educts

  ! Default value
  get_nr_reactants = 0

  select case(group)
    case(1:3,11)
      get_nr_reactants = 1
    case(4:7)
      get_nr_reactants = 2
    case(8:9)
      get_nr_reactants = 3
    case(10)
      get_nr_reactants = 4
    case default
      call raise_exception("Unknown Reaclib chapter ("//int_to_str(group)//").",&
                           "get_nr_reactants",320003)
  end select

  return
end function get_nr_reactants


!> Get number of products of the reaction based on the reaclib chapter
!!
!! Returns the number of products in the reaction.
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = get_nr_products(2)
!!~~~~~~~~~~~~~~
!! b will be 2 afterwards.
!!
!! @warning Chapter 8 may return a different value in the original reaclib format.
!!          There, reactions can also be chapter 9 reactions.
!!
!! @see [Reaclib help](https://reaclib.jinaweb.org/help.php?topic=reaclib_format&intCurrentNum=0)
!!
!! @Author M. Reichert
!! @date 01.08.22
function get_nr_products(group)
  use error_msg_class, only: raise_exception, int_to_str
  implicit none
  integer,intent(in)  :: group         !< Reaclib chapter
  integer             :: get_nr_products !< Amount of educts

  ! Default value
  get_nr_products = 0

  select case(group)
    case(1,4,8)
      get_nr_products = 1
    case(2,5,9,10)
      get_nr_products = 2
    case(3,6)
      get_nr_products = 3
    case(7,11)
      get_nr_products = 4
    case default
      call raise_exception("Unknown Reaclib chapter ("//int_to_str(group)//").",&
                           "get_nr_products",320003)
  end select

  return
end function get_nr_products



!> Total mass used to check the mass conservation
!!
!! This subroutine calculates
!! \f[ \mathrm{m\_tot} = \sum A_i \cdot Y_i \f]
!!
!! @note This subroutine used to calculate the neutron excess as well
!!       and it can be commented in again.
!!
!! \b Edited:
!!      - MR : 20.01.21 - Removed neutron excess
!! .
subroutine masscalc (Y, m_tot)
   implicit none

   real(r_kind), dimension(:), intent(in) :: Y       !< Abundances
   real(r_kind), intent(out)              :: m_tot   !< Sum of mass fractions
   ! real(r_kind), intent(out), optional    :: nexc
   ! real(r_kind)                           :: m_p, m_n

   m_tot= sum(Y(1:net_size)*isotope(1:net_size)%mass)
   ! if(present(nexc)) then
   !    m_p= sum(Y(1:net_size)*isotope(1:net_size)%p_nr)
   !    m_n= sum(Y(1:net_size)*isotope(1:net_size)%n_nr)
   !    nexc= (m_n - m_p)/m_tot
   ! endif
end subroutine masscalc




!> Analyze a string and split it with delimiter ";"
!!
!! This routine is used to analyze the parameters
!! parameter_class::detailed_balance_src_ignore, parameter_class::detailed_balance_src_q_reac,
!! parameter_class::detailed_balance_src_q_winvn. It takes a string and splits it
!! if there is a ";" in it. Afterwards it returns an array with the splitted strings.
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! a = "rath; pkrF;abc"
!! b = "ths8"
!! call analyze_src_string(a,c,d)
!! call analyze_src_string(b,e,f)
!!~~~~~~~~~~~~~~
!! Afterwards, c will be an array containing (/"rath", "pkrF", " abc"/), d will be 3,
!! e will be (/"ths8"/), and f will be 1.
!!
!! @author M. Reichert
!! @date 04.08.22
subroutine analyze_src_string(input_string,output_array,length_output)
    use parameter_class, only: max_fname_len
    implicit none
    character(len=max_fname_len),intent(in)                  :: input_string !< String with src to analyse, separated by ";"
    character(len=4),dimension(:),allocatable,intent(out)    :: output_array !< Array with sources
    integer,intent(out)                                      :: length_output!< Length of sources
    ! Internal variables
    integer             :: i                !< Loop variable
    integer             :: count            !< helper variable
    character(len=50)   :: help_string      !< helper variable
    character,parameter :: delimiter = ";"  !< Delimiter for separating sources

    INFO_ENTRY("analyze_src_string")

    ! Do nothing for an empty string
    if (trim(adjustl(input_string)) .eq. "") then
      length_output = 0
      allocate(output_array(1)) !TODO throw exception if allocation fails
      output_array(:) = ""
      return
    end if

    ! Count the amount of source files given
    length_output = 1
    do i=1,max_fname_len
      if (input_string(i:i) .eq. delimiter) length_output = length_output + 1
    end do
    ! Allocate
    allocate(output_array(length_output)) !TODO throw exception if allocation fails

    help_string = ""
    count = 0
    do i=1,max_fname_len
      if (input_string(i:i) .eq. delimiter) then
        count = count + 1
        output_array(count) = trim(adjustl(help_string))
        output_array(count) = adjustr(output_array(count))
        help_string         = ""
      else
        help_string = trim(adjustl(help_string))//input_string(i:i)
      end if
    end do
    output_array(count+1) = trim(adjustl(help_string))

    INFO_EXIT("analyze_src_string")

    return
  end subroutine analyze_src_string




end module nucstuff_class
