!> @file benam_class.f90
!!
!! The error file code for this file is ***W11***.
!! @brief  Module \ref benam_class for network initialization
!! @author Christian Winteler
!! @date   12.11.2008
!!

!> Subroutines needed to initialise the network
#include "macros.h"
module benam_class
  use error_msg_class, only: raise_exception,write_data_to_std_out,int_to_str
  implicit none
  integer,allocatable,dimension(:,:) :: minmax !< TODO: add description

  !
  ! Public and private fields and methods of the module
  !
  public:: &
     findaz, benam, load_network, get_nuclear_properties, get_minmax, &
     getcoefficients, &
     reaction_string,get_net_name,ident,convert

  private:: &
      sunet_check, read_htpf, read_winvn, read_neutrino_loss


contains

!>
!! Reads isotope names from file 'net_source' and saves them in net_names().
!! Returns number of isotopes in network.
!!
!! This function reads the sunet with the same format as the following example:
!! \file{
!!    n
!!    p
!!    d
!!    t
!!  he3
!!  he4
!!  he6
!!  li6
!!  li7
!!  li8
!!  li9
!! ...}
!!
!! \b Edited:
!!       - 11.01.14
!! .
subroutine load_network(net_source)
   use global_class, only: net_names, net_size, ihe4, ineu, ipro
   use file_handling_class

   character(*),intent(in) :: net_source  !< file where isotopes are specified
   !
   integer                 :: i,sunet,read_stat
   character(5)            :: dummy

   INFO_ENTRY("load_network")

   sunet= open_infile(net_source)

   i=1

!----- read input file to the end to determine number of isotopes in order
!----- to allocate net_names(:)
   do
      read(sunet,"(2a5)",iostat=read_stat) dummy
      if (read_stat /= 0) exit
      i=i+1
   end do
   net_size=i-1
   if (VERBOSE_LEVEL.ge.1) then
      call write_data_to_std_out("Network size",int_to_str(net_size))
   end if

   allocate(net_names(net_size))
   rewind(sunet)

!----- read input file again and write isotope names into net_names(:)
   do i=1,net_size
      read(sunet,"(2a5)") net_names(i)
   end do

   call close_io_file(sunet,net_source)

   !Setup important indices
   ihe4 = benam('  he4')
   ineu = benam('    n')
   ipro = benam('    p')

   INFO_EXIT("load_network")

   return

end subroutine load_network


!>
!! @brief Return a string to represent a given reaction
!!
!! This routine is useful for error messages in case
!! a rate is not working as intented.
!!
!! ### Example
!! For rrate(i) being the decay of ni56
!! ~~~~~~~~~~~.f90
!! a = reaction_string(rrate(i))
!! ~~~~~~~~~~~
!! will return a="ni56 => co56".
!!
!! @author Moritz Reichert
function reaction_string(reac)
   use global_class, only: net_names, reactionrate_type
   implicit none
   type(reactionrate_type),intent(in) :: reac
   character(50) :: reaction_string
   integer       :: i
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

   return
end function reaction_string



!>
!! Returns lowercase of input string. Numbers are not changed.
!!
!! ### Example:
!!~~~~~~~~~~~~.f90
!! a = 'Ti44'
!! call lowercase(a)
!! b = 'This Is A Text 123'
!! call lowercase(b)
!!~~~~~~~~~~~~
!! After this, a will be 'ti44' and b will be 'this is a text 123'
!!
!! @returns The lowercase of the input string
!!
subroutine lowercase(str)

  implicit none

  character(len=*), intent(inout) :: str !< String to be converted to lower case
  integer :: i, del

  del = iachar('a') - iachar('A')
  do i = 1, len_trim(str)
     if (lge(str(i:i),'A') .and. lle(str(i:i),'Z')) then
        str(i:i) = achar(iachar(str(i:i)) + del)
     end if
  end do

  return

end subroutine lowercase


!>
!! Converts isotope names of weak-table format to the ones in reaclib.
!!
!! ### Example:
!!~~~~~~~~~~~~.f90
!! a = '  h1'
!! call convert(a,b)
!! a = ' n01'
!! call convert(a,c)
!!~~~~~~~~~~~~
!! After this, b will be '   p' and c will be '   n'
!!
!! @returns Reaclib names of neutrons and protons for input "h1" or "n01"
!!
!! \b Edited:
!!       - 09.03.2017
!! .
subroutine convert(wnam,rnam)

  implicit none

  character*4,dimension(2)              :: wnam        !< isotope names in weak-table format
  character*5,dimension(2),intent(out)  :: rnam        !< isotope names in reaclib format
  integer                               :: i

  do i=1,2
     call lowercase(wnam(i))
     wnam(i)=adjustr(wnam(i))
     if(wnam(i).eq.'  h1')wnam(i)='   p'
     if(wnam(i).eq.'  H1')wnam(i)='   p'
     if(wnam(i).eq.' n01')wnam(i)='   n'
     if(wnam(i).eq.' N01')wnam(i)='   n'
     rnam(i)=' '//wnam(i)
  end do

  return

end subroutine convert

!> Returns the index number of isotope 'name'
!!
!! \b Edited:
!!     - 11.01.14
!!     - 20.05.19
!! .
function benam(name) result (name_index)
   use global_class, only: net_names,net_size !< names of isotopes in network
   implicit none
   character(len=5),intent(in)   :: name       !< name of the isotope to be 'indexed'
   integer                       :: name_index !< index of the isotope called 'name'
   integer                       :: i          !< Loop variable

   ! Naive implementation, loop through all:
   name_index = 0

   do i=1,net_size
     if (trim(adjustl(name)) .eq. trim(adjustl(net_names(i)))) then
        name_index = i
        exit
     end if
   end do

   ! Use internal fortran functions
   ! name_index = findloc(net_names,value=name,dim=1)
end function benam


!> Getter of net_names, translating indices to a nucleus name
!!
!! The function will also complain in case the index is out of range
!! when the verbose_level is larger than 0.
!!
!! @author Moritz Reichert
!! @date 21.01.21
function get_net_name(index,trimmed)
   use global_class, only: net_names,net_size !< names of isotopes in network
   implicit none
   integer,intent(in)                :: index        !< index of the isotope called 'name'
   logical,intent(in),optional       :: trimmed      !< If present, return a name without whitespaces
   character(:),allocatable          :: get_net_name !< name of the isotope

  ! Check bounds of the function
   if(index.lt.1 .or. (index.gt.net_size)) then
      get_net_name = '     '
   else
      if (.not. present(trimmed)) get_net_name = net_names(index)
      if (      present(trimmed)) get_net_name = trim(adjustl(net_names(index)))
   end if

end function get_net_name


!>
!! Returns the 1/n! factor where n is the number of equal isotopes
!! entering a reaction. in addition the amount by which an isotope is
!! changed in the reaction is saved (+1 / -1)
!!
!! \b Edited:
!!         - 11.01.14
!!         - 28.07.22, MR: implemented counting factor for all chapters (not only original)
!! .
subroutine getcoefficients(rate_array,length_rate_array)
   use global_class, only: reactionrate_type
   integer, intent(in)                                                :: length_rate_array !< Length of the rate array
   type(reactionrate_type),dimension(length_rate_array),intent(inout) :: rate_array !< Input rate array
   integer                      :: i, j, k, l, m !< Loop variables
   real(r_kind)                 :: dcount !< double count for chapter 10


   INFO_ENTRY("getcoefficients")

   do i=1,length_rate_array
      rate_array(i)%one_over_n_fac = 1.d0

      ! added by Marius: for fission reactions the correct coefficients
      ! are allocated directly in the fission subroutines (neutron multiplicity)
      if ((rate_array(i)%source.eq.'fiss').or.(rate_array(i)%source.eq.'ms99') &
         .or.(rate_array(i)%source.eq.'sfis')) cycle

      rate_array(i)%ch_amount      = 0.d0
      select case (rate_array(i)%group)
      case(1:3,11)
!----- reactions of group 1 to 3 only have one isotope in the initial channel
         rate_array(i)%ch_amount(1) = -1.d0
         do j=2,6
!----- if rrate(i)%parts(j)=0, there are no more isotopes in the exit channel
            if (rate_array(i)%parts(j).ne.0) rate_array(i)%ch_amount(j) = 1.d0
         end do
      case(4:7)
!----- reactions of group 4 to 7 have two isotopes in the initial channel
         rate_array(i)%ch_amount(1:2) = -1.d0
         do j=3,6
            if (rate_array(i)%parts(j).ne.0) rate_array(i)%ch_amount(j) = 1.d0
         end do
!----- if the isotopes in the initial channel are identical, 1/(n)!=1/2
         if (rate_array(i)%parts(1).eq.rate_array(i)%parts(2))             &
              rate_array(i)%one_over_n_fac = 1.d0/2.d0
      case(8:9)
!----- reactions of group 8 have 3 isotopes in the initial channel
         rate_array(i)%ch_amount(1:3) = -1.d0
         do j=4,6
            if (rate_array(i)%parts(j).ne.0) rate_array(i)%ch_amount(j) = 1.d0
         end do
!----- if all isotopes in the initial channel are identical, 1/(n)! = 1/6
         if ((rate_array(i)%parts(1).eq.rate_array(i)%parts(2)) .and.    &
              (rate_array(i)%parts(2).eq.rate_array(i)%parts(3))) then
            rate_array(i)%one_over_n_fac = 1.d0/6.d0

!----- if any two isotopes in the initial channel are identical, 1/(n)! = 1/2
!----- since the sequence of isotopes is ordered, 1 and 3 are only identical if
!----- all 3 isotopes are identical!
         else if ((rate_array(i)%parts(1).eq.rate_array(i)%parts(2)) .or.  &
              (rate_array(i)%parts(2).eq.rate_array(i)%parts(3)) .or.      &
              (rate_array(i)%parts(1).eq.rate_array(i)%parts(3))) then
            rate_array(i)%one_over_n_fac = 1.d0/2.d0
         end if
       case(10)
 !----- reactions of group 8 have 4 isotopes in the initial channel
         rate_array(i)%ch_amount(1:4) = -1.d0
         do j=5,6
            if (rate_array(i)%parts(j).ne.0) rate_array(i)%ch_amount(j) = 1.d0
         end do

         ! MR: implemented the following as if statements would become
         !     quite messy for the case with 4 reactants
         dcount = 1
         do j=1,4
           do k=j+1,4
             ! Two are equal
             if (rate_array(i)%parts(k)==rate_array(i)%parts(j)) then
               dcount = dcount+1
             end if
             do m=k+1,4
               ! Three are equal
               if ((rate_array(i)%parts(k)==rate_array(i)%parts(j)) .and. &
                   (rate_array(i)%parts(k)==rate_array(i)%parts(m))) then
                   dcount = dcount+2
               end if
               do l=m+1,4
                 ! Four are equal
                 if ((rate_array(i)%parts(k)==rate_array(i)%parts(j)) .and. &
                     (rate_array(i)%parts(k)==rate_array(i)%parts(m)) .and. &
                     (rate_array(i)%parts(k)==rate_array(i)%parts(l))) then
                     dcount = dcount+9
                 end if
               end do
             end do
           end do
         end do
         ! Set the factor
         rate_array(i)%one_over_n_fac = 1.d0/dcount
      case default
        ! TODO raise exception
        continue
        ! call raise_exception
      end select
   end do

   INFO_EXIT("getcoefficients")

end subroutine getcoefficients


!>
!! Reads nuclear properties (name,a,n,z,spin,mass excess,partition functions)
!! from file 'winvn' and htpf file and writes them into isotope(:).
!! Also the arrays of the partition function grid and temperature grid are
!! initialized here.
!!
!! \b Edited:
!!      - 11.01.14
!!      - 11.02.21 - MR: Implemented a switch for htpf functions and
!!                       moved things to additional subroutine
!! .
subroutine get_nuclear_properties()
   use parameter_class, only: isotopes_file, htpf_file,use_htpf, &
                              use_neutrino_loss_file, neutrino_loss_file
   use nucstuff_class,  only: ntgp, is_stable
   use global_class,    only: net_size, isotope, t9_data
   use file_handling_class
   implicit none
   integer                     :: winvn, htpf !< File IDs
   integer                     :: alloc_stat  !< Allocation status variable
   integer                     :: i           !< Loop variable

   INFO_ENTRY("get_nuclear_properties")

   ! Set the number of temperature grid points
   ! where the partition functions are tabulated.
   ! In case high temperature partition functions are included
   ! more grid points are used.
   if (.not. use_htpf) then
      ntgp = 24
   else
      ntgp = 72
   end if

   ! Allocate grid data array, later set by read_winvn and read_htpf
   allocate(t9_data(ntgp),stat=alloc_stat)
   if (alloc_stat .ne. 0) call raise_exception("Could not allocate partition function "//&
                                               "temperature grid","get_nuclear_properties",&
                                               110001)

   ! Allocate the grid for the partition functions in the correct size
   do i=1,net_size
      allocate(isotope(i)%part_func(ntgp),stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Could not allocate grid for partition functions."&
                                                 ,"get_nuclear_properties",&
                                                 110001)
   end do

   ! Always read winvn
   winvn= open_infile(isotopes_file)
   call read_winvn(winvn)
   call close_io_file(winvn,isotopes_file)

   ! In case high temperature partition functions are included
   if (use_htpf) then
      htpf = open_infile(htpf_file)
      call read_htpf(htpf)
      call close_io_file(htpf,htpf_file)
   end if

   ! Check if the isotope is stable
   do i=1,net_size
       isotope(i)%is_stable = is_stable(isotope(i)%p_nr,isotope(i)%n_nr)
   end do

   ! Read nuloss file
   if (use_neutrino_loss_file) then
      call read_neutrino_loss(neutrino_loss_file)
   end if


   INFO_EXIT("get_nuclear_properties")

end subroutine get_nuclear_properties



!>
!! Reads nuclear properties (name,a,n,z,spin,mass excess,partition functions)
!! from file 'winvn' and writes them into isotope(:).
!! An example entry of a winvn can look like:
!! \l_file{
!! ne33      33.000  10  23   1.5    45.997 reac1
!!  1.00000E+0  1.00000E+0  1.00000E+0  1.00000E+0  1.00000E+0  1.00000E+0  1.00000E+0  1.00000E+0
!!  1.00000E+0  1.00000E+0  1.00000E+0  1.00000E+0  1.00000E+0  1.00000E+0  1.00000E+0  1.01000E+0
!!  1.01000E+0  1.02000E+0  1.03000E+0  1.07000E+0  1.12000E+0  1.19000E+0  1.28000E+0  1.40000E+0
!! ...}
!!
!! \b Edited:
!!    - 11.02.21 - MR: Moved it from get_nuclear_properties to this subroutine
!! .
subroutine read_winvn(winvn)
   use format_class
   use global_class, only: isotope, t9_data
   implicit none
   integer,intent(in)          :: winvn !< File ID of the winvn
   integer                     :: i, j , read_stat
   character(5)                :: name, tempname
   character(5),dimension(:),allocatable    :: data_names
   integer                     :: data_cnt
   real(r_kind)                :: za
   integer                     :: zp  !< Proton number
   integer                     :: zn  !< Neutron number
   real(r_kind)                :: sp  !< Spin of the ground state
   real(r_kind)                :: bi  !< binding energy
   real(r_kind), dimension(24) :: pf  !< Partition function grid

   INFO_ENTRY("read_winvn")
   i = 1

   read (winvn,*)
   read (winvn,"(23f3.2,f3.1)") (t9_data(j),j=1,24)
!----- skip part of winvn where all isotopes are listed. the last isotope in
!----- this list appears twice.
   read (winvn,"(a5)") name                            !read isotope name
   do
      read (winvn,"(a5)",iostat = read_stat) tempname   !read next isotope name
      if (read_stat /= 0) call raise_exception("Problem reading winvn file. Check the winvn.",&
                                               "get_nuclear_properties",&
                                               110003)


!----- if last two isotope names are identical, end of list reached
      if (tempname == name) exit
      name = tempname
      i = i+1
   end do
   data_cnt = i

   allocate(data_names(data_cnt))

   rewind (winvn)
   read(winvn,'(/)')
   do i=1,data_cnt
      read(winvn,"(a5)") data_names(i)
   end do
   read(winvn,*)
!    data_names = adjustr(data_names)
   call sunet_check (data_names)
   deallocate (data_names)
!----- reading of nuclear properties starts here:
   i=0
   do
      read (winvn,my_format(4),iostat = read_stat)name,za,zp,zn,sp,bi
      if (read_stat /= 0) exit
!      read (winvn,my_format(5)) (pf(j),j=1,24)
      read (winvn,*) (pf(j),j=1,24)
      ! Special cases of aluminium isomers
      if ((trim(adjustl(name)) .eq. "al26")) then
        i = benam(" al*6")
        if  (i .ne. 0) then
          isotope(i)%name      = " al*6"
          isotope(i)%mass      = zp+zn
          isotope(i)%p_nr      = zp
          isotope(i)%n_nr      = zn
          isotope(i)%spin      = sp
          isotope(i)%mass_exc  = bi
          isotope(i)%part_func(1:24) = pf
        end if
        i = benam(" al-6")
        if  (i .ne. 0) then
          isotope(i)%name      = " al-6"
          isotope(i)%mass      = zp+zn
          isotope(i)%p_nr      = zp
          isotope(i)%n_nr      = zn
          isotope(i)%spin      = sp
          isotope(i)%mass_exc  = bi
          isotope(i)%part_func(1:24) = pf
        end if
      end if

      ! Normal cases, i.e., not Al26
      i = benam(name)
      if (i == 0) cycle
      isotope(i)%name      = name
      isotope(i)%mass      = zp+zn
      isotope(i)%p_nr      = zp
      isotope(i)%n_nr      = zn
      isotope(i)%spin      = sp
      isotope(i)%mass_exc  = bi
      isotope(i)%part_func(1:24) = pf

   end do

   ! if (VERBOSE_LEVEL .gt. 1) then
   !    do i=1,size(net_names)
   !       if(isotope(i)%name.eq.'  c12')cycle
   !       ! MR: the following may not be an error, it just means that
   !       ! the binding energy per nucleon is equal to c12? I commented out the whole
   !       ! check
   !       if(isotope(i)%mass_exc .eq. 0.d0) then
   !          write(*,*) "Warning: "//trim(adjustl(net_names(i)))//" has a mass excess"//&
   !                     "of zero. Your winvn may be faulty, check it!"
   !       end if
   !    end do
   ! end if

   INFO_EXIT("read_winvn")

end subroutine read_winvn



!>
!! Checks if sunet file contains valid isotope names
!!
!! \b Edited:
!!        - 11.01.14
!! .
subroutine sunet_check(ref_array)
   use global_class, only:net_names

   character(5),dimension(:),intent(in)  :: ref_array !< array of isotope names
   integer,dimension(:),allocatable      :: chk       !< Check variable, 1: included, 0: not included
   integer         :: ref_len, net_len
   integer         :: i,j           !< Loop variale
   character(500)  :: e_message     !< Error message
   character(5)    :: net_name_tmp  !< Temporary helper variable

   INFO_ENTRY("sunet_check")

   ref_len = size(ref_array)
   net_len = size(net_names)

   allocate(chk(net_len))
   chk=0

   do i=1,net_len

      ! Take care of isomers
      if ((trim(adjustl(net_names(i))) .ne. "al*6") .and. &
          (trim(adjustl(net_names(i))) .ne. "al-6")) then
          net_name_tmp = net_names(i)
      else
          net_name_tmp = "al26"
      end if

      in:do j=1,ref_len
         if(trim(adjustl(net_name_tmp)).eq.trim(adjustl(ref_array(j))))then
            chk(i)=1
            exit in
          end if
      end do in
   end do

   if(any(chk.eq.0)) then
      e_message = ""
      do i=1,net_len
         if(chk(i).eq.0) e_message = trim(adjustl(e_message))//net_names(i)//&
                                     " not in isotopes database (winvn)."//NEW_LINE("A")
      end do
      call raise_exception(trim(adjustl(e_message))//"Problem reading sunet file. Check the sunet.",&
                           "sunet_check",&
                           110004)
   end if

   INFO_EXIT("sunet_check")
   return

end subroutine sunet_check


!>
!! Read the file with high-temperature partition function
!! (normally datafile2.dat). An example of the file may look like:
!!
!! \l_file{
!! Title: Nuclear Partition Functions at Temperatures Exceeding 10^10^ K
!! Author: Rauscher T.
!! Table: Renormalized partition functions G(T_9_) including high temperature
!!        corrections. The values given here were calculated with level densities
!!        based on FRDM input (see text).
!! \================================================================================
!! Byte-by-byte Description of file: datafile2.txt
!! \--------------------------------------------------------------------------------
!!    Bytes Format Units Label  Explanations
!! \--------------------------------------------------------------------------------
!!    1-  5 A5     ---   Name   Nuclide name
!!    7-  8 I2     ---   Z      Charge number of nuclide
!!   10- 12 I3     ---   A      Mass number of nuclide
!!   15- 17 F3.1   ---   J0     Ground-state spin of nuclide
!!   19- 27 E9.2   ---   G12    Renormalized partition function for T = 12E9
!!   ....
!!   o13    8  13  1.5 1.00E+000 1.00E+000 9.99E-001 9.95E-001 9.90E-001 9.84E-001 9.76E-001 9.68E-001 9.59E-001 ...
!! ...}
!!
!! \b Edited: 11.01.14
subroutine read_htpf(source)
   use global_class, only: t9_data, isotope

   integer,intent(in)         :: source  !< HTPF file name
   !
   integer                    :: n, i
   character*5                :: namtmp
   character*69               :: dummy
   real*8,dimension(48)       :: pftmp
   integer                    :: stat
   integer,dimension(:),allocatable      :: check

   INFO_ENTRY("read_htpf")

   n = size(isotope)
   allocate(check(n))

   check = 0

   read(source,'(13/)')
   do i=25,72
      read(source,'(a69,es5.0)')dummy,t9_data(i)
   end do

   read(source,*)
   out:do
      read(source,'(a5,t19,48es10.2)',iostat=stat)namtmp,pftmp
      if(stat.ne.0) exit out
      namtmp = adjustr(namtmp)
      in:do i=1,n
         if(namtmp.eq.isotope(i)%name) then
            isotope(i)%part_func(25:72) = pftmp(1:48)
            check(i) = 1
            cycle out
         end if
      end do in
   end do out

   do i = 1,n
      if(check(i).ne.1) then
         isotope(i)%part_func(25:72) = isotope(i)%part_func(24)
      end if
   end do
   deallocate(check)

   INFO_EXIT("read_htpf")

   return

end subroutine read_htpf


!>
!! Finds an isotope index in the network table, given its A and Z.
!! In case it was not found, -1 is returned.
!!
!!
!! \b Edited: 11.01.14
function findaz(ai,zi) result (azind)
   use global_class, only:isotope,net_size

   integer,intent(in)  :: ai       !< A
   integer,intent(in)  :: zi       !< Z
   integer             :: azind    !< return index
   integer :: i

   azind = 0
   do i=1,net_size
      if ((isotope(i)%p_nr.eq.zi) .and. (isotope(i)%mass.eq.ai)) then
         azind = i
         exit
      end if
   end do

   if (azind.eq.0) then
      azind = -1
   end if

   return

end function findaz

!>
!! Returns Amin and Amax for each isotopic chain in minmax
!!
!! \b Edited: 11.01.14
!! .
subroutine get_minmax()
   use global_class, only: isotope, net_size
   !
   integer                       :: i
   integer                       :: ind
   integer                       :: max_p

   INFO_ENTRY("get_minmax")

   !Find maximum proton number
   max_p= 0
   do i= 1,net_size
      if (isotope(i)%p_nr .gt. max_p) max_p = isotope(i)%p_nr
   end do


   allocate(minmax(max_p,2))

   minmax = 0

   do i = 1,net_size
      if (isotope(i)%p_nr.eq.0) cycle
      ind = isotope(i)%p_nr

      if (minmax(ind,1).eq.0) then
         minmax(ind,1) = isotope(i)%mass
         minmax(ind,2) = isotope(i)%mass
      else if (isotope(i)%mass.gt.minmax(ind,2)) then
         minmax(ind,2) = isotope(i)%mass
      else if (isotope(i)%mass.lt.minmax(ind,1)) then
         minmax(ind,1) = isotope(i)%mass
      end if
   end do

   INFO_EXIT("get_minmax")

   return

end subroutine get_minmax


!> Read a file with neutrino loss energy.
!!
!! This file contains the neutrino loss energy for the beta decays in the network.
!! The information is used in case nuclear heating is turned on.
!! An example of the file could look like this:
!! \file{...
!! al34     4.803857E+00
!! al35     6.072946E+00
!! al37     5.639200E+00
!! si25     4.794861E+00
!! si26     1.968612E+00
!! si27     2.070827E+00
!! si31     8.963707E-01
!! si32     1.581000E-01
!! si34     1.778200E+00
!! si35     3.314272E+00
!! si36     3.432230E+00
!! si38     3.446900E+00
!!... }
!! The first column is the isotope name, the second column is the average neutrino loss
!! energy in MeV.
!!
!! @author Moritz Reichert
!! @date   31.03.2023
subroutine read_neutrino_loss(filename)
    use file_handling_class
    use global_class,    only: net_size,Qnuloss,isotope
    implicit none
    character(len=*), intent(in) :: filename !< Name of the file to read
    !-- Internal Variables
    real(r_kind)     :: Qnuloss_tmp  !< Temporary variable for reading
    character(len=5) :: name         !< Isotope name
    integer          :: i            !< Loop counter
    integer          :: stat         !< Status of the file read and allocation
    integer          :: file_id      !< File ID

    INFO_ENTRY("read_neutrino_loss")

    ! Allocate the array
    allocate(Qnuloss(net_size), stat=stat)

    ! Allocation status, complain if something is wrong
    if (stat .ne. 0) then
        call raise_exception("Could not allocate Qnuloss array", &
                             "read_neutrino_loss",110001)
    end if
    ! Set to empty values
    Qnuloss(:) = -1

    ! Open and read the file
    file_id = open_infile(filename)
    do
        read(file_id,*,iostat=stat) name, Qnuloss_tmp
        if (stat .ne. 0) exit
        do i = 1, net_size
            if (trim(adjustl(name)) .eq. trim(adjustl(isotope(i)%name))) then
                Qnuloss(i) = Qnuloss_tmp
                exit
            end if
        end do
    end do

    ! Close the file again
    close(file_id)

    INFO_EXIT("read_neutrino_loss")

end subroutine read_neutrino_loss


!>
!! Identifies the nuclide names and indices corresponding to z,n combinations
!! (z1,n1), (z2,n2).
!!
!! If any of the two nuclides is not part of the network,
!! err=1 and all other properties are set to 0 or "    " respectively.
!!
!!
!! \b Edited:
!!       - 12.01.14
!!       - 23.01.21, MR - moved it to benam_class
!! .
subroutine ident(z1,n1,z2,n2,nam1,nam2,ind1,ind2,err)
  use global_class

  implicit none

!----- subroutine arguments
  integer                    :: z1,z2,n1,n2 !< proton and neutron numbers
  character(5), intent(out)  :: nam1,nam2   !< nuclide names as they appear in sunet
  integer, intent(out)       :: ind1,ind2   !< nuclide index in the network
  integer, intent(out)       :: err         !< 0 if nuclides are identified successfully
                                            !! 1 if any of the two (z,n) combinations
                                            !!   could not be assigned to a nuclide
                                            !!   in the network

!----- runtime variables
  integer                    :: i,j,ini
!----- dummy variables
  integer                    :: zt,nt,indt        ! temporary values
  character(5)               :: namt
  integer                    :: inv

  err = 0
  ind1 = 0
  ind2 = 0
  nam1 = "     "
  nam2 = "     "
  if (z2.lt.z1) then
     zt = z2
     nt = n2
     z2 = z1
     n2 = n1
     z1 = zt
     n1 = nt
     inv = 1
  else
     zt=z1
     nt=n1
     inv = 0
  end if
  out: do j=1,2
     in: do i=zt,net_size
        if((i.gt.3).and.(isotope(i)%p_nr.gt.zt)) exit out
        if((isotope(i)%p_nr.eq.zt).and.(isotope(i)%n_nr.eq.nt)) then
           select case(j)
           case(1)
              ind1 = i
              nam1 = isotope(i)%name
           case(2)
              ind2 = i
              nam2 = isotope(i)%name
           end select
           exit in
        end if
     end do in
     zt = z2
     nt = n2
     ini = ind1
  end do out

  if (ind2.ne.0) then
     if (inv .eq. 1) then
        namt = nam1
        indt = ind1
        nam1 = nam2
        ind1 = ind2
        nam2 = namt
        ind2 = indt
        return
     else
        return
     end if
  else
     ind1 = 0
     nam1 = "     "
     err = 1
     return
  end if

end subroutine ident



!===============================================================================

end module benam_class
