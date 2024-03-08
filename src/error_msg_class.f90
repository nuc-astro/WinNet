!> @file error_msg_class.f90
!!
!! The error file code for this file is ***W16***.
!! @brief Module \ref error_msg_class with error handling routines
!!

!> Error handling routines
!!
!! @author Darko Mocelj, Christian Winteler
!! @date   31.07.2003
!!
!! \b Edited:
!!     - MR 11.01.2020 - Rewrote the module
!! .
#include "macros.h"
module error_msg_class
  implicit none

  private
  integer,parameter :: error_id=0          !< Default standard error unit in fortran

  logical,private   :: write_header_init=.True. !< Init flag for write_data_to_std_out
  logical,public    :: data_creation_mode=.False. !< Flag for rate creation mode
  !
  ! Public fields and methods of the module (every routine is public)
  !
  public:: &
     raise_exception, num_to_str, int_to_str, write_data_to_std_out, str_to_float,&
     str_to_int, write_final_stats_rate_creation
  private::&
     write_header
contains

!> Write the header of the standard output (usually _OUT_)
!!
!! This subroutine writes the first header of the out file.
!! This looks e.g., like:
!! \file{
!! WinNet - Nuclear reaction network
!! \n
!! \n
!! Option                                   :      Value     Unit
!! \-------------------------------------------------------------- }
!!
!! \b Edited:
!!   - 08.08.22, M.R: Implemented git version
!! .
!!
!! @author Moritz Reichert
!! @date 27.01.21
subroutine write_header()
   implicit none
   character(len=40) :: msg       !< helper string
   character(len=62) :: underline !< helper string
   character(len=10) :: value_len !< helper string
   character(len=8)  :: unit_len  !< helper string

   write(*,*) "              WinNet - Nuclear reaction network"
   write(*,*) "             ==================================="
   if (data_creation_mode) then
       write(*,*) "                     (Data creation mode)"
   end if

   write(*,*) ""
   write(*,*) ""
   msg       = "Option"
   value_len = "Value"
   unit_len  = "Unit"
   write(*,"(A)") adjustl(msg)//" : "//adjustr(value_len)//" "//adjustr(unit_len)
   underline = "--------------------------------------------------------------"
   write(*,"(A)") underline
   write_header_init = .False.

   ! Output the Git version if known
#ifdef GTAG
   call write_data_to_std_out("Release",STR(GTAG))
#else
   call write_data_to_std_out("Release","Unknown")
#endif
   ! Output the git hash if known
#ifdef GHASH
   call write_data_to_std_out("Git hash",STR(GHASH))
#else
   call write_data_to_std_out("Git hash","Unknown")
#endif

end subroutine write_header


!> Write the final stats for rate creation
!!
!! This subroutine writes the final stats when
!! only creating a folder with prepared binary data.
!!
!! @author M. Reichert
!! @date 22.07.23
subroutine write_final_stats_rate_creation
    implicit none

    INFO_ENTRY('write_final_stats_rate_creation')

    write(*,'(A)') '----------------------------------'
    write(*,'(A)') 'Finished data creation successfully.'


    INFO_EXIT('write_final_stats_rate_creation')

end subroutine write_final_stats_rate_creation



!> Write data to the standard output (usually _OUT_)
!!
!! This subroutine formats the output to a uniform design.
!! an example line is:
!! \file{
!! Network size                             :       6545         }
!!
!! @author Moritz Reichert
!! @date 27.01.21
subroutine write_data_to_std_out(str_msg,value_str,unit)
   implicit none
   character(len=*),intent(in)          :: str_msg   !< Message oriented to the left
   character(len=*),intent(in)          :: value_str !< Value represented as a string
   character(len=*),intent(in),optional :: unit      !< unit string
   character(len=40)                    :: msg       !< helper string
   character(len=10)                    :: value_len !< helper string
   character(len=8)                     :: unit_len  !< helper string
   character(len=62)                    :: tot_msg   !< Message

   if (write_header_init) call write_header()

   msg       = str_msg
   value_len = value_str
   if (present(unit)) then
      unit_len  = unit
      tot_msg   = adjustl(msg)//" : "//adjustr(value_len)//" "//adjustr(unit_len)
   else
      tot_msg   = adjustl(msg)//" : "//adjustr(value_len)
   end if

   write(*,"(A)") tot_msg

end subroutine write_data_to_std_out



!> Converts a string to an integer
!!
!! If the string is not a valid integer, an error message is raised.
!!
!! @author Moritz Reichert
!! @date 01.06.22
function str_to_int(input_string)
   implicit none
   character(len=*),intent(in) :: input_string       !< String from param file
   integer                     :: str_to_int         !< Converted integer value from input string
   integer                     :: rstat              !< iostat flag

   !< Convert string to integer
   read(input_string,'(I10)',iostat=rstat) str_to_int

   ! Raise an exception if converting does not work
   if (rstat .ne. 0) then
      call raise_exception('Could not parse string "'//trim(adjustl(input_string))//&
                           '". ', &
                           "str_to_int",&
                           160003)
   end if
end function str_to_int

!> Converts a string to a float
!!
!! If the string is not a valid integer, an error message is raised.
!!
!! @author Moritz Reichert
!! @date 01.06.22
function str_to_float(input_string)
   implicit none
   character(len=*),intent(in) :: input_string       !< String from param file
   real(r_kind)                :: str_to_float       !< Converted integer value from input string
   integer                     :: rstat              !< iostat flag

   !< Convert string to integer
   read(input_string,*,iostat=rstat) str_to_float

   ! Raise an exception if converting does not work
   if (rstat .ne. 0) then
      call raise_exception('Could not parse string "'//trim(adjustl(input_string))//&
                           '". ', &
                           "str_to_float",&
                           160004)
   end if
end function str_to_float


!>
!! Converts a given real to a string with format "(1pE10.2)".
!!
!! This function is often used in error messages as it is useful to convert
!! a number to a string in the message directly
!!
!! @author  Moritz Reichert
function num_to_str(num)
   implicit none
      real(r_kind), intent(in)     :: num        !< Input float
      character(:),allocatable     :: num_to_str !< converted float
      real(r_kind)                 :: num_h      !< Helper variable
      character(len=50)            :: out_msg    !< Helper variable
      num_h = num
      write(out_msg,"(1pE10.2)") num_h
      num_to_str = trim(adjustl(out_msg))
end function num_to_str


!>
!! Converts a given integer to a string.
!!
!! This function is often used in error messages as it is useful to convert
!! a number to a string in the message directly.
!!
!! @author  Moritz Reichert
function int_to_str(num)
   implicit none
      integer, intent(in)      :: num        !< Input integer
      character(:),allocatable :: int_to_str !< converted integer
      integer                  :: num_h      !< Helper variable
      character(len=50)        :: out_msg    !< Helper variable
      num_h = num
      write(out_msg,*) num_h
      int_to_str = trim(adjustl(out_msg))! Output a trimmed value
end function int_to_str


!>
!! Raise a exception with a given error message
!!
!! This subroutine is called when some inconsistency occured. It prints
!! a message to the OUT file that an error occured and prints the rest to the
!! standard error unit (usually file ERR in WinNet). This subroutine
!! also terminates the program.
!!
!! @author  Moritz Reichert
subroutine raise_exception(msg,sub,error_code)
   implicit none
      character(len=*), intent(in)           :: msg        !< Error message
      character(len=*), optional, intent(in) :: sub        !< in which subroutine [opt]
      integer, optional, intent(in)          :: error_code !< Errorcode, see \ref error_codes
      character(len=200)                     :: h_msg
      character(len=30)                      :: ecode_msg

      ! Write a message to the OUT file
      write(*,*) "An error occured. Check 'ERR' file for further information."

      ! Create error code
      if (present(error_code)) then
        ecode_msg="ERROR: W"//int_to_str(error_code)
      else
        ecode_msg="ERROR: W000000"
      end if

      ! Give the subroutine if possible
      if (present(sub)) then
         h_msg = "Location: "//sub//"(): "
         write(error_id,*)
      else
         h_msg = "Location: "
      end if

      ! Write the error to the standard error output
      write(error_id,"(A)") trim(adjustl(ecode_msg))
      write(error_id,"(A)")
      write(error_id,"(A)") trim(adjustl(h_msg))
      write(error_id,"(A)")
      write(error_id,"(A)") trim(adjustl(msg))
      write(error_id,*)
      ! Stop the code
      stop "Exiting."

end subroutine raise_exception

end module error_msg_class
