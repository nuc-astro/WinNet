!> @file file_handling_class.f90
!!
!! The error file code for this file is ***W18***.
!! @brief Module \ref file_handling_class
!!

!> Provide some basic file-handling routines
!!
!! @author Darko Mocelj
!! @date   07.05.04
!!
!! \b Edited:
!!   - 12.11.09: Christian Winteler
!!   - 11.01.14: Oleg Korobkin
!!   - 07.11.16: Oleg Korobkin
!!   - 22.01.21: Moritz Reichert
!! .
module file_handling_class
use error_msg_class
implicit none
integer, parameter, private :: min_unit=10      !< Minimum unit number
integer, parameter, private :: max_unit=9999    !< Maximum unit id

!
! Public and private fields and methods of the module
!
public:: &
   get_next_io_unit, open_outfile, open_infile, close_io_file, delete_io_file

contains

!> Finds the next unused unit number
!!
!! In case there is no open slot for a new file
!! this subroutine will call an error.
!! For example if there are more files open than
!! \ref max_unit, this error will get called.
function get_next_io_unit() result(next)
integer     :: next !< the next available unit number
!
integer,save:: last_unit=min_unit - 1 ! previously found open slot
integer     :: fail_count             ! number of failures
logical     :: is_open                ! file status

   fail_count= 0
   next= last_unit + 1
   forever: do
      inquire(unit=next,opened=is_open)
      if (.not. is_open) then
         last_unit=next !found it
         exit forever
      end if
      next= next + 1
      if (next > max_unit) then
         last_unit= min_unit
         if (fail_count<3) then ! attempt reset 3 times
            fail_count= fail_count + 1
            next= min_unit
         else
            ! complain and exit
            call raise_exception("Cannot find open slot to open a file.",&
                                 "get_next_io_unit",&
                                 180003)
         endif
      end if !reset try
   enddo forever

end function get_next_io_unit


!> Shorthand for opening a new file for writing (output file)
!!
!! This function is used to open a file in order
!! to write into it. It will also complain if it is
!! not possible.
function open_outfile(file_name) result(unit_no)
character(len=*), intent(in) :: file_name     !< path to open
integer                      :: unit_no       !< unit number
integer istatus

   unit_no= get_next_io_unit()
   open (unit=unit_no, file=trim(file_name), status="unknown",iostat=istatus)
   if (istatus /= 0) call raise_exception("Could not open file: "//&
                                          trim(adjustl(file_name)),"open_outfile",&
                                          180004)

end function open_outfile


!> Same for reading (input file)
!!
!! This function is used for read input files
!! and complain if it does not work.
function open_infile(file_name) result(unit_no)
character(len=*), intent(in) :: file_name     !< path to open
integer                      :: unit_no       !< unit number
integer istatus

   unit_no= get_next_io_unit()
   open (unit=unit_no, file=trim(file_name), status="old",iostat=istatus)
   if (istatus /= 0) call raise_exception("Could not open file: "//&
                                          trim(adjustl(file_name)),"open_infile",&
                                          180005)

end function open_infile


!> Close an external file
!!
!! This function is used to close a file and
!! complain if it does not work.
subroutine close_io_file(unit_no, file_name)
integer, intent(in)  :: unit_no                     !< unit number
character(len=*), optional, intent(in) :: file_name !< for reporting
integer istatus

   if (unit_no.lt.min_unit .or. unit_no.gt.max_unit) then
      ! Raise exception (out of bounds) depending on the available data
      if (present(file_name)) then
         call raise_exception("Unit value "//int_to_str(unit_no)//" is out of range"//&
                               NEW_LINE("A")//"when trying to close "//&
                               trim(adjustl(file_name)),"close_io_file",&
                               180006)
      else
         call raise_exception("Unit value "//int_to_str(unit_no)//" is out of range",&
                              "close_io_file",&
                              180006)
      end if
   endif

   ! Close the file
   close(unit=unit_no,iostat=istatus)
   ! Check if the file was closed properly
   if (istatus /= 0) then
      ! Raise exception depending on the available data
      if (present(file_name)) then
         call raise_exception("Closing a file, status =  "//int_to_str(istatus)//"."//&
                               NEW_LINE("A")//"Unit number: "//int_to_str(unit_no)//&
                               NEW_LINE("A")//"File name: "//trim(adjustl(file_name)),&
                               "close_io_file",&
                               180007)
      else
         call raise_exception("Closing a file, status =  "//int_to_str(istatus)//"."//&
                               NEW_LINE("A")//"Unit number: "//int_to_str(unit_no),&
                               "close_io_file",&
                               180007)
       endif
   endif

end subroutine close_io_file



!> Delete a file
!!
!! This function is used to delete a file. If the file does not
!! exist and raise_error is present, an error is raised.
!!
!! @author M. Reichert
subroutine delete_io_file(file_name,raise_error)
  implicit none
  character(len=*), optional, intent(in) :: file_name   !< for reporting
  integer                                :: istatus     !< IO status variable
  integer                                :: unit_no     !< Unit number
  logical,optional                       :: raise_error !< Raise error in case of failure

     unit_no= get_next_io_unit()
     open (unit=unit_no, file=trim(file_name), status="old",iostat=istatus)

     if (istatus == 0) then
       ! Delete the file if present
       close(unit_no, status='delete')
     else
       ! Raise an error if necessary
       if (present(raise_error)) then
         call raise_exception("Deleting a file, status =  "//int_to_str(istatus)//".",&
                               "delete_io_file",&
                               180008)
       end if
     end if
end subroutine delete_io_file






end module file_handling_class
