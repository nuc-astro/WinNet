!> @file readini.f90
!!
!! The error file code for this file is ***W39***.
!!
!! Contains the module \ref readini
!!

!> Subroutines for initialization and parameter parsing.
!!
!! The module is responsible for reading several input files
!! such as initial seeds via \ref read_seed_urs or the hydrodynamic trajectory
!! via \ref custom_read
!!
#include "macros.h"
module readini
  use file_handling_class
  use error_msg_class
  use parameter_class, only  : trajectory_format, trajectory_mode,track_nuclei_every,&
                               h_track_nuclei_every,&
                               nuflag,custom_snapshots,h_custom_snapshots,neutrino_mode
  implicit none

  ! Only make necessary routines public
  public  :: readini_init, readini_finalize, read_seed
  private :: read_custom_snapshots, read_track_nuclei, time_unit_conversion, &
             temp_unit_conversion, dens_unit_conversion, dist_unit_conversion, &
             nutemp_unit_conversion, nuenergy_unit_conversion, &
             nulumin_unit_conversion, consistency_check, custom_read

  contains

!>
!! Initialize the readini module and open files.
!!
!! Depending on user defined parameters
!! (e.g., \ref parameter_class::custom_snapshots)
!! this subroutine calls several subroutines to read
!! input files. After reading the files, it calls
!! \ref consistency_check to make a couple of checks and
!! raise an error if something is not correct.
!!
!! @author: Moritz Reichert
!!
!! \b Edited: 28.02.17
subroutine readini_init()
  implicit none

   ! read trajectory
   if (trim(trajectory_mode).EQ.'from_file') then
         call custom_read()
   end if

   ! Custom snapshots
   if (custom_snapshots .or. h_custom_snapshots) then
      call read_custom_snapshots()
   end if

  ! Read the tracked nuclei
  if ((track_nuclei_every .gt. 0) .or. (h_track_nuclei_every .gt. 0)) then
     call read_track_nuclei()
  end if

  ! Check for consistency
  call consistency_check()


end subroutine readini_init



!>
!! @brief Read in times for snapshots.
!!
!! These times should be given by a separate file
!! with the parameter \ref parameter_class::custom_snapshots.
!! The units in this file is days. The subroutines that are based on
!! the custom snapshots assume a sorted time array and it is therefore
!! sorted in the end of this subroutine. An example of the file could be:
!! \file{
!! 1.1574074074074073e-04
!! 0.0208333
!! 1  }
!! This would output a snapshot at 10 s, 30 min, and 1 day.
!!
!! @see timestep_module::restrict_timestep, analysis::output_iteration
!!
!! @author Moritz Reichert
!!
subroutine read_custom_snapshots()
  use parameter_class, only: snapshot_file
  use analysis,        only: snapshot_time,snapshot_amount
  implicit none
  integer                       :: c_snapshots !< File id
  integer                       :: rstat       !< Reading status
  integer                       :: file_length !< Amount of lines of the file
  integer                       :: i           !< Loop variable
  logical                       :: sorted      !< For bubblesort
  real(r_kind)                  :: helper      !< Storage for sorting array

  ! Open the file
  c_snapshots= open_infile(snapshot_file)

  ! Get the length of the file
  file_length = 0
  do
     read(c_snapshots,*,iostat=rstat)
     if (rstat .ne. 0)exit
     file_length = file_length + 1
  end do

  ! Allocate Array
  allocate(snapshot_time(file_length),stat=rstat)
  if (rstat /= 0) call raise_exception('Allocation of "snapshot_time" failed.',"read_custom_snapshots",&
                                      390001)
  ! Rewind the file to read it again
  rewind(c_snapshots)

  ! Read the file again
  do i=1,file_length
     read(c_snapshots,*,iostat=rstat) snapshot_time(i)
     if (rstat .ne. 0) call raise_exception("Unable to read custom snapshots."//NEW_LINE("A")//&
                                            'Check the input file given with the parameter "snapshot_file".',&
                                            "read_custom_snapshots",&
                                            390003)
     snapshot_time(i) = snapshot_time(i)*24.*60.*60. ! Convert Days -> Seconds
  end do

  ! Store also length of the array
  snapshot_amount = file_length
  ! Bubblesort to ensure array is sorted
  ! It is also the fastest way for already sorted arrays
  sorted = .false.
  do while (.not. sorted)
     sorted = .true.
     ! Sort the array. Use Bubblesort
     do i=1,file_length-1
        if (snapshot_time(i+1) .lt. snapshot_time(i)) then
           ! Swap entries
           helper             = snapshot_time(i+1)
           snapshot_time(i+1) = snapshot_time(i)
           snapshot_time(i)   = helper
           sorted = .false.
        end if
     end do
  end do

  ! Close the file again
  call close_io_file(c_snapshots,snapshot_file)
  return

end subroutine read_custom_snapshots






!>
!! @brief Reads in nuclei that will be tracked and where the abundance will
!! be stored in a separate file.
!!
!! If parameter_class::track_nuclei_every is larger than 0,
!! a file will be created, containing the information of the abundances of the nuclei.
!! This is based on the input parameter \ref parameter_class::track_nuclei_file
!! that contains a path to a file with the format of the nuclei being the same
!! as the format of the sunet file. An example could be:
!! \file{
!!  c13
!! ti44
!! ni56 }
!!
!! @note Verbose level greater than 2 will create the file _debug_track_nuclei.dat_
!! to debug this subroutine.
!!
!! @author Moritz Reichert
!!
!! \b Edited:
!!     - 09.06.17
!! .
subroutine read_track_nuclei()
   use parameter_class, only: track_nuclei_file
   use global_class, only: track_nuclei_indices,track_nuclei_nr
   use benam_class, only: benam
   implicit none
   integer              :: track_id       !< File ID
   integer              :: istat          !< Status of reading and allocation
   integer              :: i              !< Loop variable
   integer              :: i_tmp          !< temporary index of the loop
   integer              :: file_length    !< Loop variable
   integer              :: ind_tmp        !< Temporary index of the read nucleus
   integer              :: debug_track    !< File ID of debug file
   character*5          :: nam_tmp        !< Temporary name

   INFO_ENTRY("read_track_nuclei")

   ! Debugging
   if (VERBOSE_LEVEL .gt. 2) then
     ! Open debug file and write a header for the subroutine variance (???)
     debug_track = open_outfile('debug_track_nuclei.dat')
     write(debug_track,'(A)') &
           '# Debugging info for read_track_nuclei() from readini.f90'
     write(debug_track,'("Opened File: ",A)') track_nuclei_file
   end if

   ! Amount of nuclei to track
   track_nuclei_nr = 0
   ! File length can be larger, nuclei are not necessary contained in the list
   file_length     = 0

   ! Open the file
   track_id= open_infile(track_nuclei_file)

   do ! determine the number tracked nuclei
      read(track_id,'(a5)',iostat=istat) nam_tmp
      if (istat .ne. 0)exit

      ind_tmp = benam(nam_tmp)
      ! Check that the nucleus is also contained, otherwise don't track it
      if (ind_tmp .ne. 0) track_nuclei_nr= track_nuclei_nr + 1
      file_length = file_length + 1
   end do

   ! Allocate the array that will contain the indices of the nuclei
   allocate(track_nuclei_indices(track_nuclei_nr),stat=istat)
   if (istat /= 0) call raise_exception('Allocation of "track_nuclei_indices" failed.',"read_track_nuclei",&
                                       390001)

   ! start at the beginning of the file
   rewind(track_id)

   ! Index of the array
   i_tmp = 1
   ! Read the nuclei and store the index
   do i=1,file_length
      read(track_id,'(a5)',iostat=istat)nam_tmp
      ! Get the index of the nucleus
      ind_tmp = benam(nam_tmp)

      if (VERBOSE_LEVEL .gt. 2) then
        !Track which nucleus is read
        write(debug_track,'("Reading",A,7x,"Index: ",I10)') nam_tmp,ind_tmp
      end if
      ! write it to the array
      if (ind_tmp .ne. 0) then
         track_nuclei_indices(i_tmp) = ind_tmp
         i_tmp = i_tmp +1
      else
         ! Say something if there is a nucleus that is not known
         if (VERBOSE_LEVEL .gt. 2) then
            write(debug_track,'("Nucleus",2x,A,2x,"is not contained in network")') nam_tmp
         end if
      end if
   end do

   ! Close everything again
   call close_io_file(track_id,track_nuclei_file)

   if (VERBOSE_LEVEL .gt. 2) then
     write(debug_track,'("Tracking ",I10," nuclei")') track_nuclei_nr
     call close_io_file(debug_track,'debug_track_nuclei.dat')
   end if

   INFO_EXIT("read_track_nuclei")

   return
end subroutine read_track_nuclei



!>
!! @brief Reads in seed abundances from file \ref parameter_class::seed_file
!!
!! The input file can have a custom format specified by \ref parameter_class::seed_format.
!! For the seed_format "name X" The file would look like:
!! \file{
!! \# I'm a comment
!!    n       2.000000e-01
!!    p       1.000000e-01
!!  c12       2.000000e-01
!!  o16       3.000000e-01
!! ca40       1.500000e-01
!! tc98       0.500000e-01 }
!!
!! Note that empty lines or lines starting with "#" at the beginning of the file are skipped.
!! Another example for seed format "A Z X skip" looks like:
!! \file{
!!    1   0    2.000000e-01  unimportant
!!    1   1    1.000000e-01  unimportant
!!   12   6    2.000000e-01  unimportant
!!   16   8    3.000000e-01  unimportant
!!   40  20    1.500000e-01  unimportant
!!   98  48    0.500000e-01  unimportant }
!!
!! If the mass fractions do not sum up to one, this subroutine will rescale
!! them. In order to make this routine work, also
!! \ref parameter_class::read_initial_composition has to be set to "yes".
!!
!! @see parameter_class::seed_file, parameter_class::seed_format
!! @author Moritz Reichert
!!
!! @remark This subroutine replaced the previous "read_seed_urs" subroutine of
!!         M. Ugliano.
!!
!! \b Edited:
!!         - 01.06.22 - MR
!! .
subroutine read_seed(iniab)
  use parameter_class, only: max_fname_len,seed_format,seed_file
  use global_class,    only: isotope,net_size
  use error_msg_class, only: raise_exception,str_to_float,str_to_int,&
                             write_data_to_std_out
  implicit none
  ! The pass
  real(r_kind),dimension(net_size),intent(out) :: iniab     !< initial abundances
  ! Internal variables
  character(max_fname_len) :: current_col,col_name,col_unit !< Helper variables
  character(max_fname_len) :: help_reader                   !< Helper variable
  integer                  :: col_count                     !< Amount of columns
  integer                  :: i,j                           !< Loop variable
  integer                  :: A_i,Z_i,N_i,X_i,Y_i,Name_i    !< Column indices
  integer                  :: seed_unit                     !< Seed file identifier
  integer                  :: istat                         !< Status variable
  integer                  :: skip_count                    !< Amount of lines skipped
  logical                  :: skipping                      !< helper variable
  integer                  :: l_file                        !< length of the file
  real(r_kind)             :: mafra_norm                    !< Mass fraction norm
  integer,dimension(:),allocatable      :: A,Z,N            !< Mass number, atomic number, neutron number
  real(r_kind),dimension(:),allocatable :: Y,X              !< Abundance, mass fraction
  character*5,dimension(:),allocatable  :: Name             !< Name of nucleus
  character(max_fname_len),dimension(:),allocatable :: row_read  !< temporary array for one line

  INFO_ENTRY("read_seed")

  ! Initialize all index variables
  A_i = -1; Z_i = -1; N_i = -1; X_i = -1; Y_i = -1; Name_i = -1

  ! Convert the custom string to lower case letters
  do i = 1, max_fname_len
    j = iachar(seed_format(i:i))
    if (j>= iachar("A") .and. j<=iachar("Z") ) then
         seed_format(i:i) = achar(iachar(seed_format(i:i))+32)
    else
         continue
    end if
  end do

  ! Get correct index from the columns and count them
  current_col = ''
  col_count = 0
  do i = 2,max_fname_len
     current_col = trim(adjustl(current_col))//seed_format(i-1:i-1)
     if ((seed_format(i-1:i-1) .ne. ' ') .and. (seed_format(i:i) .eq. ' ')) then
        ! Count the columns
        col_count = col_count +1
        ! Get the name of the column
        col_name = ''
        name_loop : do j=1,max_fname_len
           if ((current_col(j:j) .eq. ':' ) .or. (current_col(j:j) .eq. ' ' )) exit name_loop
           col_name = trim(adjustl(col_name))//current_col(j:j)
        end do name_loop
        ! Get the unit of the column
        j=j+1
        col_unit = ''
        unit_loop : do j=j,max_fname_len
           if ((current_col(j:j) .eq. ' ' )) exit unit_loop
           col_unit = trim(adjustl(col_unit))//current_col(j:j)
        end do unit_loop

        ! The column and everything is known here
        ! Store the column index and units
        select case (trim(adjustl(col_name)))
        case('a')
           A_i = col_count !index of the column (mass number)
        case('z')
           Z_i = col_count !index of the column (proton number)
        case('n')
           N_i = col_count !index of the column (neutron number)
        case('x')
           X_i = col_count !index of the column (mass fraction)
        case('y')
           Y_i = col_count !index of the column (abundance)
        case('name')
           Name_i = col_count !index of the column (nucleus name)
        case('skip')
           continue
        case default
           call raise_exception('Problem when analyzing: "'//&
                                trim(adjustl(seed_format))//'". '//NEW_LINE('A')//&
                                'Can not read column "'//trim(adjustl(col_name))//'".'//&
                                NEW_LINE('A')//'Column type unknown.',"read_seeds",&
                                390021)
        end select

        ! Clear the string for the current column
        current_col = ''
      end if
    end do

    ! open the trajectory file
    seed_unit= open_infile(seed_file)
    ! Count how many lines to skip
    l_file = 0
    skip_count = 0
    skipping = .True.
    do ! determine the number of records
       read(seed_unit,"(A)",iostat=istat)help_reader
       if (istat .ne. 0)exit
       help_reader = adjustl(help_reader) ! Move everything to the left

       ! Check if the line is to skip or not
       if ((skipping .eqv. .True.) .and. (len_trim(help_reader) .eq. 0)) then
         skip_count = skip_count+1
         cycle
       elseif ((skipping .eqv. .True.) .and. (help_reader(1:1) .eq. "#")) then
         skip_count = skip_count+1
         cycle
       else
         skipping = .False.
       end if

       ! Check for empty lines in between
       if ((skipping .eqv. .False.) .and. (len_trim(help_reader) .eq. 0)) then
         if (VERBOSE_LEVEL .gt. 1) write(*,*) "Warning, seed file ended unexpectedly."
         exit
       end if

       ! Count the amount of entries
       l_file = l_file + 1
    end do

    ! Allocate arrays
    allocate(row_read(col_count),stat=istat)
    if (istat /= 0) call raise_exception('allocation failed.',"read_seeds",390001)
    allocate(A(l_file),Z(l_file),N(l_file),Y(l_file),X(l_file),Name(l_file),stat=istat)
    if (istat /= 0) call raise_exception('allocation failed.',"read_seeds",390001)

    ! Skip first rows
    rewind(seed_unit)
    do i=1,skip_count
      read(seed_unit,*)
    end do

    ! Initialize
    A(:) = 0; Z(:) = 0; N(:) = 0; Y(:) = 0; X(:) = 0 ; Name(:) = "     "
    ! Finally read the file
    do i=1,l_file
      read(seed_unit,*)row_read
      if (A_i .ne. -1) A(i) = str_to_int(row_read(A_i))
      if (Z_i .ne. -1) Z(i) = str_to_int(row_read(Z_i))
      if (N_i .ne. -1) N(i) = str_to_int(row_read(N_i))
      if (Y_i .ne. -1) Y(i) = str_to_float(row_read(Y_i))
      if (X_i .ne. -1) X(i) = str_to_float(row_read(X_i))
      if (Name_i .ne. -1) Name(i) = trim(adjustl(row_read(Name_i)))
      ! write(*,*)row_read(2)
    end do

    ! Construct other entries
    if ((A_i .eq. -1) .and. (Z_i .ne. -1) .and. (N_i .ne. -1)) then
      A(:) = Z(:)+N(:)
    elseif ((A_i .ne. -1) .and. (Z_i .eq. -1) .and. (N_i .ne. -1)) then
      Z(:) = A(:)-N(:)
    elseif ((A_i .ne. -1) .and. (Z_i .ne. -1) .and. (N_i .eq. -1)) then
      N(:) = A(:)-Z(:)
    elseif ((A_i .ne. -1) .and. (Z_i .ne. -1) .and. (N_i .ne. -1)) then
      continue
    else
      if (Name_i .ne. -1) then
          do j=1,l_file
            iso_loop: do i=1,net_size
            if(trim(adjustl(Name(j))).eq.trim(adjustl(isotope(i)%name))) then
              A(j) = isotope(i)%mass
              Z(j) = isotope(i)%p_nr
              N(j) = isotope(i)%n_nr
              exit iso_loop
            end if
          end do iso_loop
        end do
      else
        call raise_exception("Not possible to read seeds. Columns missing!",&
                             "read_seeds",390022)
      end if
    end if

    ! Store abundances / mass fractions
    if ((Y_i .ne. -1) .and. (X_i .eq. -1)) then
      X(:) = Y(:)*A(:)
    elseif ((Y_i .eq. -1) .and. (X_i .ne. -1)) then
      do i=1,l_file
        ! Avoid division by zero, filter with sunet later
        if (A(i) .eq. 0) A(i) = 1
        Y(i) = X(i)/A(i)
      end do
    elseif ((Y_i .ne. -1) .and. (X_i .ne. -1)) then
      continue
    else
      call raise_exception("Not possible to read seeds. Give either Y or X!",&
                           "read_seeds",390023)
    end if
    ! At this point, Y, X, A, Z, and N are known (or the code crashed)
    ! Renormalize
    mafra_norm = sum(X)
    call write_data_to_std_out("Mass conservation pre norm",num_to_str(mafra_norm))
    X(:) = X(:)/mafra_norm
    Y(:) = Y(:)/mafra_norm
    mafra_norm = sum(X)
    call write_data_to_std_out("Mass conservation post norm",num_to_str(mafra_norm))

    ! Also ensure initial zero Y array
    iniab(:) = 0
    ! Find correct entries
    do j=1, l_file
      inner_loop: do i=1, net_size
        if (Name_i .ne. -1) then
          ! Use the name as it is more unique
          if(trim(adjustl(Name(j))).eq.trim(adjustl(isotope(i)%name))) then
            iniab(i) = Y(j)
            exit inner_loop
          end if
        else
          ! Otherwise use the proton and mass number and hope for the best
          if ((Z(j) .eq. isotope(i)%p_nr) .and. (A(j) .eq. isotope(i)%mass)) then
            iniab(i) = Y(j)
            exit inner_loop
          end if
        end if
      end do inner_loop
    end do

    ! Close the file again
    call close_io_file(seed_unit,seed_file)
    INFO_EXIT("read_seed")

end subroutine read_seed



!>
!! Function to convert from a given time unit to seconds.
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = time_unit_conversion("ms")
!!~~~~~~~~~~~~~~
!! b will be \f$ 10^{-3} \f$ as it is the conversion from ms to s.
!!
!! @note New time units can be defined here.
!!
!! @author Moritz Reichert
function time_unit_conversion(unit)
  use parameter_class, only: trajectory_format
  implicit none
  character(len=*),intent(in)  :: unit                 !< Input unit
  real(r_kind)                 :: time_unit_conversion !< Conversion factor to obtain seconds

  select case (trim(adjustl(unit)))
     case('')
        time_unit_conversion = 1.0
     case('s')
        time_unit_conversion = 1.0
     case('ms')
        time_unit_conversion = 1.0d-3
     case('h')
        time_unit_conversion = 360
     case('min')
        time_unit_conversion = 60
     case('yrs')
        time_unit_conversion = 31536000
     case default
        call raise_exception('Problem when analyzing: "'//&
                             trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//'Time unit "'//&
                             trim(adjustl(unit))//'" unknown.',"time_unit_conversion",&
                             390004)
     end select
end function time_unit_conversion


!>
!! Function to convert from a given temperature unit to GK
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = temp_unit_conversion("K")
!!~~~~~~~~~~~~~~
!! b will be \f$ 10^{-9} \f$ as it is the conversion from K to GK.
!!
!! @note New temperature units can be defined here.
!!
!! @author Moritz Reichert
function temp_unit_conversion(unit)
  use parameter_class, only: trajectory_format
  implicit none
  character(len=*),intent(in) :: unit                 !< Input unit
  real(r_kind)                :: temp_unit_conversion !< Conversion factor to obtain GK

  select case (trim(adjustl(unit)))
     case('')
        temp_unit_conversion = 1.0
     case('-')
        temp_unit_conversion = 1.0
     case('gk')
        temp_unit_conversion = 1.0
     case('k')
        temp_unit_conversion = 1.0d-9
     case('mev')
        temp_unit_conversion = 11.604
     case('ev')
        temp_unit_conversion = 11604.52*1d-9
     case('kev')
        temp_unit_conversion = 11604.52*1d-6
     case default
        call raise_exception('Problem when analyzing: "'//&
                             trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//'Temperature unit "'//&
                             trim(adjustl(unit))//'" unknown.',"temp_unit_conversion",&
                             390005)
   end select
end function temp_unit_conversion


!>
!! Function to convert from a given density unit to g/ccm
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = dens_unit_conversion("g/m^3")
!!~~~~~~~~~~~~~~
!! b will be \f$ 10^{-6} \f$ as it is the conversion from g/m^3 to g/ccm.
!!
!! @note New density units can be defined here.
!!
!! @author Moritz Reichert
function dens_unit_conversion(unit)
  use parameter_class, only: trajectory_format
  implicit none
  character(len=*),intent(in) :: unit                 !< Input unit
  real(r_kind)                :: dens_unit_conversion !< Conversion factor to obtain g/ccm

  select case (trim(adjustl(unit)))
     case('')
        dens_unit_conversion = 1.0
     case('-')
        dens_unit_conversion = 1.0
     case('g/cm^3')
        dens_unit_conversion = 1.0
     case('g/m^3')
        dens_unit_conversion = 1d-6
     case('kg/cm^3')
        dens_unit_conversion = 1d3
     case('kg/m^3')
        dens_unit_conversion = 1d-3
     case default
        call raise_exception('Problem when analyzing: "'//&
                             trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//'Density unit "'//&
                             trim(adjustl(unit))//'" unknown.',"dens_unit_conversion",&
                             390006)
   end select
end function dens_unit_conversion


!>
!! Function to convert from a given distance unit to km
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = dist_unit_conversion("km")
!!~~~~~~~~~~~~~~
!! b will be \f$ 1 \f$ as it is the conversion from km to km.
!!
!! @note New distance units can be defined here.
!!
!! @author Moritz Reichert
function dist_unit_conversion(unit)
  use parameter_class, only: trajectory_format
  implicit none
  character(len=*),intent(in) :: unit                 !< Input unit
  real(r_kind)                :: dist_unit_conversion !< Conversion factor to obtain km

  select case (trim(adjustl(unit)))
     case('')
        dist_unit_conversion = 1.0
     case('-')
        dist_unit_conversion = 1.0
     case('km')
        dist_unit_conversion = 1.0
     case('m')
        dist_unit_conversion = 1e-3
     case('cm')
        dist_unit_conversion = 1e-5
     case default
        call raise_exception('Problem when analyzing: "'//&
                             trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//'Radius unit "'//&
                             trim(adjustl(unit))//'" unknown.',"dist_unit_conversion",&
                             390007)
     end select
end function dist_unit_conversion



!>
!! Function to convert from a given neutrino temperature unit to MeV
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = nutemp_unit_conversion("mev")
!!~~~~~~~~~~~~~~
!! b will be \f$ 1 \f$ as it is the conversion from Mev to Mev.
!!
!! @note New (anti-)neutrino temperature units can be defined here.
!!
!! @author Moritz Reichert
function nutemp_unit_conversion(unit)
  use parameter_class, only: trajectory_format
  implicit none
  character(len=*),intent(in) :: unit                   !> Input unit
  real(r_kind)                :: nutemp_unit_conversion !> Conversion factor to obtain MeV

  select case (trim(adjustl(unit)))
     case('')
        nutemp_unit_conversion = 1.0
     case('-')
        nutemp_unit_conversion = 1.0
     case('mev')
        nutemp_unit_conversion = 1.0
     case default
        call raise_exception('Problem when analyzing: "'//&
                             trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//&
                             'Neutrino temperature unit "'//&
                             trim(adjustl(unit))//'" unknown.',"nutemp_unit_conversion",&
                             390008)
  end select
end function nutemp_unit_conversion


!>
!! Function to convert from a given neutrino energy unit to neutrino temperatures in MeV
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = nuenergy_unit_conversion("mev")
!!~~~~~~~~~~~~~~
!! b will be \f$ 3.15^{-1} \f$ as it is the conversion from neutrino energies to
!! neutrino temperatures (assuming average energies).
!!
!! @note New (anti-)neutrino energy units can be defined here.
!!
!! @author Moritz Reichert
function nuenergy_unit_conversion(unit)
  use parameter_class, only: trajectory_format
  implicit none
  character(len=*),intent(in) :: unit                     !< Input unit
  real(r_kind)                :: nuenergy_unit_conversion !< Conversion factor to obtain temperatures in MeV

  select case (trim(adjustl(unit)))
     case('')
        nuenergy_unit_conversion = 1.0/3.151374374 ! Energy to temperature
     case('-')
        nuenergy_unit_conversion = 1.0/3.151374374
     case('mev')
        nuenergy_unit_conversion = 1.0/3.151374374
     case default
        call raise_exception('Problem when analyzing: "'//&
                             trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//&
                             'Neutrino energy unit "'//&
                             trim(adjustl(unit))//'" unknown.',"nuenergy_unit_conversion",&
                             390009)
 end select
end function nuenergy_unit_conversion


!>
!! Function to convert from a given neutrino luminosity unit to erg/s
!!
!! ### Example
!!~~~~~~~~~~~~~~.f90
!! b = nulumin_unit_conversion("erg/s")
!!~~~~~~~~~~~~~~
!! b will be \f$ 1 \f$ as it is the conversion from erg/s to erg/s.
!!
!! @note New (anti-)neutrino luminosity units can be defined here.
!!
!! @author Moritz Reichert
function nulumin_unit_conversion(unit)
  use parameter_class, only: trajectory_format
  implicit none
  character(len=*),intent(in) :: unit                     !< Input unit
  real(r_kind)                :: nulumin_unit_conversion  !< Conversion factor to obtain erg/s

  select case (trim(adjustl(unit)))
     case('')
        nulumin_unit_conversion = 1.0
     case('-')
        nulumin_unit_conversion = 1.0
     case('erg/s')
        nulumin_unit_conversion = 1.0
     case default
        call raise_exception('Problem when analyzing: "'//&
                             trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//&
                             'Neutrino luminosity unit "'//&
                             trim(adjustl(unit))//'" unknown.',"nulumin_unit_conversion",&
                             390010)
  end select
end function nulumin_unit_conversion


!>
!! Reads the trajectory file.
!!
!! The columns and structure of the file is determined by the user
!! defined parameter \ref parameter_class::trajectory_format.
!! The default trajectory format is given by: <b> "time temp dens rad ye" </b>.
!! This format will assume a file that looks like:
!! \file{
!! \#time[s], T9[GK], density[g/cm^3], R[km], Ye
!! \#--------------------------------------------
!! 0.0000e+00 1.0965e+01 8.7096e+12 4.9273e+01 0.03488
!! 2.5000e-05 9.9312e+00 7.4131e+12 5.0361e+01 0.03488
!! ...}
!! A more complex format, <b> "time:ms x y z log_dens log_temp:K ye" </b>, would assume a file like:
!! \file{
!! \#t [ms]       x[km]       y[km]      z[km]      log(rho[cgs]) log(T[K])   Ye
!! \#----------------------------------------------------------------------------------
!!  0.0000E+00  0.1051E+02  0.3774E+01 -0.4538E+01  0.1364E+02  0.1016E+02  0.1604E-01
!!  0.2521E-01  0.1046E+02  0.4724E+01 -0.4541E+01  0.1363E+02  0.1060E+02  0.1604E-01
!! ...}
!! Neutrino luminosities and energies are also read here. A format could look like:
!! <b> "time temp dens rad ye lnue lanue tnue tanue" </b>. \n
!! Possible column names are:
!! <pre style="line-height:.3">
!! - time   : The time (s) of the trajectory
!! - temp   : The temperature (GK) of the trajectory
!! - dens   : The density (g/ccm) of the trajectory
!! - rad    : The radius (km) of the trajectory
!! - x,y,z  : x, y, and z coordinates (km) of the trajectory
!! - ye     : The electron fraction
!! - lanue  : Anti-electron neutrino luminosities (erg/s)
!! - lnue   : Electron neutrino luminosities (erg/s)
!! - tanue  : Anti-electron neutrino temperatures (MeV)
!! - tnue   : Electron neutrino temperatures (MeV)
!! - eanue  : Anti-electron neutrino energies (MeV)
!! - enue   : Electron neutrino energies (MeV)
!! - lanux  : Anti-neutrino luminosities of heavy neutrinos (erg/s)
!! - lnux   : Neutrino luminosities  of heavy neutrinos (erg/s)
!! - tanux  : Anti-neutrino temperatures  of heavy neutrinos (MeV)
!! - tnux   : Neutrino temperatures of heavy neutrinos (MeV)
!! - eanux  : Anti-neutrino energies of heavy neutrinos (MeV)
!! - enux   : Neutrino energies of heavy neutrinos (MeV)
!! - skip   : skip a column
!! .
!! </pre>
!!
!! @note Lines starting with "#" or blank lines are skipped
!!
!! @author Moritz Reichert
!!
!! \b Created: 17.12.20
subroutine custom_read()
   use parameter_class, only: trajectory_file, trajectory_format, &
                              max_fname_len, read_initial_composition
   use nuflux_class
   use file_handling_class
   use hydro_trajectory
   implicit none
   !
   integer                                  :: traj_unit, i,j, istat
   integer                                  :: col_count
   integer                                  :: time_i,temp_i,dens_i,rad_i,ye_i  !< index of the columns
   integer                                  :: x_i,y_i,z_i                      !< index for x, y, and z if r not given
   real(r_kind)                             :: conv_time,conv_temp,conv_dens    !< unit conversion factors
   real(r_kind)                             :: conv_rad,conv_ye                 !< unit conversion factors
   real(r_kind)                             :: conv_x,conv_y,conv_z             !< unit conversion factors
   real(r_kind)                             :: x_tmp,y_tmp,z_tmp                !< Temporary storage of x,y, and z to convert them
   integer                                  :: nuetemp_i,anuetemp_i             !< index of the neutrino columns
   integer                                  :: Lnue_i,Lanue_i                   !< index of the neutrino columns
   real(r_kind)                             :: conv_Tnue,conv_Tanue             !< unit conversion factors
   real(r_kind)                             :: conv_Lnue,conv_Lanue             !< unit conversion factors
   integer                                  :: nuxtemp_i,anuxtemp_i             !< index of the neutrino columns
   integer                                  :: Lnux_i,Lanux_i                   !< index of the neutrino columns
   real(r_kind)                             :: conv_Tnux,conv_Tanux             !< unit conversion factors
   real(r_kind)                             :: conv_Lnux,conv_Lanux             !< unit conversion factors
   character(max_fname_len)                 :: current_col,col_name,col_unit
   real(r_kind),dimension(:),allocatable    :: row_read                         !< temporary array for one line
   real(r_kind)                             :: tny                              !< tiny value, preventing zeros.
   logical                                  :: log_switch                       !< flag to indicate if one columns is logarithmic
   logical                                  :: ls_time,ls_temp,ls_dens,ls_rad   !< which column is it logarithmic?
   logical                                  :: ls_x,ls_y,ls_z,ls_ye             !< which column is it logarithmic?
   logical                                  :: ls_Tnue,ls_Tanue,ls_Lnue,ls_Lanue!< which column is it logarithmic?
   logical                                  :: ls_Tnux,ls_Tanux,ls_Lnux,ls_Lanux!< which column is it logarithmic?
   integer                                  :: skip_count                       !< Check how many lines to skip
   logical                                  :: skipping                         !< Helper variable
   character(max_fname_len)                 :: help_reader                      !< Helper variable

   INFO_ENTRY("custom_read")

   tny = tiny(tny)

   ! open the trajectory file
   traj_unit= open_infile(trajectory_file)

   ! Convert the custom string to lower case letters
   do i = 1, max_fname_len
     j = iachar(trajectory_format(i:i))
     if (j>= iachar("A") .and. j<=iachar("Z") ) then
          trajectory_format(i:i) = achar(iachar(trajectory_format(i:i))+32)
     else
          continue
     end if
   end do


   ! Count how many lines to skip and how many are there
   zsteps=0
   skip_count = 0
   skipping = .True.
   do ! determine the number of records
      read(traj_unit,"(A)",iostat=istat)help_reader
      ! Go out, file ended
      if (istat .ne. 0)exit
      ! Move string to the left
      help_reader = adjustl(help_reader)
      ! Check if the line is to skip or not
      if ((skipping .eqv. .True.) .and. (len_trim(help_reader) .eq. 0)) then
        skip_count = skip_count+1
        cycle
      elseif ((skipping .eqv. .True.) .and. (help_reader(1:1) .eq. "#")) then
        skip_count = skip_count+1
        cycle
      else
        skipping = .False.
      end if
      zsteps= zsteps + 1
   end do

   ! allocate arrays
   allocate(ztime(zsteps), ztemp(zsteps), zdens(zsteps), &
            zrad(zsteps) , zYe(zsteps)  , stat=istat)
   if (istat /= 0) call raise_exception('allocation failed.',"custom_read",390001)

   ! Allocate neutrino arrays
   if (((nuflag.ge.1)) .and. (neutrino_mode .eq. 'from_file')) then
         allocate(tnue(zsteps) ,tnuebar(zsteps) , &
                  nlume(zsteps),nlumebar(zsteps), &
                  tnux(zsteps) ,tnuxbar(zsteps) , &
                  nlumx(zsteps),nlumxbar(zsteps), &
                  stat=istat)
         if (istat /= 0) call raise_exception('allocation failed.',"custom_read",390001)
   end if

   ! Set default values
   time_i = 0
   temp_i = 0
   dens_i = 0
   rad_i  = 0
   x_i    = 0
   y_i    = 0
   z_i    = 0
   ye_i   = 0
   ! Same for neutrino quantities (electron neutrinos)
   nuetemp_i  = 0; anuetemp_i = 0
   Lnue_i     = 0; Lanue_i    = 0
   ! Same for the heavy neutrinos
   nuxtemp_i  = 0; anuxtemp_i = 0
   Lnux_i     = 0; Lanux_i    = 0

   ! Get correct index from the columns and count them
   current_col = ''
   col_count = 0
   do i = 2,max_fname_len
      current_col = trim(adjustl(current_col))//trajectory_format(i-1:i-1)
      if ((trajectory_format(i-1:i-1) .ne. ' ') .and. (trajectory_format(i:i) .eq. ' ')) then
         ! Count the columns
         col_count = col_count +1
         ! Get the name of the column
         col_name = ''
         name_loop : do j=1,max_fname_len
            if ((current_col(j:j) .eq. ':' ) .or. (current_col(j:j) .eq. ' ' )) exit name_loop
            col_name = trim(adjustl(col_name))//current_col(j:j)
         end do name_loop
         ! Get the unit of the column
         j=j+1
         col_unit = ''
         unit_loop : do j=j,max_fname_len
            if ((current_col(j:j) .eq. ' ' )) exit unit_loop
            col_unit = trim(adjustl(col_unit))//current_col(j:j)
         end do unit_loop

         ! Check if the column is already the logarithm
         select case (col_name(1:4))
         case('log_')
            log_switch = .True.
         case default
            log_switch = .False.
         end select

         ! Remove the "log_" from the string to identify the column
         if (log_switch) then
            col_name = col_name(5:)
         end if

         ! Store the column index and units
         select case (trim(adjustl(col_name)))
         case('time')
            time_i = col_count !index of the column
            ls_time=log_switch !check if it should be logarithmic
            ! Get the conversion factor to seconds
            conv_time = time_unit_conversion(trim(adjustl(col_unit)))
         case('temp')
            temp_i = col_count !index of the column
            ls_temp=log_switch !check if it should be logarithmic
            ! Get the conversion factor to GK
            conv_temp = temp_unit_conversion(trim(adjustl(col_unit)))
         case('dens')
            dens_i = col_count !index of the column
            ls_dens=log_switch !check if it should be logarithmic
            ! Get the conversion factor to g/ccm
            conv_dens = dens_unit_conversion(trim(adjustl(col_unit)))
         case('rad')
            rad_i = col_count !index of the column
            ls_rad=log_switch !check if it should be logarithmic
            ! Get the conversion factor to km
            conv_rad = dist_unit_conversion(trim(adjustl(col_unit)))
         case('x')
            x_i = col_count !index of the column
            ls_x=log_switch !check if it should be logarithmic
            ! Get the conversion factor to km
            conv_x = dist_unit_conversion(trim(adjustl(col_unit)))
         case('y')
            y_i = col_count !index of the column
            ls_y=log_switch !check if it should be logarithmic
            ! Get the conversion factor to km
            conv_y = dist_unit_conversion(trim(adjustl(col_unit)))
         case('z')
            z_i = col_count !index of the column
            ls_z=log_switch !check if it should be logarithmic
            ! Get the conversion factor to km
            conv_z = dist_unit_conversion(trim(adjustl(col_unit)))
         case('ye')
            ye_i  = col_count
            ls_ye = log_switch !check if it should be logarithmic, for ye not very probable
            ! Ye does not have a unit
            conv_ye = 1.0
         case('tnue')
            nuetemp_i = col_count !index of the column
            ls_tnue   = log_switch !check if it should be logarithmic
            ! Get the conversion factor to MeV
            conv_Tnue = nutemp_unit_conversion(trim(adjustl(col_unit)))
         case('tanue')
            anuetemp_i = col_count  !index of the column
            ls_tanue   = log_switch !check if it should be logarithmic
            ! Get the conversion factor to MeV
            conv_Tanue = nutemp_unit_conversion(trim(adjustl(col_unit)))
         case('enue')
            nuetemp_i = col_count
            ls_tnue   = log_switch !check if it should be logarithmic
            ! Get the conversion factor to MeV and neutrino temperatures
            conv_Tnue = nuenergy_unit_conversion(trim(adjustl(col_unit)))
         case('eanue')
            anuetemp_i = col_count
            ls_tanue   = log_switch !check if it should be logarithmic
            ! Get the conversion factor to MeV and neutrino temperatures
            conv_Tanue = nuenergy_unit_conversion(trim(adjustl(col_unit)))
         case('lnue')
            Lnue_i = col_count
            ls_lnue= log_switch !check if it should be logarithmic
            ! Get the conversion factor to erg/s
            conv_Lnue= nulumin_unit_conversion(trim(adjustl(col_unit)))
         case('lanue')
            Lanue_i = col_count
            ls_lanue= log_switch !check if it should be logarithmic
            ! Get the conversion factor to erg/s
            conv_Lanue= nulumin_unit_conversion(trim(adjustl(col_unit)))
         case("tnux")
            nuxtemp_i = col_count !index of the column
            ls_tnux   = log_switch !check if it should be logarithmic
            ! Get the conversion factor to MeV
            conv_Tnux = nutemp_unit_conversion(trim(adjustl(col_unit)))
         case("tanux")
            anuxtemp_i = col_count  !index of the column
            ls_tanux   = log_switch !check if it should be logarithmic
            ! Get the conversion factor to MeV
            conv_Tanux = nutemp_unit_conversion(trim(adjustl(col_unit)))
         case("enux")
            nuxtemp_i = col_count
            ls_tnux   = log_switch !check if it should be logarithmic
            ! Get the conversion factor to MeV and neutrino temperatures
            conv_Tnux = nuenergy_unit_conversion(trim(adjustl(col_unit)))
         case("eanux")
            anuxtemp_i = col_count
            ls_tanux   = log_switch !check if it should be logarithmic
            ! Get the conversion factor to MeV and neutrino temperatures
            conv_Tanux = nuenergy_unit_conversion(trim(adjustl(col_unit)))
         case("lnux")
            Lnux_i = col_count
            ls_lnux= log_switch !check if it should be logarithmic
            ! Get the conversion factor to erg/s
            conv_Lnux= nulumin_unit_conversion(trim(adjustl(col_unit)))
         case("lanux")
            Lanux_i = col_count
            ls_lanux= log_switch !check if it should be logarithmic
            ! Get the conversion factor to erg/s
            conv_Lanux= nulumin_unit_conversion(trim(adjustl(col_unit)))
         case('skip')
            continue
         case default
            call raise_exception('Problem when analyzing: "'//&
                                 trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//&
                                 'Can not read column "'//trim(adjustl(col_name))//'".'//&
                                 NEW_LINE('A')//'Column type unknown.',"custom_read",&
                                 390011)
         end select
         ! Clear the string for the current column
         current_col = ''
      end if
   end do

   ! Check if all necessary hydro-columns were there
   if ((time_i .eq. 0) .or. (temp_i .eq. 0) .or. (dens_i .eq. 0)) then
       call raise_exception('Problem when analyzing: "'//&
                            trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//&
                            'Necessary column is missing.',"custom_read",&
                            390012)
   end if

   ! Ye is only necessary if no seed file is given
   if ((ye_i .eq. 0) .and. (.not. read_initial_composition)) then
       call raise_exception('Problem when analyzing: "'//&
                             trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//&
                             'Ye was missing.',"custom_read",&
                             390012)
   end if

   ! Radius is necessary for neutrinos
   if (((rad_i .eq. 0) .and. (x_i .eq. 0)) .and. (nuflag.ge.1) ) then
       call raise_exception('Problem when analyzing: "'//&
                             trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//&
                             'Radius was missing.',"custom_read",&
                             390012)
   end if

   ! Check if all neutrino-columns were there
   if (((nuflag.ge.1)) .and. (neutrino_mode .eq. 'from_file')) then
      if ((nuetemp_i .eq. 0) .or. (anuetemp_i .eq. 0) .or. (Lnue_i .eq. 0) .or. &
          (Lanue_i .eq. 0))  then
          call raise_exception('Problem when analyzing: "'//&
                              trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//&
                              'Necessary column with neutrino information is missing.',"custom_read",&
                              390013)
      end if
   end if

   ! Check that the columns for the heavies neutrinos are there in case of neutral
   ! current neutrino reactions
   if (((nuflag .eq. 4) .or. (nuflag .eq. 5)) .and. (neutrino_mode .eq. 'from_file')) then
      if ((nuxtemp_i .eq. 0) .or. (anuxtemp_i .eq. 0) .or. (Lnux_i .eq. 0) .or. &
          (Lanux_i .eq. 0))  then
          call raise_exception('Problem when analyzing: "'//&
                              trim(adjustl(trajectory_format))//'". '//NEW_LINE('A')//&
                              'Necessary column with (heavy) neutrino information is missing.',"custom_read",&
                              390013)
      end if
   end if


   ! read the file
   rewind(traj_unit)
   ! Skip the first rows again
   do i=1,skip_count
     read(traj_unit,*)
   end do
   ! allocate array to read one row
   allocate(row_read(col_count),stat=istat)
   if (istat /= 0) call raise_exception('Allocation failed.',"custom_read",390001)

   ! Read the trajectory and convert the units
   do i=1,zsteps
      read(traj_unit,*,iostat=istat)row_read
      ! Check that everything is okay
      if (istat .ne. 0) then
        call raise_exception("Could not read trajectory file."//NEW_LINE('A')//&
                             "Line "//int_to_str(skip_count+i)//" raised a problem.",&
                             "custom_read",390024)
      end if

      ! Time
      if (ls_time)       ztime(i) = 10**(row_read(time_i))
      if (.not. ls_time) ztime(i) = row_read(time_i)
      ztime(i) = ztime(i)*conv_time

      ! Density
      if (ls_dens)       zdens(i) = 10**(row_read(dens_i))
      if (.not. ls_dens) zdens(i) = row_read(dens_i)
      zdens(i) = zdens(i)*conv_dens

      ! Temperature
      if (ls_temp)       ztemp(i) = 10**(row_read(temp_i))
      if (.not. ls_temp) ztemp(i) = row_read(temp_i)
      ztemp(i) = ztemp(i)*conv_temp

      ! Electron fraction
      if (ye_i .ne. 0) then
        if (ls_ye)       zYe(i) = 10**(row_read(ye_i))
        if (.not. ls_ye) zYe(i) = row_read(ye_i)
        zYe(i) = zYe(i)*conv_ye
      else
        zYe(i) = 0.5 ! Fill with dummy value if ye will be taken from a seed file
      end if

      ! Radius is special in the case that x, y, or z was given instead of r
      if (rad_i .ne. 0) then
         if (ls_rad)       zrad(i)  = 10**(row_read(rad_i))
         if (.not. ls_rad) zrad(i)  = row_read(rad_i)
         zrad(i)  = zrad(i)*conv_rad
      elseif ((rad_i .eq. 0) .and. (x_i .ne. 0) .and. (y_i .eq. 0) .and. (z_i .eq. 0)) then
         if (ls_x)       x_tmp = 10**(row_read(x_i))
         if (.not. ls_x) x_tmp = row_read(x_i)
         zrad(i) =  x_tmp*conv_x
      elseif ((rad_i .eq. 0) .and. (x_i .ne. 0) .and. (y_i .ne. 0) .and. (z_i .eq. 0)) then
         ! check if it was converted by the logarithm
         if (ls_x)       x_tmp = 10**(row_read(x_i))
         if (.not. ls_x) x_tmp = row_read(x_i)
         if (ls_y)       y_tmp = 10**(row_read(y_i))
         if (.not. ls_y) y_tmp = row_read(y_i)
         ! Calculate radius
         zrad(i) =  sqrt((x_tmp*conv_x)**2+(y_tmp*conv_y)**2)
      elseif ((rad_i .eq. 0) .and. (x_i .ne. 0) .and. (y_i .ne. 0) .and. (z_i .ne. 0)) then
         ! check if it was converted by the logarithm
         if (ls_x)       x_tmp = 10**(row_read(x_i))
         if (.not. ls_x) x_tmp = row_read(x_i)
         if (ls_y)       y_tmp = 10**(row_read(y_i))
         if (.not. ls_y) y_tmp = row_read(y_i)
         if (ls_z)       z_tmp = 10**(row_read(z_i))
         if (.not. ls_z) z_tmp = row_read(z_i)
         zrad(i) =  sqrt((x_tmp*conv_x)**2+(y_tmp*conv_y)**2+(z_tmp*conv_z)**2)
      else
         ! The radius is not always necessary. E.g., if there is no extrapolation
         ! and no neutrinos.
         zrad(i) = 1
         ! call raise_exception('Not possible to read the radius, specify r or x!.',"custom_read",&
         !                     390014)
      end if
      ! Read neutrino quantities
      if ((nuflag.ge.1) .and. (neutrino_mode .eq. 'from_file')) then

         ! neutrino temperature
         if (ls_tnue)       tnue(i) = 10**row_read(nuetemp_i)
         if (.not. ls_tnue) tnue(i) = row_read(nuetemp_i)
         tnue(i)     = tnue(i)*conv_Tnue

         ! anti-neutrino temperature
         if (ls_tanue)       tnuebar(i) = 10**row_read(anuetemp_i)
         if (.not. ls_tanue) tnuebar(i) = row_read(anuetemp_i)
         tnuebar(i)  = tnuebar(i)*conv_Tanue

         ! neutrino luminosity
         if (ls_lnue)       nlume(i) = 10**row_read(Lnue_i)
         if (.not. ls_lnue) nlume(i) = row_read(Lnue_i)
         nlume(i)    = nlume(i)*conv_Lnue

         if (ls_lanue)       nlumebar(i) = 10**row_read(Lanue_i)
         if (.not. ls_lanue) nlumebar(i) = row_read(Lanue_i)
         nlumebar(i) = nlumebar(i)*conv_Lanue

         ! Avoid small and negative numbers
         if (nlume(i).le.0.d0) nlume(i)       = tny
         if (nlumebar(i).le.0.d0) nlumebar(i) = tny
         if (tnue(i).le.0.d0) tnue(i)         = tny
         if (tnuebar(i).le.0.d0) tnuebar(i)   = tny
      end if
      ! Do the same for the heavy neutrinos
      if ( ((nuflag .eq. 4) .or. (nuflag .eq. 5)) .and. (neutrino_mode .eq. 'from_file')) then
        ! Neutrino temperature of heavy neutrinos
        if (ls_tnux)      tnux(i) = 10**row_read(nuxtemp_i)
        if (.not. ls_tnux) tnux(i) = row_read(nuxtemp_i)
        tnux(i)     = tnux(i)*conv_Tnux

        ! Neutrino temperature of heavy anti neutrinos
        if (ls_tanux)      tnuxbar(i) = 10**row_read(anuxtemp_i)
        if (.not. ls_tanux) tnuxbar(i) = row_read(anuxtemp_i)
        tnuxbar(i)  = tnuxbar(i)*conv_Tanux

        ! Neutrino luminosity of heavy neutrinos
        if (ls_lnux)      nlumx(i) = 10**row_read(Lnux_i)
        if (.not. ls_lnux) nlumx(i) = row_read(Lnux_i)
        nlumx(i)    = nlumx(i)*conv_Lnux

        ! Neutrino luminosity of heavy anti neutrinos
        if (ls_lanux)      nlumxbar(i) = 10**row_read(Lanux_i)
        if (.not. ls_lanux) nlumxbar(i) = row_read(Lanux_i)
        nlumxbar(i) = nlumxbar(i)*conv_Lanux

        ! Avoid small and negative numbers
        if (nlumx(i).le.0.d0) nlumx(i)       = tny
        if (tnux(i).le.0.d0) tnux(i)         = tny
        if (tnuxbar(i).le.0.d0) tnuxbar(i)   = tny
        if (nlumxbar(i).le.0.d0) nlumxbar(i) = tny
      end if
   end do

   ! log representation
   ztemp = dlog(ztemp)
   zdens = dlog(zdens)
   zrad  = dlog(zrad)

   ! Same for neutrino quantities if necessary
   if ((nuflag.ge.1) .and. (neutrino_mode .eq. 'from_file')) then
     nlume    = dlog(nlume)
     nlumebar = dlog(nlumebar)
     tnue     = dlog(tnue)
     tnuebar  = dlog(tnuebar)
   end if
   ! Same for heavy neutrinos if necessary
   if ( ((nuflag .eq. 4) .or. (nuflag .eq. 5)) .and. &
        (neutrino_mode .eq. 'from_file')) then
     nlumx    = dlog(nlumx)
     nlumxbar = dlog(nlumxbar)
     tnux     = dlog(tnux)
     tnuxbar  = dlog(tnuxbar)
   end if

   ! close file
   call close_io_file(traj_unit,trajectory_file)

   INFO_EXIT("custom_read")

end subroutine custom_read


!>
!! Make tests to ensure the trajectory is okay.
!!
!! Otherwise complain. Tests are, for example, the occurence of
!! negative radii at some point or check if the electron fraction
!! is out of bounds (\f$ 0 \le y_e \le 1 \f$).
!!
!! @author Moritz Reichert
subroutine consistency_check()
   use hydro_trajectory
   implicit none
   integer :: i !< Loop variable

   ! Bunch of consistency checks for "from_file" mode
   if (trim(trajectory_mode).EQ.'from_file') then


     ! Negative time
     if (ztime(zsteps) .lt. 0) then
        call raise_exception("Last point of trajectory had negative time ("//&
                              trim(adjustl((num_to_str(ztime(zsteps)))))//&
                              ")!"//NEW_LINE('A')//'Please check your input trajectory.',&
                              "consistency_check",&
                              390015)
      ! Electron fraction out of bounds
      elseif ((maxval(zYe) .gt. 1) .or. (minval(zYe) .lt. 0)) then
          call raise_exception("The electron fraction was out of bounds (0 < Ye < 1)! "//&
                               NEW_LINE('A')//'Please check your input trajectory.',&
                               "consistency_check",&
                               390016)
      end if


      ! Check negative entries
      do i=1, zsteps
         ! Negative radius
         if (isnan(dexp(zrad(i)))) then
            call raise_exception('Problem with input radius: "'//&
                                 'The radius was below zero."'//NEW_LINE('A')//&
                                 'Please check your input trajectory.',"consistency_check",&
                                 390017)
         ! Negative temperature
         elseif (isnan(dexp(ztemp(i)))) then
            call raise_exception('Problem with input temperature: "'//&
                                 'The temperature was below zero".'//NEW_LINE('A')//&
                                 'Please check your input trajectory.',"consistency_check",&
                                 390018)
         ! Negative density
         elseif (isnan(dexp(zdens(i)))) then
            call raise_exception('Problem with input density: "'//&
                                 'The density was below zero".'//NEW_LINE('A')//&
                                 'Please check your input trajectory.',"consistency_check",&
                                 390019)
         end if
      end do

      ! Check for increasing time
      do i=1, zsteps-1
         if (ztime(i) .ge. ztime(i+1)) call raise_exception("Time in the trajectory was not monotonically increasing."//&
                                                            NEW_LINE('A')//"Check your trajectory!","consistency_check",&
                                                            390020)
      end do
   end if


end subroutine consistency_check



!>
!! Readini cleanup subroutine
!!
!! At the moment, the subroutine does nothing.
!!
!! @author Moritz Reichert
!!
!! \b Edited:
!!       - MR 28.02.2017
!!       - OK 18.06.2017
!! .
subroutine readini_finalize()
   use file_handling_class
   implicit none
   ! not much here really
end subroutine readini_finalize



end module readini
