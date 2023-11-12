!> @file driver.f90
!!
!! The error file code for this file is ***W14***.
!! @brief Contains the main \ref driver routine
!!
#include "macros.h"
program driver

  use global_class,        only: net_size
  use winnse_module,       only: pf
  use error_msg_class,     only: raise_exception, data_creation_mode, &
                                 write_final_stats_rate_creation
  use nucstuff_class,      only: el_ab
  use file_handling_class, only: delete_io_file
  use gear_module,         only: init_gear_solver
  use parameter_class
  use expansion_module
  use mergesort_module
  use nuclear_heating
  use readini
  use single_zone_vars
  use hdf5_module
  use network_init_module
  use timestep_module
  use analysis
  implicit none

  integer                    :: cnt             !< Iteration count
  character(max_fname_len)   :: parameter_file  !< Path to parameter file
  character(max_fname_len)   :: binary_out_dir  !< Path to binary dir to safe the rates
  logical                    :: init            !< Flag to indicate whether the
                                                !  time step has to be set to the
                                                !  initial value.

!---
! Parameter setup
!---
  if (iargc().eq.1) then
    call getarg(1, parameter_file)
  elseif (iargc().eq.2) then
    call getarg(1, parameter_file)
    call getarg(2, binary_out_dir)
    ! If the binary_out_dir is given,
    ! the code is in rate creation mode
    ! and the print statements will be different
    data_creation_mode = .true.
  else
    print '(A)', "ERROR: wrong number of command-line arguments"
    print '(A)', "Syntax for network calculation: ./winnet <parameter-file>"
    print '(A)', "Syntax for network data creation: ./winnet <parameter-file> <binary-out-dir>"
    stop
  endif


  call set_default_param()
  call read_param (trim(parameter_file))
  call output_param()

!---
! Initializations
!---
  if (iargc() .eq. 2) then
    ! Dont use a prepared network when trying to prepare one
    use_prepared_network = .false.
    call create_rate_folder(trim(binary_out_dir))
    ! Output final stats
    call write_final_stats_rate_creation()
    ! Dont run the network
    stop
  else
    init=.true.       !! local variable, TODO: get rid of it
    call prepare_simulation()
  end if



! initial output
  call output_initial_step(time,T9,rhob,ent,P,Rkm,Y,pf)

!
! --- MAIN EVOLUTION LOOP
!

  cnt = 0
  evolution_loop: do while (evolution_mode.ne.EM_TERMINATE)

     !-- Select timestep
     call select_optimistic_timestep(init)

     !-- Compute the next step
     !   switch between different solvers
     select case (solver)
     case (0) ! implicit Euler
        call advance_implicit_euler(cnt)
     case(1)  ! Gear's method
        call advance_gear(cnt)
     case default
        call raise_exception("Invalid solver. Choose either solver = 0 or 1.","driver",&
                             140003)
     end select

     !-- Switch between modes
     call switch_evolution_mode(init)

     !-- Rotate timelevels: save the 'previous step' quantities
     call el_ab (Y, Ye)
     dYdt(1:net_size)= (Y(1:net_size)-Y_p(1:net_size))/stepsize

     !-- Periodic output
     call output_iteration(cnt+1,time,T9,rhob,ent,P,Rkm,Y,pf)

     ! rotate timelevels
     Y_p(1:net_size)=   Y(1:net_size)
     Ye_p=   Ye
     ent_p=  ent
     P_p=  P
     T9_p=   T9
     T9h_p=  T9h
     cnt=    cnt+1
     time_p= time


  end do evolution_loop

  ! output final step
  call output_final_step(cnt,time,T9,rhob,ent,P,Rkm,Y,pf)

  ! decay to stability
  call el_ab(Y, Ye)
  !  call decay_init(Ye,T9,rhob)
  !  write (*,'(A)',advance='no') 'Decaying to stable'
  !  do i=1,40
  !     call decay(Y)
  !     call printsnap (t,T9,rhob,Y)
  !     write (*,'(A)',advance='no') '.'
  !  end do

  ! output final abundances
#ifndef USE_HDF5
  call finab(Y)
#else
  if (.not. h_finab) then
     call finab(Y)
  else
     call write_finab(Y)
  endif
#endif


  ! finalize modules
  call readini_finalize()
  call expansion_finalize()
  call analysis_finalize()
  call mergesort_finalize()
#ifdef USE_HDF5
   call hdf5_module_finalize()
#endif
  ! stop timing
  call system_clock(cl_end)

  ! output final statistics
  call output_final_stats(cnt,time)

  ! Keep the folder clean with as few files as possible.
  ! If the code reaches this line here, no error has occured.
  ! Therefore delete the error file.
  call delete_io_file("ERR")
end program driver
