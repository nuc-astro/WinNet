!> @file network_init_module.f90
!!
!! The error file code for this file is ***W30***.
!!
!! @brief \ref network_init_module contains \ref network_init and
!!         supplementary functions to set up and initialize the network
!! @author Christian Winteler
!!
!! \b Created: 23.02.09
!! \b Edited:
!!       - 15.06.17
!!       - MR 11.01.20  - introduced error_msg_class
!!       - MR 24.01.21  - splitted things to new rate modules
!! .


!> Module to prepare the reaction network
!!
!! This module contains subroutines that call other modules initialization.
!! In addition it reads reaction rates and many other things.
#include "macros.h"
module network_init_module
  use error_msg_class, only: raise_exception,int_to_str,&
                             num_to_str,write_data_to_std_out
  implicit none

  !
  ! Public and private fields and methods of the module
  !
  public:: &
      prepare_simulation, switch_evolution_mode
  private:: &
      network_init, print_reactions, print_debug, getnames

contains


!>
!! Main initialising subroutine. calls reading subroutines and fills the
!! reactionrate array.
!!
!! Reads included nuclei, Reaclib, theoretical weak rates,
!! neutrino rates, tabulated rates, fission rates (and more)
!! and prepares in this way the reaction rates
!! \ref global_class::rrate.
!!
!! @returns \ref global_class::rrate
!!
!! \b Edited:
!!      - 14.04.15
!!      - OK 15.06.17, made into a module
!! .
subroutine network_init()

  use reaclib_rate_module,    only: init_reaclib_rates, merge_reaclib_rates
  use alpha_decay_rate_module,only: init_alpha_decay_rates, merge_alpha_decays
  use tw_rate_module,         only: init_theoretical_weak_rates, merge_theoretical_weak_rates,&
                                    output_binary_weak_reaction_data
  use fission_rate_module,    only: init_fission_rates, merge_fission_rates
  use detailed_balance,       only: merge_inverse_rates, init_inverse_rates
  use parameter_class,        only: net_source
  use nuflux_class,           only: init_nuflux, merge_neutrino_rates
  use beta_decay_rate_module, only: init_ext_beta_rates,merge_beta_decays
  use benam_class,            only: get_minmax, get_nuclear_properties
  use tabulated_rate_module,  only: init_tabulated_rates, merge_tabulated_rates
  use global_class,           only: rrate, isotope, nreac, net_size
  use format_class,           only: load_format
  use nucstuff_class,         only: init_nucstuff
  use screening_module,       only: init_screening

  !
  integer                    :: alloc_stat !< Allocation state
!===========================procedure division==================================
  INFO_ENTRY("network_init")

  !-- loads custom formats my_format(:)
  call load_format()
  !-- returns nuclide data in isotope(:)
  call get_nuclear_properties()
  !-- returns Amin and Amax for each isotopic chain in minmax
  call get_minmax()
  !-- Initialize nucstuff
  call init_nucstuff()

  !-- Count the total amount of rates
  nreac = 0

  !-- Count and read reaclib rates
  call init_reaclib_rates()
  !-- Count and read alpha decay rates
  call init_alpha_decay_rates()
  !-- Count and read tabulated rates
  call init_tabulated_rates()
  !-- Count and read neutrino rates
  call init_nuflux()
  !-- Count and read external beta decays
  call init_ext_beta_rates()
  !-- Count and read fission rates
  call init_fission_rates()
  !-- Count and read theoretical weak rates
  call init_theoretical_weak_rates()
  !-- Initialize inverse reactions
  call init_inverse_rates()


  ! MR: Note that the following order is important.
  ! For example, TWR have to be merged before neutrino rates!
  ! Otherwise the results change. (Tested with Example_CCSN_wind_bliss)
  ! Furthermore, inverse rates have to be merged at the very end as
  ! they rely on already computed forward rates in the calculation
  ! of the Jacobian.
  !
  !-- Merge reaclib rates into rrate
  call merge_reaclib_rates(rrate,nreac)
  !-- Merge alpha decays
  call merge_alpha_decays(rrate,nreac)
  !-- Merge external beta decays
  call merge_beta_decays(rrate,nreac)
  !-- Merge tabulated rates into rrate
  call merge_tabulated_rates(rrate,nreac)
  !-- Merge theoretical weak rates into rrate
  call merge_theoretical_weak_rates(rrate,nreac)
  !-- Merge neutrino reactions into rrate
  call merge_neutrino_rates(rrate,nreac)
  !-- Merge fission rates into rrate
  call merge_fission_rates(rrate,nreac)
  !-- Merge inverse rates (detailed balance) into rrate
  call merge_inverse_rates(rrate,nreac)

  ! Tell how many rates there are
  if (VERBOSE_LEVEL.ge.1) then
     call write_data_to_std_out("Total reactions after merge",int_to_str(nreac))
  endif

  ! Initialize screening
  call init_screening(nreac)

  if (VERBOSE_LEVEL.gt.2) then
     !----- print list of reactions - debug_all-reactions.txt and debug_weak-reactions.txt
     call print_reactions()
     !----- print reactions for debugging - debug-react.dat
     call print_debug()
  end if


  INFO_EXIT("network_init")

end subroutine network_init


!>
!! All the necessary initializations and settings before the main loop
!!
!! To initialize it calls the init subroutines from other modules,
!! e.g., \ref readini::readini_init, \ref analysis::analysis_init.
!! It also prepares the first time depending on \ref parameter_class::initemp
!! and it fills the abundance array \ref single_zone_vars::y and
!! \ref single_zone_vars::y_p.
!!
!! \b Edited:
!!      - OK 15.06.2017
!!      - OK 27.08.2017
!! .
subroutine prepare_simulation
use parameter_class, only: initial_stepsize, flow_every, T9_analytic,&
                        Rho_analytic, Rkm_analytic, Ye_analytic,&
                        read_initial_composition, nsetemp_cold,&
                        solver, unit_define,&
                        h_flow_every,initemp_hot,initemp_cold,&
                        t_analytic, heating_mode
use global_class,    only: net_size, ineu, ipro
use flow_module,     only: flow_init
use winnse_module,   only: winnse_guess, nse_init
use inter_module,    only: inverse_interp1d, interp1d
use thermal_neutrino_module, only: thermal_neutrino_init
use hydro_trajectory
use analysis
use file_handling_class
use mergesort_module
use nucstuff_class
use pardiso_class
use ls_timmes_eos_module
use expansion_module
use nuclear_heating
use single_zone_vars
use readini
use parser_module
use gear_module, only: init_gear_solver
implicit none
integer                               :: i
real(r_kind)                          :: t_i
real(r_kind)                          :: T9_i
real(r_kind)                          :: rho_i
real(r_kind)                          :: Rkm_i
real(r_kind)                          :: Ye_i
integer                               :: alloc_stat
real(r_kind)                          :: dummy_hot, dummy_cold,dummy_temp !< Helper variables
                                                                          !< to get temperature starting point
integer                               :: dummy_index_cold,dummy_index_hot !< Helper variables
                                                                          !< to get temperature starting point
! Variables used for the Expansion/EOS----
type(timmes_eos_state)                :: inistate
integer                               :: eos_status
logical                               :: init,converged

  ! define units
    call unit_define()

  ! initializations
    evolution_mode= EM_INIT
    net_size= 0
    stepsize= initial_stepsize
  ! initialise Clenshaw-Curtis integration method
    !call init_cc(ncc)

  ! initialise mergesort module
    call mergesort_init()

  ! initialise network, returns net_size (and network data in modules)
    call network_init()

  ! initialise readini module (for opening debug files)
    call readini_init()

  ! initialise analysis module
    call analysis_init()


  ! sparse prepares the usage of Compressed Sparse Column Format
    call sparse()

  ! initialise flow subroutine, has to get called after sparse!
    if (flow_every.gt.0) then
      call flow_init()
#ifdef USE_HDF5
    elseif (h_flow_every .gt. 0) then
      call flow_init()
#endif
    end if

  ! allocate arrays
    allocate(Y_p(net_size),stat=alloc_stat)
    if ( alloc_stat /= 0) call raise_exception('Allocating "Y_p" failed.',&
                                               "prepare_simulation",300001)
    allocate(Y(net_size),stat=alloc_stat)
    if ( alloc_stat /= 0) call raise_exception('Allocating "Y" failed.',&
                                               "prepare_simulation",300001)
    allocate(dYdt(net_size),stat=alloc_stat)
    if ( alloc_stat /= 0) call raise_exception('Allocating "dYdt" failed.',&
                                               "prepare_simulation",300001)
    allocate(f(net_size),stat=alloc_stat)
    if ( alloc_stat /= 0) call raise_exception('Allocating "f" failed.',&
                                               "prepare_simulation",300001)

  ! read nuclear data needed for nse subroutine
    call nse_init()


    if (trim(trajectory_mode).EQ.'from_file') then
       !-- 'initemp_hot' parameter controls the starting point of the trajectory:
       !      0: start from the beginning;
       !     >0: start from the last T peak, at which the T exceeded initemp

       !-- 'initemp_cold' controls the minimum allowed temperature in case the
       !   trajectory is too cold. If this is the case, the trajectory will start at the
       !   peak temperature that obeys T>initemp_cold and T<initemp_hot.
       init_index= 1
       dummy_index_cold= 1
       dummy_index_hot= 1
       dummy_temp= -99
       if (initemp_hot.gt.0d0) then
          dummy_hot = log(initemp_hot)
          dummy_cold= log(initemp_cold)

          do i=1,zsteps
             ! Find peak temperature
             if ((ztemp(i) .gt. dummy_temp) .and. (ztemp(i) .gt. dummy_cold)) then
                dummy_temp= ztemp(i)
                dummy_index_cold= i
             end if

             if (ztemp(i).gt.dummy_hot) dummy_index_hot= i
          end do

          ! If there is a high temperature (initemp_hot) included, use the init index
          ! at this point in temperature, otherwise use maximum of trajectory
          if ((dummy_index_hot .eq. 1) .and. (ztemp(1) .lt. dummy_hot)) then
             init_index = dummy_index_cold
          else
             init_index = dummy_index_hot
          end if


          if(init_index.eq.zsteps) then
             call raise_exception('Temperature of the trajectory is too high (>initemp_hot).'//NEW_LINE("A")//&
                                  'The last temperature of the trajectory is '&
                                  //trim(adjustl(num_to_str(dexp(ztemp(zsteps)))))//" GK."//&
                                  NEW_LINE("A")//"Initemp is "&
                                  //trim(adjustl(num_to_str(initemp_hot)))//" GK. "//&
                                  "Check the trajectory and units.",&
                                  "prepare_simulation",300003)
          endif
       endif
       t_i=   ztime(init_index)
       T9_i=  dexp(ztemp(init_index))
       rho_i= dexp(zdens(init_index))
       Rkm_i= dexp(zrad(init_index))
       Ye_i=  zye(init_index)

       ! Check if the temperature is very far from what is expected
       ! and if initemp_hot is the temperature to start with.
       if (T9_i .gt. initemp_hot) then
         ! Is the temperature of the grid point of the trajectory
         ! within 0.0001 of the desired temperature?
         if (abs(T9_i-initemp_hot) .gt. 1d-9) then
           ! Get the correct time
           call inverse_interp1d(zsteps,ztemp,initemp_hot*(1+1d-9),ztime,t_i)
           ! Get the correct temperature, density, radius, ye
           call interp1d (zsteps,ztime,t_i,ztemp,T9_i)
           call interp1d (zsteps,ztime,t_i,zdens,rho_i)
           call interp1d (zsteps,ztime,t_i,zrad,Rkm_i,itype=itype_LINEAR)
           call interp1d (zsteps,ztime,t_i,zye,Ye_i,.True.,itype=itype_LINEAR)
         end if
       end if

       if (VERBOSE_LEVEL.ge.1) then
          call write_data_to_std_out("Thermodynamic trajectory","from file")
          call write_data_to_std_out("Starting index",int_to_str(init_index))
       end if
    else if (trim(trajectory_mode).EQ.'analytic') then

       ! Find the time of T=T_init
       call find_start_value(T9_analytic,initemp_hot,t_analytic,converged,t_i)

       ! For, e.g., constant temperature the above NR won't work
       if (converged) then
          ! take the later value in case it is too hot at the initial one
          t_i = max(t_i,t_analytic)
       else
          ! if the convergence failed, just start at the initial time
          t_i=   t_analytic
       end if

       ! get the density, temperature and radius at the initial time
       T9_i = parse_string(T9_analytic,t_i)
       rho_i= parse_string(rho_analytic,t_i)
       Rkm_i= parse_string(Rkm_analytic,t_i)
       Ye_i = parse_string(Ye_analytic,t_i)

       if (VERBOSE_LEVEL.ge.1) then
          call write_data_to_std_out("Thermodynamic trajectory","analytic")
       end if
    else
      call raise_exception("Unknown trajectory_mode '"//trim(adjustl(trajectory_mode))//"'.",&
                           "prepare_simulation",300004)
    end if

  ! initialise expansion module (for opening debug files)
    call expansion_init()

  ! start timing
    call system_clock(cl_start,cl_rate,cl_cmax)

    time=     t_i
    time_p=   t_i
    T9=       T9_i
    T9_p=     T9_i
    rhob=     rho_i
    rhob_p=   rho_i
    Rkm=      Rkm_i
    Rkm_p=    Rkm_i
    Ye=       Ye_i
    Ye_p=     Ye_i

  ! check whether we start in NSE conditions; if not read seeds from progenitor
    if (read_initial_composition) then
       if (VERBOSE_LEVEL.ge.1) call write_data_to_std_out("Initial abundance source","File")
       call read_seed(Y)
       ! Recalculate the Ye
       Ye  = sum(isotope(1:net_size)%p_nr*Y(1:net_size))
       Ye_p= Ye
       evolution_mode= EM_NETHOT
       if (T9.gt.nsetemp_cold) evolution_mode= EM_NSE
    else
       if (VERBOSE_LEVEL.ge.1) call write_data_to_std_out("Initial abundance source","NSE")
       Y(1:net_size)= 0.d0
       Y(ineu)= 1.d0-Ye  ! abundance of neutrons
       Y(ipro)= Ye       ! abundance of protons
       if (T9.gt.nsetemp_cold) then
          evolution_mode= EM_NSE
          if (VERBOSE_LEVEL.ge.1) then
             call write_data_to_std_out("Evolution mode","NSE")
          endif
          call winnse_guess(T9, rhob, Ye, Y(ineu), Y(ipro), Y)
          do i=1,net_size
             if(Y(i).lt.1.d-25) Y(i)=0.d0
          enddo
       endif
    endif
    Y_p(1:net_size)= Y(1:net_size)

  ! initialize evolution mode
    init= .true.
    call switch_evolution_mode (init)

  ! get initial entropy -> needed for expansion with timmes_eos
    inistate%abar = sum(Y(1:net_size)*isotope(1:net_size)%mass) &
                  / sum(Y(1:net_size))
    call timmes_eos(ink,T9*1.d9,rhob,Ye,inistate,eos_status)
    if(eos_status.ne.0) call raise_exception("Error when calling the EOS. Called with:"//&
                                             NEW_LINE("A")//&
                                             "T   : "//num_to_str(T9)//" GK"//NEW_LINE("A")//&
                                             "Rho : "//num_to_str(rhob)//" g/ccm"//NEW_LINE("A")//&
                                             "Ye  : "//num_to_str(Ye)//NEW_LINE("A")//&
                                             "Abar: "//num_to_str(inistate%abar),&
                                             "prepare_simulation",300005)

    ent = inistate%s
    ent_p= ent
    P = inistate%p
    P_P= P
    if (VERBOSE_LEVEL.ge.1)call write_data_to_std_out("Initial entropy",num_to_str(ent),"[kB/nuc]")
  ! initialize nuclear heating module
    if (heating_mode .ne. 0) then
       call nuclear_heating_init (T9, rhob, Ye, Y, ent)
       ! Initialize thermal neutrinos as well
       call thermal_neutrino_init()
    endif

  ! initialise Gear's method if used
  if(solver==1) call init_gear_solver(Y,t_i,initial_stepsize)


end subroutine prepare_simulation


!>
!! Create a folder with binary files that contain all necessary reaction data
!!
!! This routine creates binary files in the exact format that WinNet needs to
!! initialize the rate arrays. The reading of these files can be done much faster
!! then reading the original files. The files are stored in the folder specified
!! by the path argument. The folder is created if it does not exist.
!! Using these files is specifically useful when running many trajectories.
!!
!! @author M. Reichert
!! @date 21.07.23
subroutine create_rate_folder(path)
    use global_class,        only: net_size
    use tw_rate_module,      only: output_binary_weak_reaction_data
    use reaclib_rate_module, only: output_binary_reaclib_reaction_data
    use parameter_class,     only: unit_define, output_binary_parameter_data, &
                                   max_fname_len, output_param_prepared_network
    use mergesort_module,    only: mergesort_init
    use fission_rate_module, only: output_binary_fission_reaction_data
    use nuflux_class,        only: output_binary_neutrino_reaction_data
    use benam_class,         only: output_binary_network_data
    implicit none
    character(len=*), intent(in) :: path
    character(max_fname_len)     :: path_dir

    ! define units
    call unit_define()

    net_size= 0

    ! initialise mergesort module
    call mergesort_init()
    ! initialise network, returns net_size (and network data in modules)
    call network_init()
    ! Make sure the path ends with a slash
    path_dir = trim(path)//"/"
    ! Create folder
    call system('mkdir -p '//trim(adjustl(path_dir)))

    ! Output weak rates
    call output_binary_parameter_data(trim(adjustl(path_dir)))
    call output_binary_network_data(trim(adjustl(path_dir)))
    call output_binary_weak_reaction_data(trim(adjustl(path_dir)))
    call output_binary_reaclib_reaction_data(trim(adjustl(path_dir)))
    call output_binary_fission_reaction_data(trim(adjustl(path_dir)))
    call output_binary_neutrino_reaction_data(trim(adjustl(path_dir)))

    ! Give information
    call output_param_prepared_network(trim(adjustl(path_dir)))

end subroutine create_rate_folder



!>
!! Switches between NSE and network modes (EM_NSE / EM_NETHOT / EM_NETCOLD)
!!
!! The following conditions are included to switch the evolution mode:
!!
!! <table>
!! <caption id="multi_row">Evolution mode switches</caption>
!! <tr><th> Evolution mode <th> Condition
!! <tr><td> EM_NSE         <td> \f$T_9\f$ > parameter_class::nsetemp_hot
!! <tr><td> EM_NETHOT      <td> parameter_class::temp_reload_exp_weak_rates < \f$T_9\f$ < parameter_class::nsetemp_cold
!! <tr><td> EM_NETCOLD     <td> \f$T_9\f$ < parameter_class::temp_reload_exp_weak_rates
!! <tr><td> EM_TERMINATE   <td> Depends on termination criterion, see parameter_class::termination_criterion
!! </table>
!!
!! The following figure illustrates the different regions for an example trajectory:
!! @image html temperature_regimes.png width=400px
!!
!! @note This subroutine also changes the timestep at the transition from NSE
!!       to the network. For the gear solver, gear_module::set_timestep
!!       has to be called.
!!
!! @todo At this point it could be useful to check here the difference
!! between the saved value of dYdt from the network and the dYdt
!! calculated by the nse abchange function, right after the network->nse
!! switch
!!
!! @see single_zone_vars::evolution_mode
!!
!! \b Edited:
!!      - 17.12.16
!! .
subroutine switch_evolution_mode (init)
use single_zone_vars, only: T9,rhob,Ye,Y,Y_p,time,evolution_mode
use parameter_class,  only: nsetemp_hot,nsetemp_cold,nse_descend_t9start,&
    iwformat,termination_criterion,final_time,final_dens,final_temp, &
    temp_reload_exp_weak_rates, solver, initial_stepsize, heating_density,&
    heating_mode
use global_class,   only: ihe4,ineu,ipro,net_size, heating_switch
use tw_rate_module, only: weak, reload_exp_weak_rates
use gear_module,    only: init_gear_solver
use winnse_module,  only: winnse_descend, nse_init
use hydro_trajectory
implicit none
logical, intent(inout):: init !< Flag that indicates a necessary initialization

     if (evolution_mode.eq.EM_INIT) then
     ! initalize at the very beginning
        if (T9.gt.nsetemp_hot) then
           evolution_mode= EM_NSE
           call winnse_descend(nse_descend_t9start, T9, rhob, Ye,&
                               Y(ineu), Y(ipro), Y)
           call write_data_to_std_out("Starting evolution mode","NSE")
           ! print *, "Evolution: starting in NSE mode"
        else
           evolution_mode= EM_NETHOT
           call write_data_to_std_out("Starting evolution mode","Network")
           ! print *, "Evolution: starting in network mode"
        endif
     elseif (T9.gt.nsetemp_hot .and. evolution_mode.ne.EM_NSE) then
     ! network --> NSE
        evolution_mode= EM_NSE
        print *, "Evolution mode: switching to nse"
        call winnse_descend(nse_descend_t9start, T9, rhob, Ye,&
                            Y(ineu), Y(ipro), Y)
        init = .true. ! reset to correctly compute the timestep below
        if(solver==1) call init_gear_solver(Y,time,initial_stepsize)
     elseif (T9.lt.nsetemp_cold .and. evolution_mode.eq.EM_NSE) then
     ! NSE --> network
        evolution_mode= EM_NETHOT
        call winnse_descend(nse_descend_t9start, T9, rhob, Ye, &
                            Y(ineu), Y(ipro), Y)
        ! winnse_descend modifies abundances Y;
        ! need to copy abundances to previous timelevel, because otherwise
        ! they will be overwritten in the beginning of adapt_stepsize in
        ! advance_next_step() subroutine;
        Y_p(1:net_size)= Y(1:net_size)
        ! initialise Gear's method if used
        init = .true. ! reset to correctly compute the timestep below
        if(solver==1) call init_gear_solver(Y,time,initial_stepsize)
        print *,"Evolution mode: switching to network"
     end if

     ! Check if heating has to be enabled
     if ((.not. heating_switch) .and. (rhob .le. heating_density) .and. &
         (heating_mode .ne. 0)) then
        heating_switch = .true.
        ! Pretty output
        if (.not. init) then
            print *,"Evolution mode: switching on heating"
        end if
     end if

     !-- Problem: beta decay reaction rates are not reliable at T9 < 1.d-2
     if (weak.and.(evolution_mode.eq.EM_NETHOT).and. &
                  (T9.le.temp_reload_exp_weak_rates)) then
        ! try 1: replace theoretical weak rates with experimental ones
        call reload_exp_weak_rates() !TODO: fix and expand this subroutine
        iwformat = 0
        weak = .false.
        evolution_mode = EM_NETCOLD
        print *,"Evolution mode: switching to cold network (only experimental rates)"
     end if

     !-- Determine if simulation should terminate
     select case (termination_criterion)
     case(0) ! when trajectory ends
        if (time.ge.ztime(zsteps)) evolution_mode= EM_TERMINATE
     case(1) ! after a given final time
        if (time.ge.final_time)    evolution_mode= EM_TERMINATE
     case(2) ! when a lower temperature is hit
        if (T9.le.final_temp)      evolution_mode= EM_TERMINATE
     case(3) ! below a certain density threshold
        if (rhob.le.final_dens)    evolution_mode= EM_TERMINATE
     case(4) ! after freeze-out
        if (Y_p(ineu)/sum(Y_p(ihe4+1:net_size)).le.1.0d0) then
           evolution_mode= EM_TERMINATE
        endif
     case default
        call raise_exception('Invalid termination_criterion ('//trim(adjustl(int_to_str(termination_criterion)))//&
                             ').',"prepare_simulation", 300006)
     end select

end subroutine switch_evolution_mode








!>
!! Print an overview of the reactions involved in the calculation.
!!
!! An example file (_debug_all-reactions.txt_ ) looks like:
!! \file{...
!! O N E P A R T I C L E R E A C T I O N S
!! 1         n --->    p                              Q=  0.782 MEV     ffn
!! Reaction type: weak
!! 2         n --->    p                              Q=  0.000 MEV     nen
!! Reaction type: neutrino
!! 3         p --->    n                              Q= -0.782 MEV     ffn
!! Reaction type: weak
!! ...}
!!
!! Furthermore, a file ( _debug_weak-reactions.txt_ ) containing all
!! weak reactions is created that looks like:
!! \file{
!! OVERVIEW OF ALL WEAK-REACTIONS IN THE NETWORK
!! 1         n --->    p                              Q=  0.782 MEV     ffn
!! 2         n --->    p                              Q=  0.000 MEV     nen
!! 3         p --->    n                              Q= -0.782 MEV     ffn
!! 4         p --->    n                              Q=  0.000 MEV    nebp }
!!
!! @note Fission reactions do not obey the default reaclib format and are not
!!       visualized correct.
!!
!! \b Edited:
!!     - MR: 21.01.21 - included interface (get_net_name) to net_names instead of accessing it directly
subroutine print_reactions()

  use global_class, only: nreac, rrate, net_names
  use benam_class, only: get_net_name
  use file_handling_class
  implicit none

  integer             :: reactions,weak_unit
  integer             :: i,i1,i2,i3,cnt,wcnt
  !character(31)       :: rateSortKey

  reactions= open_outfile('debug_all-reactions.txt')
  weak_unit= open_outfile('debug_weak-reactions.txt')
  ! write(reactions,'(i1,4/)')1
  write(reactions,'(t11,a41)')  &
       'OVERVIEW OF ALL REACTIONS IN THE NETWORK'
  ! write(weak_unit,'(i1,4/)')1
  write(weak_unit,'(t11,a46)')  &
       'OVERVIEW OF ALL WEAK-REACTIONS IN THE NETWORK'
  i1=1
  i2=1
  i3=1
  cnt=1
  wcnt=1
  do i=1,nreac

     select case(rrate(i)%group)
     case(1)
        if(i1.eq.1) then
           write(reactions,'(5/,t11,a)')&
                'O N E P A R T I C L E R E A C T I O N S'
           i1=0
           cnt=1
           wcnt=1
        end if
        write(reactions,100)cnt,net_names(rrate(i)%parts(1)),  &
             get_net_name(rrate(i)%parts(2)),  &
             rrate(i)%q_value,rrate(i)%source
        cnt = cnt+1
        if(rrate(i)%is_weak) then
           write(weak_unit,100)wcnt,net_names(rrate(i)%parts(1)),  &
                get_net_name(rrate(i)%parts(2)),  &
                rrate(i)%q_value,rrate(i)%source
           wcnt = wcnt+1
        end if
     case(2)
        write(reactions,200)cnt,net_names(rrate(i)%parts(1)),  &
             get_net_name(rrate(i)%parts(2)),  &
             get_net_name(rrate(i)%parts(3)),  &
             rrate(i)%q_value,rrate(i)%source
        cnt = cnt+1
        if(rrate(i)%is_weak) then
           write(weak_unit,200)wcnt,get_net_name(rrate(i)%parts(1)),  &
                get_net_name(rrate(i)%parts(2)),  &
                get_net_name(rrate(i)%parts(3)),  &
                rrate(i)%q_value,rrate(i)%source
           wcnt = wcnt+1
        end if
     case(3)
        if (rrate(i)%parts(5).eq.0) then
           write(reactions,300)cnt,get_net_name(rrate(i)%parts(1)),  &
                get_net_name(rrate(i)%parts(2)),  &
                get_net_name(rrate(i)%parts(3)),  &
                get_net_name(rrate(i)%parts(4)),  &
                rrate(i)%q_value,rrate(i)%source
        else
           write(reactions,310)cnt,get_net_name(rrate(i)%parts(1)),  &
                get_net_name(rrate(i)%parts(2)),  &
                get_net_name(rrate(i)%parts(3)),  &
                get_net_name(rrate(i)%parts(4)),  &
                get_net_name(rrate(i)%parts(5)),  &
                rrate(i)%q_value,rrate(i)%source
        end if
        cnt = cnt+1
        if(rrate(i)%is_weak) then
           if (rrate(i)%parts(5).eq.0) then
              write(weak_unit,300)wcnt,get_net_name(rrate(i)%parts(1)),  &
                   get_net_name(rrate(i)%parts(2)),  &
                   get_net_name(rrate(i)%parts(3)),  &
                   get_net_name(rrate(i)%parts(4)),  &
                   rrate(i)%q_value,rrate(i)%source
           else
              write(weak_unit,310)cnt,get_net_name(rrate(i)%parts(1)),  &
                   get_net_name(rrate(i)%parts(2)),  &
                   get_net_name(rrate(i)%parts(3)),  &
                   get_net_name(rrate(i)%parts(4)),  &
                   get_net_name(rrate(i)%parts(5)),  &
                   rrate(i)%q_value,rrate(i)%source
           end if
           wcnt = wcnt+1
        end if
     case(4)
        if(i2.eq.1) then
           write(reactions,'(5/,t11,a)') &
                'T W O P A R T I C L E R E A C T I O N S'
           i2=0
           cnt=1
        end if
        write(reactions,400)cnt,get_net_name(rrate(i)%parts(1)),  &
             get_net_name(rrate(i)%parts(2)),  &
             get_net_name(rrate(i)%parts(3)),  &
             rrate(i)%q_value,rrate(i)%source
        !write(reactions,*) rateSortKey(rrate(i)%group,rrate(i)%parts)
        cnt = cnt+1
     case(5)
        write(reactions,500)cnt,get_net_name(rrate(i)%parts(1)),  &
             get_net_name(rrate(i)%parts(2)),  &
             get_net_name(rrate(i)%parts(3)),  &
             get_net_name(rrate(i)%parts(4)),  &
             rrate(i)%q_value,rrate(i)%source
        cnt = cnt+1
     case(6)
        write(reactions,600)cnt,get_net_name(rrate(i)%parts(1)),  &
             get_net_name(rrate(i)%parts(2)),  &
             get_net_name(rrate(i)%parts(3)),  &
             get_net_name(rrate(i)%parts(4)),  &
             get_net_name(rrate(i)%parts(5)),  &
             rrate(i)%q_value,rrate(i)%source
        cnt = cnt+1
     case(7)
        write(reactions,700)cnt,get_net_name(rrate(i)%parts(1)),  &
             get_net_name(rrate(i)%parts(2)),  &
             get_net_name(rrate(i)%parts(3)),  &
             get_net_name(rrate(i)%parts(4)),  &
             get_net_name(rrate(i)%parts(5)),  &
             get_net_name(rrate(i)%parts(6)),  &
             rrate(i)%q_value,rrate(i)%source
        cnt = cnt+1
     case(8)
        if(i3.eq.1) then
           write(reactions,'(5/,t11,a)')  &
                'T H R E E P A R T I C L E R E A C T I O N S'
           i3=0
           cnt=1
        end if
        if (rrate(i)%parts(6).ne.0) then
           write(reactions,830)cnt,get_net_name(rrate(i)%parts(1)),  &
                get_net_name(rrate(i)%parts(2)),  &
                get_net_name(rrate(i)%parts(3)),  &
                get_net_name(rrate(i)%parts(4)),  &
                get_net_name(rrate(i)%parts(5)),  &
                get_net_name(rrate(i)%parts(6)),  &
                rrate(i)%q_value,rrate(i)%source
        else if(rrate(i)%parts(5).ne.0) then
           write(reactions,820)cnt,get_net_name(rrate(i)%parts(1)),  &
                get_net_name(rrate(i)%parts(2)),  &
                get_net_name(rrate(i)%parts(3)),  &
                get_net_name(rrate(i)%parts(4)),  &
                get_net_name(rrate(i)%parts(5)),  &
                rrate(i)%q_value,rrate(i)%source
        else
           write(reactions,810)cnt,get_net_name(rrate(i)%parts(1)),  &
                get_net_name(rrate(i)%parts(2)),  &
                get_net_name(rrate(i)%parts(3)),  &
                get_net_name(rrate(i)%parts(4)),  &
                rrate(i)%q_value,rrate(i)%source
        end if
        cnt = cnt+1
     end select

     ! Write the type of reaction to the debug file
     if (     rrate(i)%reac_type  == rrt_sf) then
        write(reactions,840) "spontaneous fission"
     else if ( rrate(i)%reac_type == rrt_bf ) then
        write(reactions,840) "beta-delayed fission"
     else if ( rrate(i)%reac_type == rrt_nf  ) then
        write(reactions,840) "neutron-induced fission"
     else if ( rrate(i)%reac_type == rrt_ng  ) then
        write(reactions,840) "n-gamma"
     else if ( rrate(i)%reac_type == rrt_ag  ) then
        write(reactions,840) "alpha-gamma"
     else if ( rrate(i)%reac_type == rrt_pg  ) then
        write(reactions,840) "p-gamma"
     else if ( rrate(i)%reac_type == rrt_gn  ) then
        write(reactions,840) "gamma-n"
     else if ( rrate(i)%reac_type == rrt_ga  ) then
        write(reactions,840) "gamma-alpha"
     else if ( rrate(i)%reac_type == rrt_gp  ) then
        write(reactions,840) "gamma-p"
     else if ( rrate(i)%reac_type == rrt_pn  ) then
        write(reactions,840) "p-n"
     else if ( rrate(i)%reac_type == rrt_np  ) then
        write(reactions,840) "n-p"
     else if ( rrate(i)%reac_type == rrt_an  ) then
        write(reactions,840) "alpha-n"
     else if ( rrate(i)%reac_type == rrt_na  ) then
        write(reactions,840) "n-alpha"
     else if ( rrate(i)%reac_type == rrt_ap  ) then
        write(reactions,840) "alpha-p"
     else if ( rrate(i)%reac_type == rrt_pa  ) then
        write(reactions,840) "p-alpha"
     else if ( rrate(i)%reac_type == rrt_nu  ) then
        write(reactions,840) "neutrino"
     else if ( rrate(i)%reac_type == rrt_o ) then
        write(reactions,840) "other"
     else if ( rrate(i)%reac_type == rrt_betm  ) then
        write(reactions,840) "beta minus"
     else if ( rrate(i)%reac_type == rrt_betp  ) then
        write(reactions,840) "beta plus"
     else if ( rrate(i)%reac_type == rrt_alpd  ) then
        write(reactions,840) "alpha decay"
     else if ( rrate(i)%reac_type == rrt_pemi  ) then
        write(reactions,840) "p emission"
     else if ( rrate(i)%reac_type == rrt_nemi  ) then
        write(reactions,840) "n emission"
     else if ( rrate(i)%reac_type == rrt_ec  ) then
        write(reactions,840) "Electron capture"
     else
        write(reactions,840) "unknown"
     end if


  end do
  call close_io_file(reactions,'debug_all-reactions.txt')
  call close_io_file(weak_unit,'debug_weak-reactions.txt')
  return

100 format (i6,t12,a5,t18,'--->',t22,a5,t57,'Q=',t59,f7.3,t67,'MEV',t74,a4)
200 format (i6,t12,a5,t18,'--->',t22,a5,t27,'+',t28,a5,t57,  &
         'Q=',t59,f7.3,t67,'MEV',t74,a4)
300 format (i6,t12,a5,t18,'--->',t22,a5,t27,'+',t28,a5,t33,'+',t34,a5,t57,  &
         'Q=',t59,f7.3,t67,'MEV',t74,a4)
310 format (i6,t12,a5,t18,'--->',t22,a5,t27,'+',t28,a5,t33,'+',t34,a5,  &
         t39,'+',t40,a5,t57,'Q=',t59,f7.3,t67,'MEV',t74,a4)
400 format (i6,t12,a5,t17,'+',t18,a5,t24,'--->',t28,a5,t57,  &
         'Q=',t59,f7.3,t67,'MEV',t74,a4)
500 format (i6,t12,a5,t17,'+',t18,a5,t24,'--->',t28,a5,t33,'+',t34,a5,t57,  &
         'Q=',t59,f7.3,t67,'MEV',t74,a4)
600 format (i6,t12,a5,t17,'+',t18,a5,t24,'--->',t28,a5,t33,'+',t34,a5,  &
         t39,'+',t40,a5,t57,'Q=',t59,f7.3,t67,'MEV',t74,a4)
700 format (i6,t12,a5,t17,'+',t18,a5,t24,'--->',t28,a5,t33,'+',t34,a5,  &
         t39,'+',t40,a5,t45,'+',t46,a5,t57,'Q=',t59,f7.3,t67,'MEV',t74,a4)
810 format (i6,t12,a5,t17,'+',t18,a5,t23,'+',t24,a5,t30,'--->',t34,a5,  &
         t57,'Q=',t59,f7.3,t67,'MEV',t74,a4)
820 format (i6,t12,a5,t17,'+',t18,a5,t23,'+',t24,a5,t30,'--->',t34,a5,  &
         t39,'+',t40,a,t57,'Q=',t59,f7.3,t67,'MEV',t74,a4)
830 format (i6,t12,a5,t17,'+',t18,a5,t23,'+',t24,a5,t30,'--->',t34,a5,  &
         t39,'+',t40,a,t45,'+',t46,a,t57,'Q=',t59,f7.3,t67,'MEV',t74,a4)
840 format ('Reaction type: ',A)

end subroutine print_reactions





!> Print a debug file of the reaction
!!
!! In contrast to \ref print_reactions, this subroutine prints additional flags
!! of a reaction. The structure of a file looks like:
!!
!! <table>
!! <caption id="multi_row">File structure</caption>
!! <tr><td> Reaclib chapter                <td> isotope names <td> Q-Value <td> Source label
!! <tr><td> Amount destroyed/synthesized   <td> double counting factor
!! <tr><td>"is_weak" <td> "is_resonant" <td> "is_reverse"
!! </table>
!! An example of the file (_debug_reac.dat_) looks like:
!! \file{...
!! 5     n na32    p ne32            -16.606  rath
!! -1 -1  1  1  0  0      1.000
!!     F    F    T
!! 5     p na32    n mg32             18.317  rath
!! -1 -1  1  1  0  0      1.000
!!     F    F    F
!! ...}
!!
subroutine print_debug()

  use global_class, only: nreac, rrate
  use file_handling_class
  implicit none

  integer       :: ofile
  integer       :: i
  character*5,dimension(6)    :: rnam

  ofile= open_outfile('debug_reac.dat')
  do i =1,nreac
     call getnames(i,rnam)
     write(ofile,111)rrate(i)%group,rnam(1:6),rrate(i)%q_value,rrate(i)%source
     write(ofile,222)rrate(i)%ch_amount(1:6),rrate(i)%one_over_n_fac
     write(ofile,333)rrate(i)%is_weak,rrate(i)%is_resonant,rrate(i)%is_reverse
  end do
  call close_io_file(ofile,'debug_reac.dat')
  return

111 format(i1,t3,6a5,t35,f7.3,t44,a4)
222 format(6(i2,1x),t22,f7.3)
333 format(3l5)

end subroutine print_debug


!> Get the names of a certain reaction
!!
!! @returns all nuclei names that participate in a reaction with index ind
subroutine getnames(ind,names)
  use global_class, only: net_names, rrate
  implicit none

  integer, intent(in)                   :: ind   !< Index of the reaction
  character*5,dimension(6),intent(out)  :: names !< Names of the participating isotopes
  integer                               :: i     !< Loop variable

  do i =1,6
     if(rrate(ind)%parts(i).lt.1) then
        names(i) = '     '
     else
        names(i) = net_names(rrate(ind)%parts(i))
     end if
  end do
  return

end subroutine getnames


end module network_init_module
!----------------------------------------------------------------------------
