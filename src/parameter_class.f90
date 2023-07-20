!> @file parameter_class.f90
!!
!! The error file code for this file is ***W34***.
!!
!! @brief Module \ref parameter_class with all simulation parameters
!!

!> Contains all runtime parameters as well as phys and math constants
!!
!! Contains the following definitions:
!! \li the 'kind' attribute for all real types used throught the code;
!! \li all runtime simulation parameters;
!! \li mathematical and physical constants.
!!
!! @author Darko Mocelj
!! @date   22.01.03
!!
!! \b Editors:
!!   \li 01.12.11: Christian Winteler
!!   \li 22.07.13: Marcella Ugliano
!!   \li 11.01.14: Oleg Korobkin
!!   \li 11.01.21: Moritz Reichert
!! .
#include "macros.h"
module parameter_class
  use error_msg_class, only: raise_exception,num_to_str,int_to_str
  implicit none

  public :: unit_define

  !> constants and unit conversion factors
  type,public ::  unit_type
     real(r_kind) :: pi          !< \f$\pi\f$
     real(r_kind) :: mass_n      !< Neutron mass [MeV/c^2]
     real(r_kind) :: mass_p      !< Proton mass [MeV/c^2]
     real(r_kind) :: mass_u      !< Atomic mass unit [MeV/c^2]
     real(r_kind) :: mass_e      !< Electron mass [MeV/c^2]
     real(r_kind) :: me          !< Electron mass [g]
     real(r_kind) :: n_a         !< Avogadro constant
     real(r_kind) :: hbc         !< hbar*c [Mev*fm]
     real(r_kind) :: k_b         !< Boltzman constant [J/K]
     real(r_kind) :: k_mev       !< Boltzmann constant [MeV/K]
     real(r_kind) :: kerg        !< Boltzmann constant [erg/K]
     real(r_kind) :: conv_ev     !< conversion factor eV to J
     real(r_kind) :: clight      !< speed of light in vacuum [cm/s]
     real(r_kind) :: h           !< Planck constant [J*s]
     real(r_kind) :: h_mevs      !< Planck constant [MeV*s]
     real(r_kind) :: hbar_mev    !< Planck constant reduced [MeV*s]
     real(r_kind) :: hix         !< Atomic mass unit [MeV/(cm/s)^2]
     real(r_kind) :: amu         !< Atomic mass unit [g]
     real(r_kind) :: grav        !< Gravitational Constant [cm^3/(g*s^2)]
     real(r_kind) :: msol        !< solar mass [g]
     real(r_kind) :: ergtomev    !< conversion from erg to MeV
  end type unit_type
  type(unit_type),public :: unit !< constants and unit conversion factors

!>-- hardcoded parameters
  integer,public,parameter  :: param_name_len = 200   !< maximum length of a parameter name
  integer,public,parameter  :: max_fname_len = 400    !< maximum length of filenames

!>-- runtime parameters set when reading parameter file
  integer                 :: out_every                      !< main output frequency
  integer                 :: snapshot_every                 !< snapshot output frequency
  integer                 :: h_snapshot_every               !< snapshot output frequency in hdf5 format
  integer                 :: flow_every                     !< flow output frequency
  integer                 :: h_flow_every                   !< flow output frequency in hdf5 format
  integer                 :: timescales_every               !< timescales output frequency
  integer                 :: h_timescales_every             !< timescales output frequency in hdf5 format
  integer                 :: engen_every                    !< Energy generation output frequency
  integer                 :: h_engen_every                  !< Energy generation output frequency in hdf5 format
  integer                 :: nrdiag_every                   !< frequency of NR loop diagnostic messages (wrt iteration counter)
  integer                 :: mainout_every                  !< frequency of mainout output
  integer                 :: h_mainout_every                !< HDF5 output frequency of the mainout
  integer                 :: track_nuclei_every             !< frequency of track nuclei output
  integer                 :: h_track_nuclei_every           !< frequency of track nuclei output in hdf5 format
  integer                 :: top_engen_every                !< frequency of energy generators toplist
  logical                 :: use_alpha_decay_file           !< Switch for using additional alpha decay rates
  integer                 :: alpha_decay_zmin               !< Minimum Z for additional alpha decay rates
  integer                 :: alpha_decay_zmax               !< Maximum Z for additional alpha decay rates
  character(max_fname_len):: alpha_decay_src_ignore         !< Source flag(s) to ignore within the alpha decay rates
  logical                 :: alpha_decay_ignore_all         !< Flag whether rates should actually get replaced or only added
  integer                 :: iwformat                       !< defines format of the weak rates (0 = tabulated, 1 = log<ft>)
  integer                 :: iwinterp                       !< defines the interpolation for the weak rates (0 = bilinear, 1 = bicubic)
  real(r_kind)            :: temp_reload_exp_weak_rates     !< temperature below which one should not use theoretical weak rates so they are replaced with exp. weak rates (min 1.d-2)
  real(r_kind)            :: freeze_rate_temp               !< Tmperature at which rates get frozen (for reacl. rates this should be 1d-2GK)
  integer                 :: nuflag                         !< defines type of neutrino reactions used
  integer                 :: fissflag                       !< defines type of fission fragment distribution used
  integer                 :: expansiontype                  !< defines prescription used for parametrized expansion after the last timestep of the hydro input
  integer                 :: extrapolation_width            !< how many points from the end of trajectory to use when computing residual expansion
  integer                 :: nse_calc_every                 !< Compute NSE abundances every x step.
  character*20            :: trajectory_mode                !< determines how trajectory is calculated (from_file|analytic|expand)
  logical                 :: read_initial_composition       !< specify whether initial distribution of abundances should be read from file
  logical                 :: use_tabulated_rates            !< switch for using tabulated rates (e.g. talysNGrates.dat)
  logical                 :: use_beta_decay_file            !< switch for using different format for beta decays
  character(max_fname_len):: beta_decay_src_ignore          !< Source flag(s) to ignore within the beta decay file
  logical                 :: use_timmes_mue                 !< Use electron chemical potentials from timmes EOS for theoretical weak rates
  logical                 :: use_detailed_balance           !< Calculate the inverse reactions via detailed balance rather than using them form file
  logical                 :: use_detailed_balance_q_reac    !< Use Q-value from reaclib for the calculation of detailed balance
  logical                 :: use_thermal_nu_loss            !< Whether to include thermal neutrino loss or not.
  integer                 :: nu_loss_every                  !< Output neutrino loss and gain.
  integer                 :: h_nu_loss_every                !< Output neutrino loss and gain in hdf5 format.
  character(max_fname_len):: detailed_balance_src_ignore    !< Source flag(s) to ignore within calculated detailed balance
  character(max_fname_len):: detailed_balance_src_q_reac    !< Source flag(s) to use q-value from rate file for inverse reaction
  character(max_fname_len):: detailed_balance_src_q_winvn   !< Source flag(s) to use q-value from winvn file for inverse reaction
  logical                 :: use_neutrino_loss_file         !< Use a file with Qnu values?
  character(max_fname_len):: neutrino_loss_file             !< Path to a file containing Qnu values
  logical                 :: custom_snapshots               !< If true, a file must be provided with numbers in days. Snapshots will be created for these points in time
  logical                 :: h_custom_snapshots             !< Same, but in hdf5 format
  real(r_kind)            :: engen                          !< total energy generation [MeV/s]
  real(r_kind)            :: engen_alpha                    !< energy generation from alpha-decays [MeV/s]
  real(r_kind)            :: engen_beta                     !< energy generation from beta-decays [MeV/s]
  real(r_kind)            :: engen_fiss                     !< energy generation from fission [MeV/s]
  real(r_kind)            :: initemp_cold                   !< T [GK] lowest allowed temperature to start the calculation from
  real(r_kind)            :: initemp_hot                    !< T [GK] for the starting point of the trajectory: =0: from the beginning; >0: from the last T>initemp
  real(r_kind)            :: nse_delt_t9min                 !< Minimum temperature [GK] when descending to desired temperature in NSE
  real(r_kind)            :: nse_nr_tol                     !< Tolerance for the NR loop in the NSE calculation
  integer                 :: nse_max_it                     !< Maximum amount of NSE iterations
  integer                 :: nse_solver                     !< Solver for calculating NSE. 0: Newton-Raphson, 1: Powell's hybrid method
  real(r_kind)            :: nsetemp_cold                   !< T [GK] for the nse->network switch
  real(r_kind)            :: nsetemp_hot                    !< T [GK] for the nse<-network switch
  logical                 :: calc_nsep_energy               !< calculate neutron separation energy?
  logical                 :: h_engen_detailed               !< Output the energy per parent nucleus and reaction type
  real(r_kind)            :: heating_frac                   !< use this fraction of nuclear-generated energy for heating
  real(r_kind)            :: heating_density                !< Density at which nuclear heating will be switched on (-1) to always include heating
  real(r_kind)            :: nse_descend_t9start            !< high initial temperature in GK for winnse_descend subroutine
  real(r_kind)            :: t_analytic                     !< for parameteric trajectories: initial time
  integer                 :: termination_criterion          !< condition to terminate the simulation ([0]=trajectory_file, 1=final_time, 2=final_temp, 3=final_dens, 4=neutron freeze-out)
  real(r_kind)            :: initial_stepsize               !< this value is used as a stepsize at initial step
  real(r_kind)            :: final_time                     !< termination time in seconds
  real(r_kind)            :: final_temp                     !< termination temperature [GK]
  real(r_kind)            :: final_dens                     !< termination density [g/cm3]
  real(r_kind)            :: timestep_max                   !< Maximum factor for the change of the timestep. The new timestep is only allowed to be timestep_max * old_timestep. Default value is 2.
  real(r_kind)            :: timestep_factor                !< Factor for the change of the timestep (see nu in Winteler 2012 Eq. 2.49). Default value is 1.0d-1
  real(r_kind)            :: timestep_hydro_factor          !< Factor for the maximum change of the hydrodynamic quantities (density and temperature)
  real(r_kind)            :: timestep_Ymin                  !< Lower limit of the abundance to contribute to the timestep calculation, default value is 1.0d-10
  real(r_kind)            :: gear_eps                       !< Abundance accuracy for gear solver
  real(r_kind)            :: gear_escale                    !< Normalization cutoff for gear solver, similar to timestep_Ymin for Euler
  real(r_kind)            :: gear_cFactor                   !< Conservative timestep factor for gear solver [0.1, ... , 0.4]
  real(r_kind)            :: gear_nr_eps                    !< Convergence criterion for the newton-raphson of the gear solver
  integer                 :: gear_nr_maxcount               !< Maximum newton-raphson iterations for gear solver
  integer                 :: gear_nr_mincount               !< Minimum newton-raphson iterations for gear solver
  logical                 :: gear_ignore_adapt_stepsize     !< Flag whether gear should ignore the adapt stepsize loop
  logical                 :: timestep_traj_limit            !< Should the timestep be limited by the timestep of the trajectory
  logical                 :: use_htpf                       !< Use high temperature partition functions or not
  logical                 :: h_finab                        !< Store the finab in hdf5 format rather than in ascii format
  integer                 :: solver                         !< solver flag (0 - implicit Euler, 1 - Gear's method, ...), is integer as it is faster than comparing strings every timestep
  integer                 :: heating_mode                   !< Mode for heating: 0 - no heating, 1 - heating using an entropy equation, 2 - heating from the energy generation and trajectory changes
  integer                 :: screening_mode                 !< Mode for coulomb corrections: 0 - no screening, 1 - screening using the prescription of Kravchuk & Yakovlev 2014
  integer                 :: interp_mode                    !< Mode for interpolation of temperature and density
  character(max_fname_len):: trajectory_file                !< name of trajectory data file
  character(max_fname_len):: seed_file                      !< name of file with initial seeds for trajectory run
  character(max_fname_len):: seed_format                    !< Seed format
  character(max_fname_len):: net_source                     !< list of isotopes included in the network (sunet)
  character(max_fname_len):: isotopes_file                  !< properties of all isotopes in the network: masses, partition functions etc. (winvn)
  character(max_fname_len):: htpf_file                      !< high-temperature partition functions (htpf.dat)
  character(max_fname_len):: reaclib_file                   !< reaction rate library (reaclib)
  character(max_fname_len):: fission_rates                  !< reaction library for fission rates (same format as reaclib_file)
  character(max_fname_len):: weak_rates_file                !< weak rates library (twr.dat)
  character(max_fname_len):: tabulated_rates_file           !< tabulated rates library (e.g. talysNGrates.dat)
  character(max_fname_len):: chem_pot_file                  !< tabulated chemical potential of electron gas
  character(max_fname_len):: nsep_energies_file             !< neutron separation energies
  character(max_fname_len):: nunucleo_rates_file            !< neutrino reaction rates on nucleons
  character(max_fname_len):: nuchannel_file                 !< Contains neutrino channel information as in Sieverding et al. 2018
  character(max_fname_len):: nurates_file                   !< Neutrino reactions on heavy nuclei as in Sieverding et al. 2018
  character(max_fname_len):: snapshot_file                  !< File that contains days, where a snapshot should be made
  character(max_fname_len):: nfission_file                  !< Fission table for neutron-induced fission
  character(max_fname_len):: bfission_file                  !< Fission table for beta-delayed and spontaneous fission
  character(max_fname_len):: track_nuclei_file              !< File of nuclei to track. Gives an output similar to mainout.dat
  character(max_fname_len):: alpha_decay_file               !< File with additional alpha decays
  character(max_fname_len):: beta_decay_file                !< File for reading in beta decays in different format
  character(max_fname_len):: trajectory_format              !< Format to read the trajectory
  character(max_fname_len):: neutrino_mode                  !< Similar to trajectory mode
  character(max_fname_len):: T9_analytic                    !< analytic temperature [T9]
  character(max_fname_len):: rho_analytic                   !< analytic density [g/cm3]
  character(max_fname_len):: Rkm_analytic                   !< analytic radial scale [km]
  character(max_fname_len):: Ye_analytic                    !< analytic electron fraction
  character(max_fname_len):: Le                             !< electron-neutrino luminosities [erg/s]
  character(max_fname_len):: Lebar                          !< electron-antineutrino luminosities [erg/s]
  character(max_fname_len):: Enue                           !< average electron-neutrino energies [MeV]
  character(max_fname_len):: Enuebar                        !< average electron-antineutrino energies [MeV]
  character(max_fname_len):: Lx                             !< Muon and Tauon neutrino luminosities [erg/s]
  character(max_fname_len):: Lxbar                          !< Muon and Tauon  antineutrino luminosities [erg/s]
  character(max_fname_len):: Enux                           !< average Muon and Tauon neutrino energies [MeV]
  character(max_fname_len):: Enuxbar                        !< average Muon and Tauon antineutrino energies [MeV]

!>-- Newton-Raphson iterative loop parameters
  integer                  :: nr_maxcount             !< no more that this many iterations in NR
  integer                  :: nr_mincount             !< Minimum iterations in NR
  real(r_kind)             :: nr_tol                  !< exit NR if tolerance less than this value
  integer                  :: adapt_stepsize_maxcount !< max. iterations in adapting the stepsize

!>-- parameters for efficient numerical integration of effphase in the interval [1,infinity]
  real(r_kind), dimension(:), allocatable :: weights  !< weights for the numerical integration
  real(r_kind), dimension(:), allocatable :: xnodes   !< corresponding nodes "
  logical :: weightsCalculated = .false.              !< switch to calculated weights and nodes for [1,infinity]

  integer :: ncc = 256                                  !< nr of points for Clenshaw-Curtis integration
  real(r_kind), dimension(:), allocatable :: dcc,MATcc  !< matrices for Clenshaw-Curtis integration
  real(r_kind), dimension(:,:), allocatable :: Mcc

contains

!>
!! Declares values for the elements in unit_type
!!
subroutine unit_define()
  unit%pi      = 3.14159265358979323846d0
  unit%mass_n  = 939.56533d0     ! mev/c^2
  unit%mass_p  = 938.271998d0    ! mev/c^2
  unit%mass_u  = 931.494013d0    ! mev/c^2
  unit%mass_e  = 0.510998910d0   ! mev
  unit%me      = 9.1093826d-28   ! g
  unit%n_a     = 6.02214179d23   ! mol^(-1)
  unit%hbc     = 197.327053d0    ! mev*fm
  unit%k_b     = 1.380662d-23    ! j/k
  unit%k_mev   = 8.617343d-11    ! mev/k
  unit%conv_ev = 1.602189246d-19 ! 1ev = xxx j
  unit%kerg    = 1.3806504d-16   ! erg
  unit%clight  = 2.99792458d10   ! cm/s
  unit%h       = 6.62606876d-34  ! Js
  unit%hbar_mev= 6.582122d-22    ! MeVs
  unit%h_mevs  = 4.135667273d-21 ! MeVs
  unit%hix     = 1.036427d-18    ! amu in MeV/(cm/s)^2  (convert from [erg/g] to [MeV/baryon])
  unit%amu     = 1.66053873d-24  ! g
  unit%grav    = 6.67384d-8      ! cm^3/(g*s^2)
  unit%msol    = 1.9891d33       ! g
  unit%ergtomev= 0.62415d6       ! 1 erg = xxx MeV
end subroutine unit_define


!>
!! This function reads the parameter file
!!
!! Blank lines as well as comment lines are skipped.
!! The parameters are read into the proper variable
!! that are declarated at the beginning of the file.
!!
!! \b Edited:
!!   - 01.02.14
!! .
subroutine read_param(parfile)
   use file_handling_class
   implicit none

   character(*),intent(in)   :: parfile
   !
   character(1000)           :: line
   character(param_name_len) :: param_name
   character(2000)           :: param_value
   character(2), parameter   :: blanks = " "//achar(9)
   integer :: parfile_unit, istat, ieq, i1,i2,ln

   parfile_unit= open_infile(parfile)

   ln= 1
   read_loop: do
      read (parfile_unit,'(A)',iostat=istat) line
      if (istat .ne. 0) exit     ! end of file
      i1= verify(line,blanks)
      if(i1.eq.0) cycle          ! skip blank lines
      i2= verify(line,blanks,back=.true.)
      line= line(i1:i2)          ! trim tabs and spaces
      if(line(1:1).eq.'#') cycle ! skip comments, which start with '#'
      ieq= index(line,'=')       ! look for an '=' sign
      if(ieq.eq.0) then
         call raise_exception("Could not read parameter file in line # "//&
                              trim(adjustl(int_to_str(ln)))//" :"//NEW_LINE("A")//&
                              trim(adjustl(line)) ,"read_param",340003)
      endif
      i2= verify(line(1:ieq-1),blanks,back=.true.)
      if(i2.eq.0) then
         call raise_exception("Could not read parameter file in line # "//&
                              trim(adjustl(int_to_str(ln)))//" :"//NEW_LINE("A")//&
                              trim(adjustl(line)) ,"read_param",340003)
      endif
      param_name= line(1:i2)     ! parse param_name and param_value
      i2= verify(line,blanks,back=.true.)
      i1= ieq-1 + verify(line(ieq:i2),blanks//"=")
      param_value= line(i1:i2)
      !print '(A," = :",A,";")', trim(param_name), trim(param_value)
      call set_param(param_name,param_value)
      ln= ln + 1
   enddo read_loop

   ! Check the parameter for consistency
   call check_param

   ! close the file
   close(parfile_unit)

   return

end subroutine read_param


!>
!! Sets a global parameter param_name to the value, given by its string
!! representation param_value
!!
!! \b Edited:
!!    - M.R. 02.11.22, Implemented more meaningfull error msg, made use of *_params
!!.
subroutine set_param(param_name,param_value)

   implicit none
   character(*), intent(in) :: param_name
   character(*), intent(in) :: param_value
   !
   character(9999)         :: all_possible_par
   character(*), parameter :: integer_params =  &
      ":out_every:snapshot_every:nrdiag_every:mainout_every:iwformat:timescales_every" // &
      ":nuflag:fissflag:termination_criterion:flow_every:expansiontype:h_snapshot_every" // &
      ":track_nuclei_every:nr_maxcount:adapt_stepsize_maxcount:extrapolation_width:solver" // &
      ":nse_calc_every:engen_every:top_engen_every:h_mainout_every:h_track_nuclei_every"//&
      ":h_timescales_every:h_flow_every:h_engen_every:gear_nr_maxcount:iwinterp:heating_mode"//&
      ":nr_mincount:gear_nr_mincount:alpha_decay_zmin:alpha_decay_zmax:nse_max_it:screening_mode"//&
      ":nu_loss_every:h_nu_loss_every:interp_mode:nse_solver"
   character(*), parameter :: real_params =  &
      ":temp_reload_exp_weak_rates:engen:initemp_cold:initemp_hot:nsetemp_cold" // &
      ":nsetemp_hot:heating_frac:nse_descend_t9start"    // &
      ":t_analytic:gear_eps:gear_escale:gear_cFactor:gear_nr_eps"// &
      ":timestep_max:timestep_factor:timestep_Ymin"// &
      ":nr_tol:timestep_hydro_factor:final_time:final_temp:final_dens"// &
      ":initial_stepsize:freeze_rate_temp:nse_nr_tol:nse_delt_t9min"//&
      ":heating_density"
   character(*), parameter :: logical_params =  &
      ":read_initial_composition:use_htpf:h_finab" // &
      ":gear_ignore_adapt_stepsize" // &
      ":calc_nsep_energy:timestep_traj_limit"    // &
      ":custom_snapshots:h_custom_snapshots" // &
      ":h_engen_detailed:use_detailed_balance:use_timmes_mue" // &
      ":use_detailed_balance_q_reac" // &
      ":use_tabulated_rates:use_beta_decay_file" //&
      ":use_alpha_decay_file:alpha_decay_ignore_all"//&
      ":use_neutrino_loss_file:use_thermal_nu_loss"
   character(*), parameter :: string_params =  &
      ":trajectory_file:seed_file:net_source:isotopes_file" // &
      ":htpf_file:reaclib_file:fission_rates:weak_rates_file" // &
      ":chem_pot_file:nsep_energies_file:alpha_decay_src_ignore" // &
      ":nunucleo_rates_file:nuchannel_file"      // &
      ":nfission_file:bfission_file"          // &
      ":trajectory_mode:trajectory_format"            // &
      ":track_nuclei_file:nurates_file"          // &
      ":snapshot_file:beta_decay_file:neutrino_mode"  // &
      ":T9_analytic:rho_analytic:Rkm_analytic:Ye_analytic" // &
      ":Le:Lebar:Enue:Enuebar:seed_format" // &
      ":Lx:Lxbar:Enux:Enuxbar:alpha_decay_file" // &
      ":detailed_balance_src_ignore:detailed_balance_src_q_reac"//&
      ":detailed_balance_src_q_winvn" // &
      ":tabulated_rates_file:beta_decay_src_ignore" // &
      ":neutrino_loss_file"

   logical         :: lparam_value
   integer         :: i2
   real(r_kind)    :: score
   character(999)  :: cl_par
   character(2000) :: str_value
   character(500)  :: h_err_msg !< Helper error message


   i2= index(param_value, "#")
   if ((param_value(1:1).eq."'") .or.(param_value(1:1).eq.'"')) &
   then
     str_value= trim(param_value(2:len_trim(param_value)-1))
   elseif(i2.gt.0) then
     str_value= trim(param_value(1:i2-1))
   else
     str_value= trim(param_value)
   endif

   if(len_trim(str_value).ge.5) then
     lparam_value= (str_value(1:5).eq."'yes'") &
               .or.(str_value(1:5).eq.'"yes"')
   elseif(len_trim(str_value).ge.3) then
     lparam_value= (str_value(1:3).eq."yes")
   else
     lparam_value= .false.
   endif

!--- integer parameters
   if(param_name.eq."out_every") then
     out_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."snapshot_every") then
     snapshot_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."h_snapshot_every") then
     h_snapshot_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."h_mainout_every") then
     h_mainout_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."flow_every") then
     flow_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."h_flow_every") then
     h_flow_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."timescales_every") then
     timescales_every = read_integer_param(str_value,param_name)
  elseif(param_name.eq."h_timescales_every") then
     h_timescales_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."engen_every") then
     engen_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."h_engen_every") then
     h_engen_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."nrdiag_every") then
     nrdiag_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."mainout_every") then
     mainout_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."iwformat") then
     iwformat = read_integer_param(str_value,param_name)
   elseif(param_name.eq."iwinterp") then
     iwinterp = read_integer_param(str_value,param_name)
   elseif(param_name.eq."nuflag") then
     nuflag = read_integer_param(str_value,param_name)
   elseif(param_name.eq."fissflag") then
     fissflag = read_integer_param(str_value,param_name)
   elseif(param_name.eq."termination_criterion") then
     termination_criterion = read_integer_param(str_value,param_name)
   elseif(param_name.eq."expansiontype") then
     expansiontype = read_integer_param(str_value,param_name)
   elseif(param_name.eq."nr_maxcount") then
     nr_maxcount = read_integer_param(str_value,param_name)
   elseif(param_name.eq."nr_mincount") then
     nr_mincount = read_integer_param(str_value,param_name)
   elseif(param_name.eq."adapt_stepsize_maxcount") then
     adapt_stepsize_maxcount = read_integer_param(str_value,param_name)
   elseif(param_name.eq."track_nuclei_every") then
     track_nuclei_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."h_track_nuclei_every") then
     h_track_nuclei_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."top_engen_every") then
     top_engen_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."extrapolation_width") then
     extrapolation_width = read_integer_param(str_value,param_name)
   elseif(param_name.eq."solver") then
     solver = read_integer_param(str_value,param_name)
   elseif(param_name.eq."heating_mode") then
     heating_mode = read_integer_param(str_value,param_name)
   elseif(param_name.eq."screening_mode") then
     screening_mode = read_integer_param(str_value,param_name)
   elseif(param_name.eq."interp_mode") then
     interp_mode = read_integer_param(str_value,param_name)
   elseif(param_name.eq."nse_calc_every") then
     nse_calc_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."gear_nr_maxcount") then
     gear_nr_maxcount = read_integer_param(str_value,param_name)
   elseif(param_name.eq."gear_nr_mincount") then
     gear_nr_mincount = read_integer_param(str_value,param_name)
   elseif(param_name.eq."alpha_decay_zmin") then
     alpha_decay_zmin = read_integer_param(str_value,param_name)
   elseif(param_name.eq."alpha_decay_zmax") then
     alpha_decay_zmax = read_integer_param(str_value,param_name)
   elseif(param_name.eq."nse_max_it") then
     nse_max_it = read_integer_param(str_value,param_name)
   elseif(param_name.eq."nse_solver") then
     nse_solver = read_integer_param(str_value,param_name)
   elseif(param_name.eq."nu_loss_every") then
     nu_loss_every = read_integer_param(str_value,param_name)
   elseif(param_name.eq."h_nu_loss_every") then
     h_nu_loss_every = read_integer_param(str_value,param_name)

!--- real parameters
   elseif(param_name.eq."temp_reload_exp_weak_rates") then
     temp_reload_exp_weak_rates= read_float_param(str_value,param_name)
   elseif(param_name.eq."engen") then
     engen= read_float_param(str_value,param_name)
   elseif(param_name.eq."initemp_cold") then
     initemp_cold= read_float_param(str_value,param_name)
  elseif(param_name.eq."initemp_hot") then
     initemp_hot= read_float_param(str_value,param_name)
   elseif(param_name.eq."nsetemp_cold") then
     nsetemp_cold= read_float_param(str_value,param_name)
   elseif(param_name.eq."nse_nr_tol") then
     nse_nr_tol= read_float_param(str_value,param_name)
   elseif(param_name.eq."nse_delt_t9min") then
     nse_delt_t9min= read_float_param(str_value,param_name)
   elseif(param_name.eq."nsetemp_hot") then
     nsetemp_hot= read_float_param(str_value,param_name)
   elseif(param_name.eq."heating_frac") then
     heating_frac= read_float_param(str_value,param_name)
   elseif(param_name.eq."heating_density") then
     heating_density= read_float_param(str_value,param_name)
   elseif(param_name.eq."nse_descend_t9start") then
     nse_descend_t9start= read_float_param(str_value,param_name)
   elseif(param_name.eq."initial_stepsize") then
     initial_stepsize= read_float_param(str_value,param_name)
   elseif(param_name.eq."final_time") then
     final_time= read_float_param(str_value,param_name)
   elseif(param_name.eq."final_temp") then
     final_temp= read_float_param(str_value,param_name)
   elseif(param_name.eq."final_dens") then
     final_dens= read_float_param(str_value,param_name)
   elseif(param_name.eq."t_analytic") then
     t_analytic= read_float_param(str_value,param_name)
   elseif(param_name.eq."gear_eps") then
     gear_eps= read_float_param(str_value,param_name)
   elseif(param_name.eq."gear_escale") then
     gear_escale= read_float_param(str_value,param_name)
   elseif(param_name.eq."gear_cFactor") then
     gear_cFactor= read_float_param(str_value,param_name)
   elseif(param_name.eq."gear_nr_eps") then
     gear_nr_eps= read_float_param(str_value,param_name)
   elseif(param_name.eq."timestep_max") then
     timestep_max= read_float_param(str_value,param_name)
   elseif(param_name.eq."timestep_factor") then
     timestep_factor= read_float_param(str_value,param_name)
   elseif(param_name.eq."timestep_hydro_factor") then
     timestep_hydro_factor= read_float_param(str_value,param_name)
   elseif(param_name.eq."timestep_Ymin") then
     timestep_Ymin= read_float_param(str_value,param_name)
   elseif(param_name.eq."nr_tol") then
     nr_tol= read_float_param(str_value,param_name)
   elseif(param_name.eq."freeze_rate_temp") then
     freeze_rate_temp= read_float_param(str_value,param_name)
!--- logical type parameters
   elseif(param_name.eq."read_initial_composition") then
     read_initial_composition= lparam_value
   elseif(param_name.eq."calc_nsep_energy") then
     calc_nsep_energy= lparam_value
   elseif(param_name.eq."h_engen_detailed") then
     h_engen_detailed= lparam_value
   elseif(param_name.eq."timestep_traj_limit") then
     timestep_traj_limit= lparam_value
   elseif(param_name.eq."custom_snapshots") then
     custom_snapshots= lparam_value
   elseif(param_name.eq."h_custom_snapshots") then
     h_custom_snapshots= lparam_value
   elseif(param_name.eq."use_htpf") then
     use_htpf= lparam_value
   elseif(param_name.eq."h_finab") then
     h_finab= lparam_value
   elseif(param_name.eq."use_timmes_mue") then
     use_timmes_mue= lparam_value
   elseif(param_name.eq."use_tabulated_rates") then
     use_tabulated_rates= lparam_value
   elseif(param_name.eq."use_beta_decay_file") then
     use_beta_decay_file= lparam_value
   elseif(param_name.eq."use_alpha_decay_file") then
     use_alpha_decay_file= lparam_value
   elseif(param_name.eq."use_detailed_balance") then
     use_detailed_balance= lparam_value
   elseif(param_name.eq."use_detailed_balance_q_reac") then
     use_detailed_balance_q_reac= lparam_value
   elseif(param_name.eq."use_thermal_nu_loss") then
    use_thermal_nu_loss= lparam_value
   elseif(param_name.eq."use_neutrino_loss_file") then
     use_neutrino_loss_file= lparam_value
   elseif(param_name.eq."gear_ignore_adapt_stepsize") then
     gear_ignore_adapt_stepsize= lparam_value
   elseif(param_name.eq."alpha_decay_ignore_all") then
     alpha_decay_ignore_all= lparam_value
!--- string parameters
   elseif(param_name.eq."trajectory_mode") then
     trajectory_mode= trim(str_value)
   elseif(param_name.eq."trajectory_file") then
     trajectory_file= trim(str_value)
  elseif(param_name.eq."T9_analytic") then
     T9_analytic= trim(str_value)
  elseif(param_name.eq."rho_analytic") then
     rho_analytic= trim(str_value)
  elseif(param_name.eq."Rkm_analytic") then
     Rkm_analytic= trim(str_value)
  elseif(param_name.eq."Ye_analytic") then
     Ye_analytic= trim(str_value)
   elseif(param_name.eq."alpha_decay_file") then
     alpha_decay_file= trim(str_value)
   elseif(param_name.eq."beta_decay_file") then
     beta_decay_file= trim(str_value)
   elseif(param_name.eq."seed_file") then
     seed_file= trim(str_value)
   elseif(param_name.eq."seed_format") then
     seed_format= trim(str_value)
   elseif(param_name.eq."snapshot_file") then
     snapshot_file= trim(str_value)
   elseif(param_name.eq."net_source") then
     net_source= trim(str_value)
   elseif(param_name.eq."isotopes_file") then
     isotopes_file= trim(str_value)
   elseif(param_name.eq."htpf_file") then
     htpf_file= trim(str_value)
   elseif(param_name.eq."reaclib_file") then
     reaclib_file= trim(str_value)
   elseif(param_name.eq."fission_rates") then
     fission_rates= trim(str_value)
   elseif(param_name.eq."weak_rates_file") then
     weak_rates_file= trim(str_value)
   elseif(param_name.eq."tabulated_rates_file") then
     tabulated_rates_file= trim(str_value)
   elseif(param_name.eq."chem_pot_file") then
     chem_pot_file= trim(str_value)
   elseif(param_name.eq."nsep_energies_file") then
     nsep_energies_file= trim(str_value)
   elseif(param_name.eq."nuchannel_file") then
     nuchannel_file= trim(str_value)
   elseif(param_name.eq."nurates_file") then
     nurates_file= trim(str_value)
   elseif(param_name.eq."nunucleo_rates_file") then
     nunucleo_rates_file= trim(str_value)
   elseif(param_name.eq."nfission_file") then
     nfission_file= trim(str_value)
   elseif(param_name.eq."bfission_file") then
     bfission_file= trim(str_value)
   elseif(param_name.eq."trajectory_format") then
     trajectory_format= trim(str_value)
   elseif(param_name.eq."track_nuclei_file") then
     track_nuclei_file= trim(str_value)
   elseif(param_name.eq."neutrino_mode") then
     neutrino_mode= trim(str_value)
   elseif(param_name.eq."Le") then
     Le = trim(str_value)
   elseif(param_name.eq."Lebar") then
     Lebar = trim(str_value)
   elseif(param_name.eq."Enue") then
     Enue = trim(str_value)
   elseif(param_name.eq."Enuebar") then
     Enuebar = trim(str_value)
   elseif(param_name.eq."Lx") then
     Lx = trim(str_value)
   elseif(param_name.eq."Lxbar") then
     Lxbar = trim(str_value)
   elseif(param_name.eq."Enux") then
     Enux = trim(str_value)
   elseif(param_name.eq."Enuxbar") then
     Enuxbar = trim(str_value)
   elseif(param_name.eq."beta_decay_src_ignore") then
     beta_decay_src_ignore = trim(str_value)
   elseif(param_name.eq."alpha_decay_src_ignore") then
     alpha_decay_src_ignore = trim(str_value)
   elseif(param_name.eq."detailed_balance_src_ignore") then
     detailed_balance_src_ignore = trim(str_value)
   elseif(param_name.eq."detailed_balance_src_q_reac") then
     detailed_balance_src_q_reac = trim(str_value)
   elseif(param_name.eq."detailed_balance_src_q_winvn") then
     detailed_balance_src_q_winvn = trim(str_value)
   elseif(param_name.eq."neutrino_loss_file") then
     neutrino_loss_file = trim(str_value)
   else
     ! Give a detailed error message
     ! For this first merge all possible parameters:
     all_possible_par = integer_params//real_params//logical_params//string_params

     ! Calculate score for similarity
     call search_close_parameter(trim(adjustl(all_possible_par)),trim(adjustl(param_name)),score,cl_par)
     ! Tell the user if there is a close parameter
     ! Probably he meant this?
     if (score .lt. max(LEN_TRIM(adjustl(param_name)),&
                        LEN_TRIM(adjustl(cl_par)))/2) then
       ! I found something close
       h_err_msg = NEW_LINE('A')//"Did you mean '"//trim(adjustl(cl_par))//"'?"
     else
       ! No clue what the user was thinking, dont give more information
       h_err_msg =""
     end if
     ! Finally raise the exception
     call raise_exception('Unknown parameter: '//trim(adjustl(param_name))//"."&
                          //trim(adjustl(h_err_msg)),"set_param",340004)
   endif
end subroutine set_param




!> Converts a string to an integer
!!
!! If the string is not a valid integer, an error message is raised.
!!
!! @author Moritz Reichert
!! @date 22.01.21
function read_integer_param(input_string,param_name)
   implicit none
   character(len=*),intent(in) :: input_string       !< String from param file
   character(len=*),intent(in) :: param_name         !< Parameter name
   integer                     :: read_integer_param !< Converted integer value from input string
   integer                     :: rstat              !< iostat flag

   !< Convert string to integer
   read(input_string,'(I10)',iostat=rstat) read_integer_param

   ! Raise an exception if converting does not work
   if (rstat .ne. 0) then
      call raise_exception('Could not parse parameter "'//trim(adjustl(param_name))//&
                           '". '//NEW_LINE("A")//'The value "'//&
                           trim(adjustl(input_string))//&
                           '" is not valid for this parameter. '//NEW_LINE("A")//&
                           'This parameter assumes an integer.', &
                           "read_integer_param",&
                           340005)
   end if
end function read_integer_param


!> Converts a string to an float
!!
!! If the string is not a valid float, an error message is raised.
!!
!! @author Moritz Reichert
!! @date 22.01.21
function read_float_param(input_string,param_name)
   implicit none
   character(len=*),intent(in) :: input_string       !< String from param file
   character(len=*),intent(in) :: param_name         !< Parameter name
   real(r_kind)                :: read_float_param   !< Converted float value from input string
   integer                     :: rstat              !< iostat flag

   !< Convert string to float
   read(input_string,*,iostat=rstat) read_float_param

   ! Raise an exception if converting does not work
   if (rstat .ne. 0) then
      call raise_exception('Could not parse parameter "'//trim(adjustl(param_name))//&
                           '". '//NEW_LINE("A")//'The value "'//&
                           trim(adjustl(input_string))//&
                           '" is not valid for this parameter. '//NEW_LINE("A")//&
                           'This parameter assumes a float.', &
                           "read_float_param",&
                           340005)
   end if
end function read_float_param




!>
!! Sets default parameters
!! Parameters are sorted in alphabetical order
!!
subroutine set_default_param
   implicit none
   !
   alpha_decay_ignore_all      = .True.
   alpha_decay_file            = "alpha_decays.dat"
   alpha_decay_src_ignore      = "wc07;wc12;wc17"
   alpha_decay_zmax            = 184
   alpha_decay_zmin            = 50
   adapt_stepsize_maxcount     = 20
   beta_decay_file             = "beta_decays.dat"
   beta_decay_src_ignore       = "wc07;wc12;wc17"
   bfission_file               = "BFISSION"
   calc_nsep_energy            = .false.
   chem_pot_file               = "chem_table.dat"
   custom_snapshots            = .false.
   h_custom_snapshots          = .false.
   engen_every                 = 0
   h_engen_every               = 0
   Enue                        = "4.0"
   Enuebar                     = "5.4"
   Enux                        = "4.0"
   Enuxbar                     = "5.4"
   expansiontype               = 1
   extrapolation_width         = 2
   final_dens                  = 1.e0
   final_temp                  = 1.e-2
   final_time                  = 0.e0
   fissflag                    = 0
   flow_every                  = 0
   freeze_rate_temp            = 1.d-2
   gear_eps                    = 1.0d-3
   gear_escale                 = 1.0d-12
   gear_cFactor                = 0.25
   gear_nr_maxcount            = 10
   gear_nr_mincount            = 1
   gear_nr_eps                 = 1.0d-6
   gear_ignore_adapt_stepsize  = .true.
   h_engen_detailed            = .false.
   h_flow_every                = 0
   h_finab                     = .false.
   h_mainout_every             = 0
   heating_frac                = 0.4d0
   heating_density             = 1d11
   heating_mode                = 0
   htpf_file                   = "htpf.dat"
   initemp_cold                = 9.e0
   initemp_hot                 = 9.e0
   initial_stepsize            = 1.d-12
   interp_mode                 = 5
   isotopes_file               = "winvn"
   iwformat                    = 0
   iwinterp                    = 0
   Le                          = "2.0e51"
   Lebar                       = "2.7e51"
   Lx                          = "2.0e51"
   Lxbar                       = "2.7e51"
   mainout_every               = 1
   net_source                  = "sunet"
   neutrino_loss_file          = "nuloss.dat"
   nfission_file               = "NFISSION"
   nrdiag_every                = 0
   nr_tol                      = 1.e-5
   nr_maxcount                 = 3
   nr_mincount                 = 2
   nse_calc_every              = 1
   nse_delt_t9min              = 1d-16
   nse_max_it                  = 25
   nse_nr_tol                  = 1d-6
   nse_descend_t9start         = 100.0
   nse_solver                  = 0
   nsep_energies_file          = "frdm_sn.dat"
   nsetemp_hot                 = 8.e0
   nsetemp_cold                = 7.e0
   nuflag                      = 0
   neutrino_mode               = 'analytic'
   nuchannel_file              = "nu_channels.dat"
   nu_loss_every               = 0
   h_nu_loss_every             = 0
   nunucleo_rates_file         = "neunucleons.dat"
   nurates_file                = "nucross.dat"
   out_every                   = 10
   reaclib_file                = "Reaclib"
   fission_rates               = "fissionrates_frdm"
   read_initial_composition    = .false.
   rho_analytic                = "1.e12"
   Rkm_analytic                = "50.e0"
   screening_mode              = 1
   seed_file                   = "seed"
   seed_format                 = "Name X"
   snapshot_every              = 0
   h_snapshot_every            = 0
   snapshot_file               = "snapshot_freq.dat"
   solver                      = 0
   t_analytic                  = 0.e0
   T9_analytic                 = "10.e0"
   tabulated_rates_file        = "talysNGrates.dat"
   temp_reload_exp_weak_rates  = 1.d-2
   termination_criterion       = 0
   timescales_every            = 0
   h_timescales_every          = 0
   timestep_max                = 2.0e0
   timestep_factor             = 1.0e-1
   timestep_hydro_factor       = 5.0e-2
   timestep_Ymin               = 1.0e-10
   timestep_traj_limit         = .true.
   track_nuclei_every          = 0
   h_track_nuclei_every        = 0
   top_engen_every             = 0
   track_nuclei_file           = "track_nuclei_file"
   trajectory_format           = "time temp dens rad ye"
   trajectory_mode             = "from_file"
   trajectory_file             = "shock.dat"
   detailed_balance_src_ignore = ""
   detailed_balance_src_q_reac = ""
   detailed_balance_src_q_winvn= ""
   use_htpf                    = .false.
   use_tabulated_rates         = .false.
   use_beta_decay_file         = .false.
   use_alpha_decay_file        = .false.
   use_thermal_nu_loss         = .True.
   use_timmes_mue              = .True.
   use_detailed_balance        = .false.
   use_detailed_balance_q_reac = .false.
   use_neutrino_loss_file      = .false.
   weak_rates_file             = "twr.dat"
   Ye_analytic                 = "0.1e0"

end subroutine set_default_param

!>
!! Output parameters to a file
!!
!! Parameters are sorted in alphabetical order
!!
subroutine output_param
   use file_handling_class
   implicit none
   !
   integer :: ofile
   character(3) :: yesno

   if (VERBOSE_LEVEL .ge. 2) then

     ofile= open_outfile("param.out")
         write(ofile,'(A,I5)') 'adapt_stepsize_maxcount     = ' , adapt_stepsize_maxcount
           write(ofile,'(2A)') 'alpha_decay_ignore_all      = ' , yesno(alpha_decay_ignore_all)
           write(ofile,'(3A)') 'alpha_decay_src_ignore      = "', trim(alpha_decay_src_ignore),'"'
           write(ofile,'(3A)') 'alpha_decay_file            = "', trim(alpha_decay_file),'"'
         write(ofile,'(A,I5)') 'alpha_decay_zmax            = ' , alpha_decay_zmax
         write(ofile,'(A,I5)') 'alpha_decay_zmin            = ' , alpha_decay_zmin
           write(ofile,'(3A)') 'beta_decay_file             = "', trim(beta_decay_file),'"'
           write(ofile,'(3A)') 'beta_decay_src_ignore       = "', trim(beta_decay_src_ignore),'"'
           write(ofile,'(3A)') 'bfission_file               = "', trim(bfission_file),'"'
           write(ofile,'(2A)') 'calc_nsep_energy            = ' , yesno(calc_nsep_energy)
           write(ofile,'(3A)') 'chem_pot_file               = "', trim(chem_pot_file),'"'
           write(ofile,'(2A)') 'custom_snapshots            = ' , yesno(custom_snapshots)
           write(ofile,'(3A)') 'detailed_balance_src_ignore = "', trim(detailed_balance_src_ignore),'"'
           write(ofile,'(3A)') 'detailed_balance_src_q_reac = "', trim(detailed_balance_src_q_reac),'"'
           write(ofile,'(3A)') 'detailed_balance_src_q_winvn= "', trim(detailed_balance_src_q_winvn),'"'
         write(ofile,'(A,I5)') 'engen_every                 = ' , engen_every
           write(ofile,'(3A)') 'Enue                        ="' , trim(Enue),'"'
           write(ofile,'(3A)') 'Enuebar                     ="' , trim(Enuebar),'"'
           write(ofile,'(3A)') 'Enux                        ="' , trim(Enux),'"'
           write(ofile,'(3A)') 'Enuxbar                     ="' , trim(Enuxbar),'"'
         write(ofile,'(A,I1)') 'expansiontype               = ' , expansiontype
         write(ofile,'(A,I4)') 'extrapolation_width         = ' , extrapolation_width
     write(ofile,'(A,es14.7)') 'final_dens                  ='  , final_dens
     write(ofile,'(A,es14.7)') 'final_temp                  ='  , final_temp
     write(ofile,'(A,es14.7)') 'final_time                  ='  , final_time
         write(ofile,'(A,I1)') 'fissflag                    = ' , fissflag
           write(ofile,'(3A)') 'fission_rates               = "', trim(fission_rates),'"'
         write(ofile,'(A,I5)') 'flow_every                  = ' , flow_every
     write(ofile,'(A,es14.7)') 'freeze_rate_temp            ='  , freeze_rate_temp
     write(ofile,'(A,es14.7)') 'gear_cFactor                ='  , gear_cFactor
     write(ofile,'(A,es14.7)') 'gear_eps                    ='  , gear_eps
     write(ofile,'(A,es14.7)') 'gear_escale                 ='  , gear_escale
           write(ofile,'(2A)') 'gear_ignore_adapt_stepsize  = ' , yesno(gear_ignore_adapt_stepsize)
     write(ofile,'(A,es14.7)') 'gear_nr_eps                 ='  , gear_nr_eps
         write(ofile,'(A,I5)') 'gear_nr_maxcount            = ' , gear_nr_maxcount
         write(ofile,'(A,I5)') 'gear_nr_mincount            = ' , gear_nr_mincount
           write(ofile,'(2A)') 'h_custom_snapshots          = ' , yesno(h_custom_snapshots)
           write(ofile,'(2A)') 'h_engen_detailed            = ' , yesno(h_engen_detailed)
         write(ofile,'(A,I5)') 'h_engen_every               = ' , h_engen_every
           write(ofile,'(2A)') 'h_finab                     = ' , yesno(h_finab)
         write(ofile,'(A,I5)') 'h_flow_every                = ' , h_flow_every
         write(ofile,'(A,I5)') 'h_mainout_every             = ' , h_mainout_every
         write(ofile,'(A,I5)') 'h_nu_loss_every             = ' , h_nu_loss_every
         write(ofile,'(A,I5)') 'h_snapshot_every            = ' , h_snapshot_every
         write(ofile,'(A,I5)') 'h_track_nuclei_every        = ' , h_track_nuclei_every
         write(ofile,'(A,I5)') 'h_timescales_every          = ' , h_timescales_every
     write(ofile,'(A,es14.7)') 'heating_density             ='  , heating_density
     write(ofile,'(A,es14.7)') 'heating_frac                ='  , heating_frac
         write(ofile,'(A,I1)') 'heating_mode                = ' , heating_mode
           write(ofile,'(3A)') 'htpf_file                   = "', trim(htpf_file),'"'
     write(ofile,'(A,es14.7)') 'initemp_cold                ='  , initemp_cold
     write(ofile,'(A,es14.7)') 'initemp_hot                 ='  , initemp_hot
     write(ofile,'(A,es14.7)') 'initial_stepsize            ='  , initial_stepsize
           write(ofile,'(3A)') 'isotopes_file               = "', trim(isotopes_file),'"'
         write(ofile,'(A,I1)') 'interp_mode                 = ' , interp_mode
         write(ofile,'(A,I1)') 'iwformat                    = ' , iwformat
         write(ofile,'(A,I1)') 'iwinterp                    = ' , iwinterp
           write(ofile,'(3A)') 'Le                          ="' , trim(Le),'"'
           write(ofile,'(3A)') 'Lebar                       ="' , trim(Lebar),'"'
           write(ofile,'(3A)') 'Lx                          ="' , trim(Lx),'"'
           write(ofile,'(3A)') 'Lxbar                       ="' , trim(Lxbar),'"'
         write(ofile,'(A,I5)') 'mainout_every               = ' , mainout_every
           write(ofile,'(3A)') 'net_source                  = "', trim(net_source),'"'
           write(ofile,'(3A)') 'nfission_file               = "', trim(nfission_file),'"'
           write(ofile,'(3A)') 'neutrino_loss_file          = "', trim(neutrino_loss_file),'"'
           write(ofile,'(3A)') 'neutrino_mode               = ' , trim(neutrino_mode)
         write(ofile,'(A,I5)') 'nr_maxcount                 = ' , nr_maxcount
         write(ofile,'(A,I5)') 'nr_mincount                 = ' , nr_mincount
     write(ofile,'(A,es14.7)') 'nr_tol                      ='  , nr_tol
         write(ofile,'(A,I5)') 'nrdiag_every                = ' , nrdiag_every
         write(ofile,'(A,I5)') 'nse_calc_every              = ' , nse_calc_every
     write(ofile,'(A,es14.7)') 'nse_delt_t9min              ='  , nse_delt_t9min
     write(ofile,'(A,es14.7)') 'nse_descend_t9start         ='  , nse_descend_t9start
         write(ofile,'(A,I5)') 'nse_max_it                  = ' , nse_max_it
     write(ofile,'(A,es14.7)') 'nse_nr_tol                  ='  , nse_nr_tol
         write(ofile,'(A,I5)') 'nse_solver                  = ' , nse_solver
           write(ofile,'(3A)') 'nsep_energies_file          = "', trim(nsep_energies_file),'"'
     write(ofile,'(A,es14.7)') 'nsetemp_cold                ='  , nsetemp_cold
     write(ofile,'(A,es14.7)') 'nsetemp_hot                 ='  , nsetemp_hot
         write(ofile,'(A,I1)') 'nuflag                      = ' , nuflag
           write(ofile,'(3A)') 'nuchannel_file              = "', trim(nuchannel_file),'"'
         write(ofile,'(A,I5)') 'nu_loss_every               = ' , nu_loss_every
           write(ofile,'(3A)') 'nunucleo_rates_file         = "', trim(nunucleo_rates_file),'"'
           write(ofile,'(3A)') 'nurates_file                = "', trim(nurates_file),'"'
         write(ofile,'(A,I5)') 'out_every                   = ' , out_every
           write(ofile,'(3A)') 'reaclib_file                = "', trim(reaclib_file),'"'
           write(ofile,'(2A)') 'read_initial_composition    = ' , yesno(read_initial_composition)
           write(ofile,'(3A)') 'rho_analytic                = "', trim(rho_analytic),'"'
           write(ofile,'(3A)') 'Rkm_analytic                = "', trim(Rkm_analytic),'"'
         write(ofile,'(A,I1)') 'screening_mode              = ' , screening_mode
           write(ofile,'(3A)') 'seed_file                   = "', trim(seed_file),'"'
           write(ofile,'(3A)') 'seed_format                 = "', trim(seed_format),'"'
         write(ofile,'(A,I5)') 'snapshot_every              = ' , snapshot_every
           write(ofile,'(3A)') 'snapshot_file               = "', trim(snapshot_file),'"'
         write(ofile,'(A,I1)') 'solver                      = ' , solver
     write(ofile,'(A,es14.7)') 't_analytic                  ='  , t_analytic
           write(ofile,'(3A)') 'T9_analytic                 = "', trim(T9_analytic),'"'
           write(ofile,'(3A)') 'tabulated_rates_file        = "', trim(tabulated_rates_file),'"'
     write(ofile,'(A,es14.7)') 'temp_reload_exp_weak_rates  = ' , temp_reload_exp_weak_rates
         write(ofile,'(A,I1)') 'termination_criterion       = ' , termination_criterion
         write(ofile,'(A,I5)') 'timescales_every            = ' , timescales_every
     write(ofile,'(A,es14.7)') 'timestep_max                ='  , timestep_max
     write(ofile,'(A,es14.7)') 'timestep_factor             ='  , timestep_factor
     write(ofile,'(A,es14.7)') 'timestep_hydro_factor       ='  , timestep_hydro_factor
     write(ofile,'(A,es14.7)') 'timestep_Ymin               ='  , timestep_Ymin
           write(ofile,'(2A)') 'timestep_traj_limit         = ' , yesno(timestep_traj_limit)
         write(ofile,'(A,I5)') 'top_engen_every             = ' , top_engen_every
         write(ofile,'(A,I5)') 'track_nuclei_every          = ' , track_nuclei_every
           write(ofile,'(3A)') 'track_nuclei_file           = "', trim(track_nuclei_file),'"'
           write(ofile,'(3A)') 'trajectory_file             = "', trim(trajectory_file),'"'
           write(ofile,'(3A)') 'trajectory_format           = "', trim(trajectory_format),'"'
           write(ofile,'(3A)') 'trajectory_mode             = "', trajectory_mode,'"'
           write(ofile,'(2A)') 'use_alpha_decay_file        = ' , yesno(use_alpha_decay_file)
           write(ofile,'(2A)') 'use_beta_decay_file         = ' , yesno(use_beta_decay_file)
           write(ofile,'(2A)') 'use_detailed_balance        = ' , yesno(use_detailed_balance)
           write(ofile,'(2A)') 'use_detailed_balance_q_reac = ' , yesno(use_detailed_balance_q_reac)
           write(ofile,'(2A)') 'use_htpf                    = ' , yesno(use_htpf)
           write(ofile,'(2A)') 'use_neutrino_loss_file      = ' , yesno(use_neutrino_loss_file)
           write(ofile,'(2A)') 'use_tabulated_rates         = ' , yesno(use_tabulated_rates)
           write(ofile,'(2A)') 'use_thermal_nu_loss         = ' , yesno(use_thermal_nu_loss)
           write(ofile,'(2A)') 'use_timmes_mue              = ' , yesno(use_timmes_mue)
           write(ofile,'(3A)') 'weak_rates_file             = "', trim(weak_rates_file),'"'
           write(ofile,'(3A)') 'Ye_analytic                 = "', trim(Ye_analytic),'"'


     close(ofile)
  end if

end subroutine output_param


!>
!! Check for consistency of parameter
!!
!! Some parameters might not work together or might cause inconsistencies.
!! This routine checks the user input and complains.
!!
!! @author: M. Reichert
!!
!! \b Edited:
!!      - 21.05.19, MR
!! .
subroutine check_param
   implicit none

   ! NSE temperatures are not set correct
   if (nsetemp_cold .gt. nsetemp_hot) then
      call raise_exception('The parameter "nsetemp_cold" ('//trim(adjustl(num_to_str(nsetemp_cold)))//&
                           ' GK) should be smaller than "nsetemp_hot" ('//&
                           trim(adjustl(num_to_str(nsetemp_cold)))//' GK).' ,&
                           "check_param", 340006)
   end if

   ! Init temperatures are correct
   if (initemp_cold  .gt. initemp_hot) then
      call raise_exception('The parameter "initemp_cold" ('//trim(adjustl(num_to_str(initemp_cold)))//&
                           ' GK) should be smaller than "initemp_hot" ('//&
                           trim(adjustl(num_to_str(initemp_hot)))//' GK).' ,&
                           "check_param",340007)
   endif

   if (nr_mincount .gt. nr_maxcount) then
      call raise_exception('The parameter "nr_mincount" ('//trim(adjustl(int_to_str(nr_mincount)))//&
                           ') should be smaller than "nr_maxcount" ('//&
                           trim(adjustl(int_to_str(nr_maxcount)))//').' ,&
                           "check_param",340008)
   end if

   if (gear_nr_mincount .gt. gear_nr_maxcount) then
      call raise_exception('The parameter "gear_nr_mincount" ('//trim(adjustl(int_to_str(gear_nr_mincount)))//&
                           ') should be smaller than "gear_nr_maxcount" ('//&
                           trim(adjustl(int_to_str(gear_nr_maxcount)))//').' ,&
                           "check_param",340008)
   end if

   ! Check that termination criterium can be fullfilled
   if ((termination_criterion .eq. 0) .and. (trim(adjustl(trajectory_mode)) .eq. "analytic")) then
    call raise_exception('The parameter "termination_criterion" is set to "0" (after end of trajectory)'//&
                         'but the trajectory mode is set to "analytic".',&
                         "check_param",340009)
   end if

end subroutine check_param




!> Search for similar existing parameters.
!!
!! Takes a list of parameters that are separated by ":" and a comparison
!! parameter and returns the most similar parameter from the list. Furthermore,
!! a score is returned that gives a measure of the similarity of the strings.
!! This function is only called if a user specified parameter was not included
!! in the list of possible parameters.
!!
!! @author M. Reichert
!! @date 02.11.22
subroutine search_close_parameter(possible_pars,input_par,score,cl_par)
  implicit none
  character(*),intent(in)  :: possible_pars !< All possible parameter, long string
                                            ! that is separated by ":"
  character(*),intent(in)  :: input_par     !< Parameter to compare
  character(*),intent(out) :: cl_par        !< Most similar parameter
  real(r_kind),intent(out) :: score         !< Minimum distance between parameters
  character(999)           :: helper_char   !< Helping string
  character(999)           :: ms_par        !< Most similar parameter
  real(r_kind)             :: cl_dist       !< Closest distance
  real(r_kind)             :: h_dist        !< helping distance
  integer                  :: length        !< length of closest word
  integer                  :: k             !< Loop variable

  ms_par  = ""  ! Initialize empty
  cl_dist = 999 ! Initialize with large number
  length  = 999 ! Initialize with large number
  helper_char = "" ! Initialize empty
  ! Start from 2 since first character is a ":"
  do k=2,len(possible_pars)

    ! Separate the parameters that are separated by ":"
    if (possible_pars(k:k) .ne. ":") then
      helper_char = trim(adjustl(helper_char))//possible_pars(k:k)
    end if

    if ((possible_pars(k:k) .eq. ":") .or. (k .eq. len(possible_pars))) then
      ! Calculate Levenshtein distance and check similarity, only take the
      ! most similar one
      h_dist = LevenshteinDistance(trim(adjustl(helper_char)),trim(adjustl(input_par)))
      if (h_dist .lt. cl_dist) then
        cl_dist = h_dist
        ms_par  = trim(adjustl(helper_char))
      end if
      ! Reset the parameter
      helper_char = ""
    end if

  end do

  ! Fill the output
  score = cl_dist
  cl_par = trim(adjustl(ms_par))

end subroutine search_close_parameter


!> Calculates the Levenshtein distance between two strings
!!
!! The Levenshtein distance gives a measure of the similarity between
!! two strings. The similarity is expressed in terms of a float afterwards.
!! This routine follows the pseudocode [here](https://en.wikipedia.org/wiki/Levenshtein_distance)
!!
!! @author M. Reichert
!! @date 02.11.22
function LevenshteinDistance(par1,par2) result(dist)
  implicit none
  character(*),intent(in)                 :: par1      !< Input string 1
  character(*),intent(in)                 :: par2      !< Input string 2
  real(r_kind)                            :: dist      !< Measure of the similarity
  character(1)                            :: upper1    !< Helper character
  character(1)                            :: upper2    !< Helper character
  integer                                 :: i, j, h   !< Loop variables
  integer                                 :: m, n      !< Lengths of par1 and par2
  integer                                 :: stat      !< Allocation status
  real(r_kind)                            :: scost     !< Substitution cost
  real(r_kind),dimension(:,:),allocatable :: d         !< Helper matrix

  ! Get the length of the strings
  m = len(par1)
  n = len(par2)

  allocate(d(0:m,0:n),stat=stat)
  if (stat .ne. 0) call raise_exception('Allocation of "d" failed.',&
                                        "LevenshteinDistance",320001)

  ! Initialize with zeros
  d(:,:) = 0

  ! Fill the first column
  do i=1, m
    d(i, 0) = i
  end do
  ! Fill the first row
  do j=1, n
    d(0, j) = j
  end do

  do j=1,n
    do i=1,m
      if (par1(i:i) .eq. par2(j:j)) then
        scost = 0
      else
        ! Check if it is only a case insensitivity
        h = iachar(par1(i:i))
        if (h>= iachar("a") .and. h<=iachar("z") ) then
             upper1 = achar(iachar(par1(i:i))-32)
        else
             upper1 = par1(i:i)
        end if
        h = iachar(par2(j:j))
        if (h>= iachar("a") .and. h<=iachar("z") ) then
             upper2 = achar(iachar(par2(j:j))-32)
        else
             upper2 = par2(j:j)
        end if

        ! It is not as bad when its only a mismatch in the upper/lower case
        if (upper2 .eq. upper1) then
          scost = 0.2
        else
          scost = 1
        end if
      end if
      ! Calculate the total cost of substitutions and exchanges
      d(i,j) = min(d(i-1,j)+1, d(i, j-1) + 1, d(i-1, j-1) + scost)
    end do
  end do

  ! Return value
  dist = d(m,n)

end function LevenshteinDistance




end module parameter_class


!>
!! Converts .TRUE.->"yes", .FALSE. -> "no"
!!
function yesno(claim) result(quote)
   implicit none
   character(3)       :: quote
   logical,intent(in) :: claim
   if(claim) then
      quote= "yes"
   else
      quote= "no"
   endif
end function
