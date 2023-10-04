!> @file global_class.f90
!!
!! The error file code for this file is ***W23***.
!! @brief Modules: \ref global_class, \ref hydro_trajectory
!!

!> Contains types and objects shared between multiple modules
!!
!! This module is reserved as a container for objects used by
!! multiple modules. In general, it is a good practice to define
!! new variables or types according to their functionality
!! within their respective specialized modules, rather than as
!! globals.
!!
!! @author Oleg Korobkin
!! @date   11.01.2014
!!
#include "macros.h"
module global_class

implicit none

!> data fields for the nuclides contained in the network
type,public                                :: isotope_type
   character(5)                            :: name       !< isotope name
   integer                                 :: mass       !< number of nucleons
   integer                                 :: p_nr       !< number of protons
   integer                                 :: n_nr       !< number of neutrons
   logical                                 :: is_stable  !< true if isotope is stable
   real(r_kind)                            :: spin       !< ground state spin
   real(r_kind)                            :: mass_exc   !< mass excess [MeV]
   real(r_kind),dimension(:),allocatable   :: part_func  !< tabulated partition functions (given by Rauscher XY)
end type isotope_type
type(isotope_type),dimension(:),allocatable,public :: isotope !< all nuclides used in the network

!> Variables related to tracking individual nuclei
integer,dimension(:),allocatable,public :: track_nuclei_indices !< indices of tracked nuclei
integer,public                          :: track_nuclei_nr      !< amount of tracked nuclei

!> Variables related to nuclear heating
logical,public :: heating_switch !< factor to switch on/off nuclear heating

!> reaction rate type
type,public :: reactionrate_type
   integer                   :: group          !< group index of reaction [1:11]
   integer,dimension(6)      :: parts          !< isotope index of reaction participants (0 if empty)
   character(4)              :: source         !< source of reaction rate
   integer,dimension(6)      :: ch_amount      !< defines whether participant is destroyed (-1) or created (+1)
   integer,dimension(4,6)    :: cscf_ind       !< cscf_ind (i,j) gives the position of the entry \f$\frac{d\dot{Y}_{j}}{dY_{i}}\f$ in the cscf data array
   real(r_kind)              :: one_over_n_fac !< 1/n! factor, where n is the number of equal isotopes entering reaction i. Needed to prevent multiple counts of a single reaction.
   logical                   :: is_weak        !< true if reaction is a weak reaction
   logical                   :: is_resonant    !< true if reaction fit corresponds to a resonance [unused]
   logical                   :: is_reverse     !< true if reaction rate is calculated via detailed balance from corresponding forward reaction
   logical                   :: is_const       !< true if reaction rate is constant (param(2:7) .eq. 0.d0)
   integer                   :: reac_type      !< "rrt_sf": spontaneous fission, "rrt_bf": beta-delayed fission, "rrt_nf": neutron-induced fission,
                                               !< "rrt_ng": n-gamma, "rrt_ag":alpha-gamma, "rrt_gn": gamma-n, "rrt_ga": gamma-alpha, "rrt_nu":neutrino ,"rrt_o":other
                                               !< "rrt_betm: beta-minus, "rrt_betp": beta-plus, "rrt_alpd": alpha-decay, "rrt_nemi": neutron emission, "rrt_pemi": proton emission
   integer                   :: reac_src       !< "rrs_reacl": Reaclib, "rrs_tabl": Tabulated, "rrs_nu": neutrino, "rrs_detb": detailed balance , "rrs_twr": theoretical weak rate
   real(r_kind)              :: q_value        !< reaction Q-value [MeV]
   real(r_kind)              :: nu_frac        !< Energy fraction of neutrinos radiated away
   real(r_kind),dimension(9) :: param          !< REACLIB fit parameters
   real(r_kind)              :: cached         !< computed rate
   real(r_kind)              :: cached_p       !< computed rate without density, one_over_nfrac, pf, and screening
end type reactionrate_type
type(reactionrate_type),dimension(:),allocatable,public :: rrate !< array containing all reaction rates used in the network
type(reactionrate_type),dimension(:),allocatable,public :: rrate_weak_exp !< array saving the exp. weak rates from reaclib that are replaced by theo weak rates


!> data fields for neutrino rates given in Langanke&Kolbe 2001
type,public                     :: nurate_type
   integer,dimension(6)         :: parts            !< isotope index of participants
   integer,dimension(6)         :: ch_amount        !< defines whether participant is destroyed
   character(4)                 :: source           !< source of reaction rate
   real(r_kind),dimension(7)    :: cs               !< tabulated reaction cross section
   real(r_kind),dimension(7)    :: avE              !< tabulated average energy of the absorped neutrino
   real(r_kind)                 :: ravE             !< interpolated average energy of the absorped neutrino
   real(r_kind)                 :: rcs_e            !< interpolated cross section for electron neutrinos
   real(r_kind)                 :: rcs_x            !< interpolated cross section for muon and tauon neutrinos
   integer                      :: kind             !< kind of neutrino reaction 1=nue,2=anue,3=numx,4=anumx
end type nurate_type
type(nurate_type),dimension(:),allocatable,public :: nurate !< neutrino rates


!> Flow type to store flows
type,public :: flow_vector
   integer               :: iin !< Index of ingoing nucleus
   integer               :: iout!< Index of outgoing nucleus
   real(r_kind)          :: fwd !< Forward flow
   real(r_kind)          :: bwd !< Backward flow
end type flow_vector

! -- network-related variables
integer,public                               :: net_size          !< total number of isotopes (network size)
integer,public                               :: ihe4,ineu,ipro    !< index of alphas, neutrons and protons
character*5,dimension(:),allocatable,public  :: net_names         !< list of isotopes contained in the network
real(r_kind),dimension(:),allocatable,public :: Qnuloss           !< Qnu for decay of each isotope [MeV]
real(r_kind),dimension(:),allocatable,public :: t9_data           !< temperatures at which partition functions are given [GK]
integer,public                               :: nreac             !< total number of reactions

!> @todo Move the following to tw_rate_module
integer                     :: common_weak_rates !< Counter for rates that are included in Reaclib and theoretical weak rates
integer                     :: only_theo_weak_rates !< Counter for rates that are not included in Reaclib, but in theoretical weak rates

!---- debug
real(r_kind),dimension(5)   :: nag_state  !< t,t9,rho_b for debugging
!---- end debug

end module global_class
