!> Simulation variables for a single zone model
!!
!! The error file code for this file is ***W41***.
!!
!! It contains values of the current step and of the previous one (_p).
!! Additionally it contains \ref evolution_mode that keeps track of
!! the current network state.
!!
!! @author: O. Korobkin
#include "macros.h"
module single_zone_vars
   implicit none
   real(r_kind)                          :: time, time_p           !< Time
   real(r_kind)                          :: stepsize               !< Stepsize
   real(r_kind)                          :: T9, T9_p               !< Temperature [GK]
   real(r_kind)                          :: T9h, T9h_p             !< Temperature [GK] storage for heating mode 2
   real(r_kind)                          :: rhob, rhob_p           !< Density [g/cm^3]
   real(r_kind)                          :: Rkm, Rkm_p             !< Radius of the outflow [km]
   real(r_kind)                          :: Ye, Ye_p               !< Electron fraction [mol/g]
   real(r_kind)                          :: ent, ent_p             !< Entropy [kB/baryon]
   real(r_kind),dimension(:),allocatable :: Y, Y_p                 !< Abundances
   real(r_kind),dimension(:),allocatable :: f, dYdt                !< Time derivative of the abundances
   integer                               :: evolution_mode         !< NSE, network hot/cold, etc.
end module single_zone_vars
