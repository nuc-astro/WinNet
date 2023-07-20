!> @file hydro_trajectory.f90
!!
!! The error file code for this file is ***W25***.
!! This file contains the module \ref hydro_trajectory

!> Contains arrays representing thermodynamic conditions from hydro trajectory file
!!
!! @author O. Korobkin
!!
!! \b Edited:
!!     - MR 26.12.20  - implement init_index
!! .
#include "macros.h"
module hydro_trajectory
   implicit none

   integer                                      :: zsteps     !< number of timesteps in the hydro trajectory
   integer                                      :: init_index !< Initial index in the trajectory
   real(r_kind),dimension(:),allocatable,public :: ztime      !< time information from trajectory
   real(r_kind),dimension(:),allocatable,public :: ztemp      !< temperature information from trajectory
   real(r_kind),dimension(:),allocatable,public :: zdens      !< density information from trajectory
   real(r_kind),dimension(:),allocatable,public :: zye        !< electron fraction information from trajectory
   real(r_kind),dimension(:),allocatable,public :: zrad       !< radii from trajectory
   real(r_kind),dimension(:),allocatable,public :: zvel       !< velocities from trajectory

end module hydro_trajectory
