!> @file thermal_neutrino_module.f90
!!
!! Contains the module \ref thermal_neutrino_module
!!

!>
!! @brief The thermal_neutrino_module serves as interface to
!!        the neutrino emission routines from the sneut5.f90 file.
!!
!! The thermal_neutrino_module is based on the fit of
!! [Itoh et al. 1996](https://ui.adsabs.harvard.edu/abs/1996ApJS..102..411I/abstract)
!! accessed via [Cococubed](https://cococubed.com/code_pages/nuloss.shtml).
!!
!! @author M. Reichert
!! @date   12.04.2023
include 'sneut5.f90'
include 'funct_fermi1.f90'
#include "../macros.h"
module thermal_neutrino_module

  ! Variables for output
  real(r_kind),public, save :: snu   = 0 !< Neutrino emissivity (erg/g)
  real(r_kind),public, save :: spair = 0 !< Pair production emissivity (erg/g)
  real(r_kind),public, save :: splas = 0 !< Plasmon emission emissivity (erg/g)
  real(r_kind),public, save :: sphot = 0 !< Photo neutrino emissivity (erg/g)
  real(r_kind),public, save :: sbrem = 0 !< Bremsstrahlung emissivity (erg/g)
  real(r_kind),public, save :: sreco = 0 !< Recombination neutrino emissivity (erg/g)

  !
  ! Public and private fields and methods of the module.
  !
  public:: &
     thermal_neutrinos, thermal_neutrino_init
!   private:: &


contains


!> @brief Initializes the thermal_neutrino_module.
!!
!! This subroutine is called in the initialization of the code.
!! It only contains debug statements.
!!
!! @author M. Reichert
!! @date   12.04.2023
subroutine thermal_neutrino_init()
  use parameter_class, only: use_thermal_nu_loss, &
                             heating_mode
  use file_handling_class
  implicit none
  real(r_kind)  :: temp
  real(r_kind)  :: rho
  real(r_kind)  :: dummy
  real(r_kind)  :: abar
  real(r_kind)  :: ye
  integer       :: id_debug
  real(r_kind)  :: i

  ! Only do something if the parameters are enabled
  if ((use_thermal_nu_loss) .and. (heating_mode .gt. 0)) then
    if (VERBOSE_LEVEL .gt. 1) then
        ! Make everyone aware that thermal neutrinos are used
        write(*,*) "Using thermal neutrino loss"
    end if

    if (VERBOSE_LEVEL .gt. 2) then
        ! Write debug file to check if the neutrino loss is working correctly
        id_debug = open_outfile('debug_thermal_neutrino_loss.dat')

        ! Write header
        write(id_debug,"(A)") "# Thermal neutrino losses for constant temperature and a 50 percent mixture of C12 and O16"
        write(id_debug,"(A1,(A14,3x),*(A15,3x))") "#", "T(GK)",  "rho(g cm^-3)", "snu(erg/g/s)",  "spair(erg/g/s)", &
                                                  "splas(erg/g/s)",  "sphot(erg/g/s)",  "sbrem(erg/g/s)", "sreco(erg/g/s)"
        ! The following will produce the same data as plotted on
        ! https://cococubed.com/code_pages/nuloss.shtml
        ! Set the constant conditions
        temp = 5e-1
        abar = 14
        ye   = 0.5
        do i=0,12,0.01
            rho = 10**i
            call thermal_neutrinos(abar,ye,temp,rho,1d0,dummy)
            write(id_debug,"(*(E15.5E3,3x))") temp,rho,snu,spair,splas,sphot,sbrem,sreco
        end do

        ! Close the file again
        call close_io_file(id_debug,'debug_thermal_neutrino_loss.dat')
    end if
  end if
end subroutine thermal_neutrino_init


!> @brief Calculates the neutrino emissivity.
!!
!! This subroutine is an interface to the subroutine
!! sneut5_aa from the sneut5.f90 file.
!! It calculates the neutrino emissivity and its derivatives
!! with respect to temperature, density, average mass number and
!! average atomic number.
!! The derivatives are not used in the current version of the code.
!! The file sneut5.f90 is based on the fit of
!! [Itoh et al. 1996](https://ui.adsabs.harvard.edu/abs/1996ApJS..102..411I/abstract)
!! accessed via [Cococubed](https://cococubed.com/code_pages/nuloss.shtml).
!!
!! @author M. Reichert
!! @date   12.04.2023
subroutine thermal_neutrinos(abar,Ye,temp,den,timestep,neutrino_loss)
    implicit none
    real(r_kind),intent(in)  :: abar     !< Average mass number
    real(r_kind),intent(in)  :: Ye       !< Electron fraction
    real(r_kind),intent(in)  :: temp     !< Temperature (GK)
    real(r_kind),intent(in)  :: den      !< Density (g/cm^3)
    real(r_kind),intent(in)  :: timestep !< Timestep (s)
    real(r_kind),intent(out) :: neutrino_loss !< Neutrino emissivity (erg/s/cm^3)
    ! Internal variables
    real(r_kind) :: dsnudt,dsnudd,dsnuda,dsnudz !< Derivatives of neutrino emissivity
    real(r_kind) :: T_K !< Temperature in Kelvin
    real(r_kind) :: zbar !< Average atomic number

    INFO_ENTRY("thermal_neutrinos")

    ! Calculate average atomic number
    zbar = abar*Ye

    ! Convert to Kelvin
    T_K = temp*1.0e9
    call sneut5_aa(T_K,den,abar,zbar, &
                   snu,dsnudt,dsnudd,dsnuda,dsnudz, &
                   spair,splas,sphot,sbrem,sreco)

    spair = spair*timestep !erg/g/s -> erg/g
    splas = splas*timestep !erg/g/s -> erg/g
    sphot = sphot*timestep !erg/g/s -> erg/g
    sbrem = sbrem*timestep !erg/g/s -> erg/g
    sreco = sreco*timestep !erg/g/s -> erg/g

    ! Set output variable
    neutrino_loss = snu*timestep

    INFO_EXIT("thermal_neutrinos")

end subroutine thermal_neutrinos


end module thermal_neutrino_module
