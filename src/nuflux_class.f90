!> @file nuflux_class.f90
!!
!! The error file code for this file is ***W33***.
!!
!! @brief Module \ref nuflux_class and nu-flux related stuff
!!


!> Contains variables and parameters related to neutrino fluxes
!!
!! The neutrino temperature grid that is used for the reaction rates
!! can be changed here.
#include "macros.h"
module nuflux_class
   use parameter_class, only: neutrino_mode
   use parameter_class, only: Le,Lebar,Enue,Enuebar,unit
   use parameter_class, only: Lx,Lxbar,Enux,Enuxbar
   use parser_module,   only: parse_string
   use global_class,    only: reactionrate_type, nurate_type
   use hydro_trajectory
   use file_handling_class
   implicit none


!> Channel type for neutrino reactions according to
!! [Sieverding et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...865..143S/abstract)
   type,public :: nu_channel_type
   integer                   :: id             !< Channel ID
   integer                   :: n_n            !< number of neutrons emitted
   integer                   :: n_p            !< number of protons emitted
   integer                   :: n_a            !< number of alpha particles emitted
   integer                   :: Zdiff          !< change in atomic number
   integer                   :: Adiff          !< change in mass number
   integer                   :: type           !< 1: Charged current nue, 2: Charged current nuebar,
                                               !  3: Neutral current nux, 4: Neutral current nuxbar
   end type nu_channel_type
   type(nu_channel_type), dimension(:), allocatable :: nu_channels !< Array of neutrino channels

   integer   :: nnu   !< Amount of neutrino reactions

   type(reactionrate_type), dimension(:), allocatable :: rrate_nu !< Neutrino rate array


! --- parameters
   real(r_kind),dimension(4), public   :: fluxnu  !< @brief Neutrino fluxes
                                                  !< This variable is calculated and
                                                  !< set in \ref nuflux.
                                                  !< The dimension 4 accounts
                                                  !< for the different neutrino flavors,
                                                  !< however, only (anti-)electron neutrinos
                                                  !< are implemented at the moment.

   real(r_kind),dimension(4), private  :: tempnu  !< @brief Neutrino temperatures.
                                                  !< This variable is set in \ref nutemp
                                                  !< either from the values of a
                                                  !< trajectory or from an analytic
                                                  !< expression. The dimension 4 accounts
                                                  !< for the different neutrino flavors,
                                                  !< however, only (anti-)electron neutrinos
                                                  !< are implemented at the moment.

   integer, parameter   :: nt_nugrid=7 !< Length of the neutrino rate grid

   real(r_kind),private :: sigtempnu(nt_nugrid) = &
   (/2.8d0,3.5d0,4.0d0,5.0d0,6.4d0,8.0d0,10.0d0/) !< temperature grid for neutrino cross sections
                                                  !< @warning This neutrino temperature grid has to
                                                  !< match with the one in the files neunuclei.dat,
                                                  !< neunucleon.dat, and aneunuclei.dat

   real(r_kind), private                :: rs     !< Schwarzschild radius for M=mns [km], calculated in init_nuflux

   ! Helper variables for reading neutrino rates on heavier nuclei
   type(nurate_type),dimension(:),allocatable,private :: anunuc !< anti-neutrino reactions on nuclides
   type(nurate_type),dimension(:),allocatable,private :: nunuc  !< neutrino reactions on nuclides


!--Debugging
   integer,private             :: lumin_debugfile !< Debug file to write neutrino luminosities,
                                                  !< neutrinospheres, and temperatures.


!...total luminosity
!   cross sections are of order 10^-42 [cm^2]
!   luminosity is of order 10^51 [erg]
!   ==> overall factor of 10^9
!      data lumtot /17.9d9,17.7d9,0.0d0,0.0d0/

! --- variables
   real(r_kind),dimension(:),allocatable,public :: tnue,tnuebar        !< (anti-)electron neutrino temperatures from trajectory
   real(r_kind),dimension(:),allocatable,public :: tnux,tnuxbar        !< (anti-)mu and tau neutrino temperatures from trajectory
   real(r_kind),dimension(:),allocatable,public :: nlume,nlumebar      !< (anti-)electron neutrino number luminosities from trajectory
   real(r_kind),dimension(:),allocatable,public :: nlumx,nlumxbar      !< (anti-)mu and tau neutrino number luminosities from trajectory

   logical,private :: include_nc_reactions !< Flag to include neutrino reactions that are not charged current
   integer,private :: n_nuclei,n_anuclei,n_nucleo, n_cc, n_nc


   character(len=*), parameter, private:: nu_binary_name = "nu_binary.windat" !< Name of the neutrino binary file

   !
   ! Public and private fields and methods of the module.
   !
   public:: &
         init_nuflux, nutemp, nuflux, nucs, merge_neutrino_rates, &
         output_binary_neutrino_reaction_data, calculate_nu_rate
   private:: &
         write_reac_verbose_out, set_nutype, &
         read_binary_neutrino_reaction_data, &
         read_reactions_sieverding, read_neutrino_rates

contains


!>
!! Initialize nuflux module
!!
!! Calculates the Schwarzschild radius
!! and opens debug files and writes headers to them.
!! Furthermore it reads the neutrino reactions
!!
!! \b Edited:
!!    - 12.04.17
!!    - 25.01.21, MR - Moved the reading of the nu-rates to this place
!!  .
subroutine init_nuflux()
  use parameter_class, only: nuflag, use_prepared_network, &
                             prepared_network_path
  implicit none

  ! Count individual reaction rates
  n_nuclei=0; n_anuclei=0; n_nucleo=0; n_cc=0; n_nc=0;

  if (nuflag .ge. 1) then

    if ((nuflag .eq. 3) .or. (nuflag .eq. 4)) then
        include_nc_reactions = .true.
    else
        include_nc_reactions = .false.
    end if

     ! Open debug file
     if (VERBOSE_LEVEL .gt. 2.) then
        lumin_debugfile = open_outfile('debug_lumin.dat')
        write(lumin_debugfile,*) 'Debug file of nuflux_class.f90'
        ! TODO Add better description
        write(lumin_debugfile,'(A)') '   Time [s],  NLume [1/s],   NLumebar [1/s], Fnue, Fnuebar, R[cm]'
     endif

     if (use_prepared_network) then
        call read_binary_neutrino_reaction_data(prepared_network_path)
     else
        !-- Read neutrino reactions
        call read_neutrino_rates()
     end if

     !-- Give output statistics
     call write_reac_verbose_out()
  end if

end subroutine init_nuflux


!>
!! Determines neutrino flux for current time (in units of cm^-2 s^-1).
!!
!! The flux is hereby either based on an analytic forumlas of the
!! neutrino luminosities and energies or on the trajectory file.
!!
!! @note In previous versions, also GR corrections were applied here, c.f.,
!! [Otsuki et al. 2000](https://ui.adsabs.harvard.edu/abs/2000ApJ...533..424O/abstract).
!!
!! \b Edited:
!!         - 11.01.14
!!         - MR 18.12.20
!! .
subroutine nuflux(time, rkm)
  use parameter_class, only: unit
  use inter_module
  implicit none

  real(r_kind), intent(in)   :: time   !< current time
  real(r_kind)               :: rkm    !< radius in km
  !
  integer                    :: i
  integer                    :: itime    !index for interpolation
  integer                    :: inuf     !neutrino flavor
  real(r_kind),dimension(4)  :: c0nu     !constant for all 4 nu-flavors
  real(r_kind)               :: frac     !
  real(r_kind)               :: rcm      !radius in cm
  real(r_kind)               :: r2       !4*pi*r*r
  real(r_kind), dimension(4) :: nlum     !current neutrino luminosities
  real(r_kind), dimension(4) :: fluxcom  !current neutrino fluxes for all 4 flavours

  ! Interpolate (or convert from energies to) neutrino temperatures
  call nutemp(time)

  ! account for cross section units of 10^-42 cm^2
  c0nu = 1.0d-42

  !...radius in units of cm; include factor 4*pi
  rcm = rkm*1.d5
  r2=4.d0*unit%pi*rcm*rcm


  select case (neutrino_mode)

   ! for use with analytic neutrino quantities
   case ("analytic")
!...initialize variables

     nlum(1) = parse_string(Le,time)     ![erg/s]
     nlum(2) = parse_string(Lebar,time)  ![erg/s]

     ! Also include heavier neutrinos
     if (include_nc_reactions) then
        nlum(3) = parse_string(Lx,time)     ![erg/s]
        nlum(4) = parse_string(Lxbar,time)  ![erg/s]
     end if

   ! Use data from luminosity file
   case ('from_file')

    ! Interpolate the neutrino luminosities
     itime = 0
     if ( time > ztime(1) .and. time <= ztime(zsteps) ) then

        do i=1,zsteps
           if ( time <= ztime(i) ) exit
        end do
        itime=i
        frac=(time-ztime(itime-1))/(ztime(itime)-ztime(itime-1))
        nlum(1)=nlume(itime-1)+frac*(nlume(itime)-nlume(itime-1))
        nlum(2)=nlumebar(itime-1)+frac*(nlumebar(itime)-nlumebar(itime-1))
        ! Also include other neutrino flavors
        if (include_nc_reactions) then
           nlum(3)=nlumx(itime-1)+frac*(nlumx(itime)-nlumx(itime-1))
           nlum(4)=nlumxbar(itime-1)+frac*(nlumxbar(itime)-nlumxbar(itime-1))
        end if

     else if ( time > ztime(zsteps) ) then
        nlum(1) = nlume(zsteps)
        nlum(2) = nlumebar(zsteps)
        ! Also include other neutrino flavors
        if (include_nc_reactions) then
           nlum(3) = nlumx(zsteps)
           nlum(4) = nlumxbar(zsteps)
        end if
     else if ( time .eq. ztime(1) ) then
        nlum(1) = nlume(1)
        nlum(2) = nlumebar(1)
        ! Also include other neutrino flavors
        if (include_nc_reactions) then
           nlum(3) = nlumx(1)
           nlum(4) = nlumxbar(1)
        end if
     end if
     nlum(1)=dexp(nlum(1))
     nlum(2)=dexp(nlum(2))
     ! Also include other neutrino flavors
     if (include_nc_reactions) then
        nlum(3)=dexp(nlum(3))
        nlum(4)=dexp(nlum(4))
      end if

     if (nlum(1).lt.1e-20) nlum(1)=0.d0
     if (nlum(2).lt.1e-20) nlum(2)=0.d0
     if (include_nc_reactions) then
        if (nlum(3).lt.1e-20) nlum(3)=0.d0
        if (nlum(4).lt.1e-20) nlum(4)=0.d0
     end if


  end select

  ! Calculate number luminosities
  ! Avoid neutrino reactions with very small energies
  do inuf=1,2
    if (tempnu(inuf) .le. 1e-20) then
        nlum(inuf) = 0.d0
    else
        nlum(inuf) = nlum(inuf)*unit%ergtomev/(tempnu(inuf)*3.151374374)
    end if
  end do

  if (include_nc_reactions) then
    do inuf=3,4
      if (tempnu(inuf) .le. 1e-20) then
          nlum(inuf) = 0.d0
      else
          nlum(inuf) = nlum(inuf)*unit%ergtomev/(tempnu(inuf)*3.151374374)
      end if
    end do
  end if


  fluxcom=0.d0
  fluxnu=0.d0
  do inuf=1,2
     fluxcom(inuf) = c0nu(inuf) * nlum(inuf) / r2
     fluxnu(inuf)  = fluxcom(inuf)
  end do

  if (include_nc_reactions) then
     do inuf=3,4
        fluxcom(inuf) = c0nu(inuf) * nlum(inuf) / r2
        fluxnu(inuf)  = fluxcom(inuf)
     end do
  end if

  ! Verbose output
  if (VERBOSE_LEVEL .gt. 2) then
    write(lumin_debugfile,'(6(es23.14,5x))') &
                           time, nlum(1), nlum(2),fluxnu(1)/c0nu(1),fluxnu(2)/c0nu(2),rcm
  endif

  return
end subroutine nuflux



!> Calculates the cross section times the neutrino flux.
!!
!! This subroutine calculates the cross section times
!! the neutrino flux for a given neutrino reaction. Here,
!! the different neutrino types are distinguised.
!! 1: CC reaction (electron neutrino)
!! 2: CC reaction (anti-electron neutrino)
!! 3: NC reaction (all neutrinos)
!! 4: NC reaction (all anti-neutrinos)
!! .
!!
!! \b Edited:
!!        - 12.02.23 MR: Moved this from the jacobian class here
!! .
subroutine calculate_nu_rate(rrate,rat_calc)
  use global_class, only: nurate
  implicit none
  ! Declare the pass
  type(reactionrate_type),intent(inout)  :: rrate    !< rate instance
  real(r_kind),intent(out)               :: rat_calc !< rate value
  ! Internal variables
  integer                             :: wind      !< Weak index
  type(nurate_type)                   :: nurate_tmp!< temporary nurate instance
  real(r_kind)                        :: av_E      !< Average energy

  ! Save the source flag
  wind       = int(rrate%param(1))
  nurate_tmp = nurate(wind)

  rrate%nu_frac = 0

  ! Get the rate
  if (nurate_tmp%kind .eq. 1) then
    rat_calc = fluxnu(1)*nurate(wind)%rcs_e
    rrate%nu_frac = -1d0*sign(1d0,rrate%q_value)*abs(nurate(wind)%ravE/rrate%q_value)
  else if (nurate_tmp%kind .eq. 2) then
    rat_calc = fluxnu(2)*nurate(wind)%rcs_e
    rrate%nu_frac = -1d0*sign(1d0,rrate%q_value)*abs(nurate(wind)%ravE/rrate%q_value)
  else if (nurate_tmp%kind .eq. 3) then
    rat_calc = fluxnu(1)*nurate(wind)%rcs_e + &
               fluxnu(3)*nurate(wind)%rcs_x
    rrate%nu_frac = -1d0*sign(1d0,rrate%q_value)*abs(nurate(wind)%ravE/rrate%q_value)
  else if (nurate_tmp%kind .eq. 4) then
    rat_calc = fluxnu(2)*nurate(wind)%rcs_e + &
               fluxnu(4)*nurate(wind)%rcs_x
    rrate%nu_frac = -1d0*sign(1d0,rrate%q_value)*abs(nurate(wind)%ravE/rrate%q_value)
  else
    call raise_exception("Neutrino kind not implemented yet. "// &
                         "Got: "//trim(int_to_str(nurate_tmp%kind))//". "// &
                         "Only 1, 2 (CC) and 3, 4 (NC) are supported.", &
                         "calculate_nu_rate", 330003)
  end if


  return
end subroutine calculate_nu_rate



!>
!! Interpolate neutrino cross sections
!!
!! The subroutine lin-log interpolates the neutrino cross sections
!! depending on the neutrino temperature (stored in \ref tempnu).
!! For this it also uses the temperature grid \ref sigtempnu.
!! The final interpolated cross sections are stored in
!! \ref global_class::nurate\%rcs.
!!
!! @returns reaction cross sections for neutrino reactions into
!! \ref global_class::nurate\%rcs.
!!
!! \b Edited:
!!       - 13.01.14
!!       - 27.01.21
!!       - 19.07.22, M.R., removed small bug when at the interpolation boundary
!!       - 22.02.23, M.R., removed possible division by zero
!! .
subroutine nucs()
   use global_class, only: nurate
   implicit none

   integer                    :: i,j,n,k, typ
   integer, dimension(4)      :: isave
   real(r_kind)               :: rcs     !cross section for neutrino reaction
   real(r_kind), dimension(4) :: frac
   real(r_kind)               :: cs1,cs2,t1,t2,avE,avE1,avE2
   real(r_kind), dimension(2) :: csg
   integer                    :: start
   integer                    :: count

   n = size(nurate)
!...interpolate neutrino cross sections
   do j=1,4
      do i=1,nt_nugrid
         if ( tempnu(j) .le. sigtempnu(i) ) exit
      end do
!.....isave gives the upper value for interpolation.....................
      isave(j)=i
      if ((isave(j).le.1).or.(isave(j).gt.nt_nugrid)) cycle
      t2=dlog(sigtempnu(isave(j)))
      t1=dlog(sigtempnu(isave(j)-1))
      frac(j)=(dlog(tempnu(j))-t1)/(t2-t1)
   end do


   do i=1,n
    typ = nurate(i)%kind
    start = typ
    ! We also want to calculate the cs for electron (anti-)neutrinos
    if (typ .gt. 2) start = typ-2
    count = 0
    ! Calculate reactions for electron neutrinos, and muon/tauon neutrinos
    inner_loop: do k=start,4,2
        count = count+1
        if ( tempnu(k) .eq. 0.d0) then
            csg(count) = 0.d0
            avE = 0.d0
        else if (maxval(tempnu(:)) .lt. sigtempnu(1)) then
            ! If all temperatures are below the grid, set the cross section to zero
            ! Otherwise calculate the cross section at the lowest grid point
            csg(count) = 0
            avE = 0
        else if ( isave(k) .eq. 1 ) then
            csg(count) = nurate(i)%cs(1)
            avE = nurate(i)%avE(1)
        else if ( isave(k) .gt. nt_nugrid ) then
            csg(count) = nurate(i)%cs(nt_nugrid)
            avE = nurate(i)%avE(nt_nugrid)
        else
            if (nurate(i)%cs(isave(k)) .eq. 0.0d0 ) then
                cs2=dlog(tiny(rcs))                      ! 2009-12-13
            else
                cs2=dlog(nurate(i)%cs(isave(k)))
            end if

            if (nurate(i)%cs(isave(k)-1) .eq. 0.0d0 ) then  ! 2009-12-13
               cs1=dlog(tiny(rcs))                      ! 2009-12-13
            else                                      ! 2009-12-13
               cs1=dlog(nurate(i)%cs(isave(k)-1))
            end if                                    ! 2009-12-13
            rcs=cs1+frac(k)*(cs2-cs1)
            csg(count)=dexp(rcs)

            ! Calculate average energy of absorbed neutrinos
            if ((nurate(i)%avE(isave(k)-1) .eq. 0.0d0 ) .or. &
                (nurate(i)%avE(isave(k)) .eq. 0.0d0 )) then
                avE = 0.0d0
            else
                avE1 = nurate(i)%avE(isave(k)-1)
                avE2 = nurate(i)%avE(isave(k))
                avE = avE1+frac(k)*(avE2-avE1)
            end if
        end if

        if (typ .le. 2)  exit inner_loop
    end do inner_loop
    nurate(i)%rcs_e = csg(1)
    nurate(i)%rcs_x = csg(2)
    nurate(i)%ravE  = avE
   end do

   return

end subroutine nucs


!>
!! @brief Calculates (anti-) neutrino temperatures.
!!
!! The determination depends on the \ref parameter_class::neutrino_mode.
!! In case of an analytic expression for the neutrino temperatures (or energies)
!! it evaluates and parses a string. In case of a trajectory determining the
!! temperatures, this subroutine lin-log interpolates (inside, e.g.,
!! \ref ztime and \ref tnue). The subroutine also
!! ensures that the temperature does not become negative.
!!
!! \b Edited:
!!         - 13.01.14
!!         - MR: 18.12.20
!! .
subroutine nutemp(time)
  implicit none
  integer, parameter :: iflux=0

  real(r_kind), intent(in)   :: time   !< current time
  !
  integer      :: i
  integer      :: itime     !index of current time
  real(r_kind) :: frac

  ! Initialize neutrino temperatures
  tempnu(:) = 0d0

  ! The determination of neutrino temperatures depends on the
  ! neutrino mode
  select case (neutrino_mode)
   case ('analytic')

     tempnu(1) = parse_string(Enue,time)/3.151374374
     tempnu(2) = parse_string(Enuebar,time)/3.151374374

     ! In case neutral current reactions are included
     if (include_nc_reactions) then
       tempnu(3) = parse_string(Enux,time)/3.151374374
       tempnu(4) = parse_string(Enuxbar,time)/3.151374374
     end if

   !...neutrino temperature (interpolate only if t < t_end(AB)
   case ('from_file')

     tempnu(1)=tnue(1)
     tempnu(2)=tnuebar(1)

     ! Check if heavier neutrinos are included
     if (include_nc_reactions) then
       tempnu(3)=tnux(1)
       tempnu(4)=tnuxbar(1)
     end if

    !...log interpolation
     itime = 0
     if ( time > ztime(1) .and. time < ztime(zsteps) ) then
        do i=1,zsteps
           if ( time <= ztime(i) ) exit
        end do
        itime=i
        frac=(time-ztime(itime-1))/(ztime(itime)-ztime(itime-1))
        tempnu(1)=tnue(itime-1)+frac*(tnue(itime)-tnue(itime-1))
        tempnu(2)=tnuebar(itime-1)+frac*(tnuebar(itime)-tnuebar(itime-1))
        if (include_nc_reactions) then
          tempnu(3)=tnux(itime-1)+frac*(tnux(itime)-tnux(itime-1))
          tempnu(4)=tnuxbar(itime-1)+frac*(tnuxbar(itime)-tnuxbar(itime-1))
        end if
     else if ( time >= ztime(zsteps) ) then
        tempnu(1)=tnue(zsteps)
        tempnu(2)=tnuebar(zsteps)
        if (include_nc_reactions) then
          tempnu(3)=tnux(zsteps)
          tempnu(4)=tnuxbar(zsteps)
        end if
     end if
     tempnu(1)=dexp(tempnu(1))
     tempnu(2)=dexp(tempnu(2))
     if (include_nc_reactions) then
        tempnu(3)=dexp(tempnu(3))
        tempnu(4)=dexp(tempnu(4))
     end if

  end select

  ! Ensure positive temperatures
  if (tempnu(1).lt.1e-20) tempnu(1)=0.d0
  if (tempnu(2).lt.1e-20) tempnu(2)=0.d0

  if (include_nc_reactions) then
      if (tempnu(3).lt.1e-20) tempnu(3)=0.d0
      if (tempnu(4).lt.1e-20) tempnu(4)=0.d0
  end if

  return
end subroutine nutemp


!> Write the amount of individual reactions to the out
!!
!! The rates are always counted, for a certain verbose level they
!! are also printed to the OUT file
!!
!! @author M. Reichert
!! @date 27.01.21
subroutine write_reac_verbose_out()
   use error_msg_class, only: int_to_str,write_data_to_std_out
   implicit none
   character(len=7) :: tmp !< temporary character for pretty output

   if (VERBOSE_LEVEL .ge. 1) then
      call write_data_to_std_out("Amount neutrino rates",int_to_str(nnu))
   elseif (VERBOSE_LEVEL .ge. 2) then
      if (nnu .gt. 0) write(*,"(A)") ""
      if (nnu .gt. 0)       write(*,"(A)") "    Neutrino rates:  "
      if (nnu .gt. 0)       write(*,"(A)") "   |------------------------|"
      tmp = int_to_str(nnu)
      if (nnu .gt. 0)      write(*,"(A)")  "   | Total         :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_nucleo)
      if (n_nucleo .gt. 0) write(*,"(A)")  "   | on nucleons   :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_nuclei)
      if (n_nuclei .gt. 0) write(*,"(A)")  "   | nue on nuclei :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_anuclei)
      if (n_anuclei .gt. 0) write(*,"(A)") "   | anue on nuclei:"//adjustr(tmp)//" |"
      tmp = int_to_str(n_cc)
      if (n_cc .gt. 0)      write(*,"(A)") "   | CC  on nuclei :"//adjustr(tmp)//" |"
      tmp = int_to_str(n_nc)
      if (n_nc .gt. 0)      write(*,"(A)") "   | NC  on nuclei :"//adjustr(tmp)//" |"
      if (nnu .gt. 0)       write(*,"(A)") "   |------------------------|"
      if (nnu .gt. 0) write(*,"(A)") ""
   end if

end subroutine write_reac_verbose_out



!> Read neutrino reactions and fill the global_class::nurate array
!!
!! Depending on the nuflag parameter, different files will be opened and
!! read.
!!
!! @see readnucs
!!
!! \b Edited:
!!     -  23.01.21, MR - Moved it to this subroutine from network_init_module
!!  .
subroutine read_neutrino_rates()
   use parameter_class,     only: nuflag,nunucleo_rates_file,&
                                  nuchannel_file,nurates_file
   use mergesort_module,    only: nurate_ms,rrate_sort
   use global_class,        only: nurate, nurate_type, isotope
   use benam_class,         only: getcoefficients
   use file_handling_class
   implicit none
   ! Internal variables
   integer    :: nunucleo !< Id for neutrino reaction file
   integer    :: i,j      !< Loop variable


   ! Count the amount of neutrino reactions
   nnu = 0
   !----- if neutrino rates are used read them from the respective files
   select case(nuflag)
      ! No neutrinos, just do nothing
      case(0)
        nnu = 0
        continue
      ! nuflag = 1...neutrinos only on neutrons and protons
      case(1)
        nunucleo= open_infile(nunucleo_rates_file)
        call readnucs(nunucleo,1)
        close(nunucleo)
        nnu = size(nurate)
        ! Count individual reaction rates
        n_nucleo = nnu
    ! nuflag = 2, 3, 4 neutrino reactions on heavy nuclei from Sieverding et al. (2018)
    ! 3: Only charged current reactions
    ! 4: Only neutral current reactions
    ! 5: Both charged and neutral current reactions
     case(2:4)
        nunucleo= open_infile(nunucleo_rates_file)
        call readnucs(nunucleo,1)
        close(nunucleo)
        ! Count individual reaction rates
        n_nucleo = size(nurate)

        call read_channels(nuchannel_file)
        call read_reactions_sieverding(nurates_file, nuflag-2)
        n_nuclei = size(nunuc)
        call nurate_ms(nurate,size(nurate),nunuc,size(nunuc),0)
        nnu = size(nurate)
    case default
        call raise_exception("Not implemented neutrino mode (nuflag: "// &
                             int_to_str(nuflag)//").", "read_neutrino_rates",&
                             330005)
  end select

  ! set the neutrino type
  call set_nutype(nnu)

  !----- create rrate-type array of neutrino rates, needed to merge them later
    allocate(rrate_nu(nnu))
    do i=1,nnu
       rrate_nu(i)%group         = 1 ! Put all in chapter 1 even when this may be a lie
       rrate_nu(i)%parts         = 0
       rrate_nu(i)%parts(1:6)    = nurate(i)%parts(1:6)
       rrate_nu(i)%source        = nurate(i)%source
       rrate_nu(i)%ch_amount(:)  = nurate(i)%ch_amount(:)
       rrate_nu(i)%is_reverse    = .false.
       rrate_nu(i)%is_resonant   = .false.
       rrate_nu(i)%cached        = -1
       rrate_nu(i)%reac_type     = rrt_nu
       rrate_nu(i)%reac_src      = rrs_nu
       rrate_nu(i)%param(1)      = dble(i)
       rrate_nu(i)%is_weak       = (nurate(i)%kind .eq. 1) .or. &
                                   (nurate(i)%kind .eq. 2)
       rrate_nu(i)%one_over_n_fac= 1
       ! Calculate the Q-value
       rrate_nu(i)%q_value = 0
       do j=1,6
        if (rrate_nu(i)%parts(j) .eq. 0) exit

        rrate_nu(i)%q_value = rrate_nu(i)%q_value - &
                  isotope(rrate_nu(i)%parts(j))%mass_exc*&
                  rrate_nu(i)%ch_amount(j)
       end do

    end do

    ! sort them (just to be sure)
    call rrate_sort(nnu,rrate_nu)


    ! Dont call getcoefficients. The ch_amount is already
    ! calculated in the nuflux_class and one_over_n_fac
    ! is always 1 for neutrino reactions
    ! call getcoefficients(rrate_nu,nnu)

end subroutine read_neutrino_rates


!> Read the reactions from a file in binary format
!!
!! This subroutine reads neutrino reactions and all
!! neutrino data from a binary file. This file has
!! to be previously prepared.
!!
!! @author M. Reichert
!! @date 21.06.2023
subroutine read_binary_neutrino_reaction_data(path)
    use parameter_class,     only: nuflag, max_fname_len
    use global_class,        only: nurate
    use file_handling_class, only: open_unformatted_infile
    implicit none
    ! Internal variables
    character(len=*), intent(in) :: path !< Path to the output directory
    integer                      :: file_id !< File ID
    integer                      :: alloc_stat !< Allocation status

    if (nuflag .gt. 0) then
        ! Open an unformatted file
        file_id = open_unformatted_infile(trim(adjustl(path))//trim(adjustl(nu_binary_name)))

        read(file_id) nnu
        read(file_id) n_nuclei, n_anuclei, n_nucleo, n_cc, n_nc
        ! Allocate the arrays
        allocate(nurate(nnu),stat=alloc_stat)
        if (alloc_stat .ne. 0) then
            call raise_exception("Could not allocate nurate array.", &
                                 "read_binary_neutrino_reaction_data", 330001)
        end if
        read(file_id) nurate

        close(file_id)
    end if

 end subroutine read_binary_neutrino_reaction_data


!> Write the reactions to a file in binary format
!!
!! This subroutine writes the reactions to a file in binary format.
!! This is done for preparing a folder with all reaction rates that
!! can be read in later on from many tracers.
!!
!! @author M. Reichert
!! @date 21.06.2023
subroutine output_binary_neutrino_reaction_data(path)
    use parameter_class,     only: nuflag
    use global_class,        only: nurate
    use file_handling_class, only: open_unformatted_outfile
    implicit none
    ! Internal variables
    character(len=*), intent(in) :: path !< Path to the output directory
    integer                      :: file_id !< File ID

    if (nuflag .gt. 0) then
        ! Open an unformatted file
        file_id = open_unformatted_outfile(trim(adjustl(path))//trim(adjustl(nu_binary_name)))

        write(file_id) nnu
        write(file_id) n_nuclei, n_anuclei, n_nucleo, n_cc, n_nc
        write(file_id) nurate

        close(file_id)
    end if

 end subroutine output_binary_neutrino_reaction_data


!> Read the reactions from a file in the format of Sieverding et al 2018
!!
!! This subroutine reads the reactions from a file in the format
!! of [Sieverding et al 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...865..143S/abstract).
!! The file contains the reaction rates for charged and neutral
!! current neutrino reactions on heavy nuclei.
!! An example of the file looks like:
!! \file{
!!  000  Z = 2    A = 4   channels:  6
!!  002  1.099778E-03  7.960983E-03  2.308440E-02  1.163089E-01  5.736357E-01  2.100826E+00  6.820661E+00
!!  023  1.145175E-03  7.713879E-03  2.115540E-02  9.541441E-02  4.077819E-01  1.295820E+00  3.633845E+00
!!  ...
!! }
!! Here the first line contains the atomic number, the mass number and the amount of reactions.
!! The second line contains the channel specified in the channel file (see \ref read_channels)
!! followed by the reaction rates. The reaction rates are tabulated on a (neutrino) temperature grid
!! of Tnu (MeV) = 2.800 3.500 4.000 5.000 6.400 8.000 10.000.
!!
!! @author M. Reichert
!! @date 12.02.2023
subroutine read_reactions_sieverding(reaction_file_path, reactype)
    use global_class, only: nurate, nurate_type, ipro, ineu, ihe4
    use benam_class, only : findaz
    use file_handling_class
    implicit none
    character(len=*), intent(in) :: reaction_file_path !< Path to the neutrino reaction file
    integer, intent(in)          :: reactype !< 0: charged current only,
                                             !  1: neutral current only,
                                             !  2: both
    ! Internal variables
    integer :: file_id              !< File ID
    integer :: nreactions           !< Number of reactions
    integer :: reac_count           !< Reaction counter
    integer :: i, j, k              !< Loop variable
    integer :: stat                 !< Status variable
    integer :: istat                !< Status variable
    integer :: nchannels            !< Number of channels
    integer :: Z                    !< Atomic number
    integer :: A                    !< Mass number
    integer :: Z_out                !< Proton number of the product
    integer :: A_out                !< Mass number of the product
    integer :: ch_id                !< Channel ID
    real(r_kind),dimension(7) :: cs !< Cross section
    integer :: nucl_id              !< Network index of reacting nucleus
    integer :: nucl_id_out          !< Network index of product nucleus
    integer :: ind                  !< Helper variable

    ! Open the file
    file_id = open_infile(reaction_file_path)
    ! Initialize the amount of reactions
    nreactions = 0
    ! Count the reactions
    do
        read(file_id,'(9x,i3,6x,i3,12x,i3)',iostat=stat) &
            Z, A, nchannels
        if (stat.ne.0) exit
        ! Check if nucleus is included in the network
        nucl_id = findaz(A,Z)

        if (nucl_id .ne. -1) then
            ! Loop through the reactions
            do j = 1, nchannels
                read(file_id,'(I3, 2x, 7(E12.6,2x))') ch_id, cs(:)
                ! Skip neutral current reactions

                if (((reactype .eq. 0 ) .and. (nu_channels(ch_id)%type .gt. 2))) cycle
                ! Skip charged current reactions
                if (((reactype .eq. 1 ) .and. (nu_channels(ch_id)%type .le. 2))) cycle

                ! Check if the reaction is included
                Z_out = Z-nu_channels(ch_id)%Zdiff
                A_out = A-nu_channels(ch_id)%Adiff
                nucl_id_out = findaz(A_out,Z_out)
                if (nucl_id_out .ne. -1) nreactions = nreactions + 1
            end do
        else
            ! Skip the reactions
            do j = 1, nchannels
                read(file_id,'(I3, 2x, 7(E12.6,2x))')
            end do
        end if
    end do

    ! Allocate the nurate array
    allocate(nunuc(nreactions),stat=istat)
    if (istat.ne.0) call raise_exception('Allocation of "nunuc" failed.',&
                                         "read_reactions_sieverding",330001)

    ! Rewind the file
    rewind(file_id)

    ! Read the reactions
    reac_count = 1
    do
        read(file_id,'(9x,i3,6x,i3,12x,i3)',iostat=stat) &
            Z, A, nchannels
        if (stat.ne.0) exit
        ! Check if nucleus is included in the network
        nucl_id = findaz(A,Z)

        if (nucl_id .ne. -1) then
            ! Loop through the reactions
            do j = 1, nchannels
                read(file_id,'(I3, 2x, 7(E12.6,2x))') ch_id, cs(:)

                ! Skip neutral current reactions
                if (((reactype .eq. 0 ) .and. (nu_channels(ch_id)%type .gt. 2))) cycle
                ! Skip charged current reactions
                if (((reactype .eq. 1 ) .and. (nu_channels(ch_id)%type .le. 2))) cycle

                ! Check if the reaction is included
                Z_out = Z-nu_channels(ch_id)%Zdiff
                A_out = A-nu_channels(ch_id)%Adiff
                nucl_id_out = findaz(A_out,Z_out)

                if (nucl_id_out .ne. -1) then
                    ! Store the reaction
                    nunuc(reac_count)%parts(:)     = 0
                    nunuc(reac_count)%ch_amount(:) = 0
                    nunuc(reac_count)%parts(1)     = nucl_id
                    nunuc(reac_count)%parts(2)     = nucl_id_out
                    nunuc(reac_count)%ch_amount(1) = -1
                    nunuc(reac_count)%ch_amount(2) =  1

                    ! Put protons, neutrons, and alphas in the reaction
                    ! parts array. Use the channel amount to specify
                    ! the number as it can be more particles than allowed
                    ! in the reaclib format
                    ind = 3
                    if (nu_channels(ch_id)%n_p .gt. 0) then
                        nunuc(reac_count)%parts(ind)     = ipro
                        nunuc(reac_count)%ch_amount(ind) = nu_channels(ch_id)%n_p
                        ind = ind + 1
                    endif

                    if (nu_channels(ch_id)%n_n .gt. 0) then
                        nunuc(reac_count)%parts(ind)     = ineu
                        nunuc(reac_count)%ch_amount(ind) = nu_channels(ch_id)%n_n
                        ind = ind + 1
                    endif

                    if (nu_channels(ch_id)%n_a .gt. 0) then
                        nunuc(reac_count)%parts(ind)     = ihe4
                        nunuc(reac_count)%ch_amount(ind) = nu_channels(ch_id)%n_a
                        ind = ind + 1
                    endif

                    ! Set the cross section, source, and type of the reaction
                    nunuc(reac_count)%source = "siev"
                    nunuc(reac_count)%cs(:)  = cs(:)
                    nunuc(reac_count)%avE(:) = 0
                    nunuc(reac_count)%kind   = nu_channels(ch_id)%type

                    ! Count reactions for the output
                    if (nunuc(reac_count)%kind .le. 2) then
                        n_cc = n_cc + 1
                    else
                        n_nc = n_nc + 1
                    end if

                    ! Count the reaction to know the index
                    reac_count = reac_count + 1
                end if
            end do
        else
            ! Skip the entries of rates that are irrelevant
            ! because the nucleus is not contained in the network
            do j = 1, nchannels
                read(file_id,'(I3, 2x, 7(E12.6,2x))')
            end do
        end if

    end do

    ! Say someting in the output
    call write_data_to_std_out('Neutrino reactions on (A,Z)',int_to_str(nreactions))

    ! Debug statement, TODO: Make this within verbose level
    ! do i=1, nreactions
    !     write(44,*) nurate_string(nunuc(i))
    ! end do

    ! Close the file
    close(file_id)

end subroutine read_reactions_sieverding



!> This function returns a string with the reaction information
!!
!! The string starts with the kind of reaction
!! (1: cc with nue, 2: cc wit anue, 3: nc with nux, 4: nc with anux)
!! followed by the reactant and the products seperated by "=>".
!!
!! The function is only used for debugging purposes.
!!
!! @author Moritz Reichert
!! @date 12.02.2023
function nurate_string(nurate) result(reaction_string)
    use global_class,  only: net_names, nurate_type
    use benam_class,   only: get_net_name
    type(nurate_type), intent(in) :: nurate


    character(50) :: reaction_string
    integer       :: i
    logical       :: s_prod

    reaction_string = "("//int_to_str(nurate%kind)//") "//&
                      trim(adjustl(net_names(nurate%parts(1))))

    do i =2, 6
       s_prod = .false.
       ! Make a cut between educts and products
        if (i .eq. 2) then
            reaction_string = trim(adjustl(reaction_string))//" =>"
            s_prod = .true.
        end if

       ! Write the names of the isotopes
       if (nurate%parts(i) .ne. 0) then
          if (.not. s_prod) then
             reaction_string = trim(adjustl(reaction_string))//" + "//&
                               int_to_str(nurate%ch_amount(i))//" "//&
                               trim(adjustl(get_net_name(nurate%parts(i))))
          else
             reaction_string = trim(adjustl(reaction_string))//" "//&
                               int_to_str(nurate%ch_amount(i))//" "//&
                               trim(adjustl(get_net_name(nurate%parts(i))))
          end if
       end if

    end do

    return


end function nurate_string





!> Read the channels from a file in the format of Sieverding et al 2018
!!
!! This subroutine reads in the channels from a file in the format of
!! Sieverding et al 2018. The file is read in line by line and the
!! information is stored in the global variable nuchannels.
!!
!! An example of the file looks like:
!! \file{
!! id    type     particle emission
!! \n
!! 001 nue    (cc)  0p 0n 0a
!! 002 nue    (cc)  1p 0n 0a
!! 003 nue    (cc)  0p 1n 0a
!! 004 nue    (cc)  0p 0n 1a
!! 005 nue    (cc)  2p 0n 0a
!! ...}
!! where the type "cc" stands for charged current and "nc" for neutral
!! current. The particle emission is the number of particles emitted
!! in the reaction. The first line is a header line and is ignored.
!!
!! @see [Sieverding et al 2018](https://iopscience.iop.org/article/10.3847/1538-4357/aadd48),
!!      \ref read_reactions_sieverding
!!
!! @author M. Reichert
!! @date 12.02.2023
subroutine read_channels(channel_file_path)
    use file_handling_class
    implicit none
    character(len=*), intent(in) :: channel_file_path !< path to the channel file
    integer                      :: nchannels   !< amount of channels
    integer                      :: id          !< channel id
    character(len=6)             :: nutype      !< neutrino type string (nue, nuebar, nux, nuxbar)
    character(len=2)             :: reac_type   !< reaction type string (cc or nc)
    integer                      :: dec_p       !< amount of protons 0, +1, or -1
    integer                      :: amount_p    !< amount of protons released
    integer                      :: amount_n    !< amount of neutrons released
    integer                      :: amount_a    !< amount of alpha particles released
    integer                      :: type_id     !< Neutrino reaction type
    integer                      :: filename_id !< file id
    integer                      :: i           !< loop variable
    integer                      :: stat        !< status variable

    ! Open the file and get the id
    filename_id= open_infile(channel_file_path)

    ! Count the amount of channels / lines in the file
    nchannels = 0
    ! Skip the header of the file
    read(filename_id,*)
    read(filename_id,*)
    do
        read(filename_id,*,iostat=stat)
        if (stat.ne.0) exit
        nchannels = nchannels + 1
    end do

    ! Allocate the array
    allocate(nu_channels(nchannels), stat=stat)
    if (stat.ne.0) call raise_exception('Allocation of "nu_channels" failed.', &
                                        "read_channels", 330001)

    ! Read the file again and store the information
    rewind(filename_id)
    read(filename_id,*)
    read(filename_id,*)

    do i=1,nchannels
        read(filename_id,"(i3,x,a6,2x,a2,3x,i1,2x,i1,2x,i1)",iostat=stat) &
             id, nutype, reac_type, amount_p, amount_n, amount_a
        if (stat.ne.0) exit

        ! Check that everything is read in correctly
        if (i .ne. id) then
            call raise_exception("Channel id does not match the line number. "// &
                                 "Expected: "//int_to_str(i)//", got: "//int_to_str(id)//".", &
                                 "read_channels", 330006)
        end if

        ! Charged current reaction
        if (reac_type.eq."cc") then
            if (trim(adjustl(nutype)).eq."nue") then
                type_id = 1
                dec_p   = -1
            elseif (trim(adjustl(nutype)).eq."nuebar") then
                type_id = 2
                dec_p   = 1
            else
                call raise_exception("Unknown neutrino type in nu_channels file. "// &
                                     "Only 'nue' and 'nuebar' are allowed for CC reactions. "// &
                                     "Got: '"//trim(adjustl(nutype))//"'.", &
                                     "read_channels", 330003 )
            end if
        ! Neutral current reaction
        else if (reac_type.eq."nc") then
            if (trim(adjustl(nutype)).eq."nux") then
                type_id = 3
                dec_p   = 0
            elseif (trim(adjustl(nutype)).eq."nuxbar") then
                type_id = 4
                dec_p   = 0
            else
                call raise_exception("Unknown neutrino type in nu_channels file. "// &
                                     "Only 'nux' and 'nuxbar' are allowed for NC reactions. "// &
                                     "Got: '"//trim(adjustl(nutype))//"'.", &
                                     "read_channels", 330003 )
            end if
        else
            call raise_exception("Unknown reaction type in nu_channels file. "// &
                                 "Only 'cc' and 'nc' are allowed. "// &
                                 "Got: '"//trim(adjustl(reac_type))//"'.", &
                                 "read_channels", 330004 )
        end if

        ! Store the information
        nu_channels(i)%id   = id
        nu_channels(i)%type = type_id
        nu_channels(i)%n_n  = amount_n ! Released neutrons
        nu_channels(i)%n_p  = amount_p ! Released protons
        nu_channels(i)%n_a  = amount_a ! Released alphas
        ! Difference of protons and mass number in the products
        nu_channels(i)%Zdiff  = amount_p+dec_p
        nu_channels(i)%Adiff  = amount_n+amount_p+amount_a*4
    end do

    ! Close the file
    close(filename_id)

end subroutine read_channels




!>
!! Reads neutrino-nuclei reaction file
!!
!! This subroutine reads in neutrino reactions depending on the
!! input typ. The type of the call is dependent on parameter_class::nuflag.
!! For nuflag = 1, this subroutine is called with typ 1 only, for e.g.,
!! nuflag = 2 it is called 3 times with types 1, 2, and 3.
!!
!! <table>
!! <caption id="multi_row">typ translation</caption>
!! <tr><th> typ <th> Example filename  <th> Meaning
!! <tr><td> 1   <td> _neunucleons.dat_ <td> Reads (anti-)neutrino reactions on neutrons and protons
!! <tr><td> 2   <td> _neunuclei.dat_   <td> Reads neutrino reactions on nuclei
!! <tr><td> 3   <td> _aneunuclei.dat_  <td> Reads Antineutrino reactions on nuclei
!! </table>
!! The file _neunucleons.dat_ looks like:
!! \file{
!! n    p                             nen
!! 12.37     18.61     23.85     36.32     58.26     89.82    139.06
!! p    n                            nebp
!! 6.96     11.14     14.63     22.82     36.66     55.40     82.45  }
!!
!! The file _neunuclei.dat_:
!!
!! \file{
!! Cross sections for (nue,e-) reactions on nuclei for alpha=0.0
!! \=============================================================
!! \n
!!   T [MeV]              2.8          3.5          4.0         5.0          6.4          8.0          10.0
!! \n
!!  N =  2  Z =  2    (4He  )
!!  Integrated           0.00         0.01         0.03         0.15         0.71         2.50         7.91
!! ... }
!!
!! and the file _aneunuclei.dat_:
!!
!! \file{
!! T [MeV]              2.8          3.5          4.0         5.0          6.4          8.0          10.0
!! \n
!! N =  2  Z =  2    (4He  )
!! Integrated           0.00         0.01         0.02         0.09         0.41         1.29         3.64
!! ...}
!!
!!
!! @returns Array of neutrino reactions, depending on input typ either
!!          global_class::nurate (typ 1), global_class::nunuc (typ 2),
!!          or global_class::anunuc (typ 3).
!!
!! @see [Bruenn et al. 2002](https://arxiv.org/pdf/astro-ph/0211404.pdf),
!!      [Froehlich et al. 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...637..415F/abstract)
!!      [Langanke & Kolbe 2001](https://ui.adsabs.harvard.edu/abs/2001ADNDT..79..293L/abstract)
!!
!! \b Edited:
!!      - 12.01.14
!!      - 23.01.21, MR - moved it from network_init_module
!!
!! @remark This routine is used to read the reactions on nucleons only (typ 1).
!! .
subroutine readnucs(nufile,typ)
  use global_class,    only: nurate, nurate_type
  use benam_class ,    only: benam, ident
  use error_msg_class, only: write_data_to_std_out, int_to_str
  implicit none

  integer, intent(in)    :: nufile    !< neutrino cross section file unit
  integer, intent(in)    :: typ       !< type of nu-reaction: 1 -> nu/anu on n/p;
                                      !!                      2 -> nu on nuclide;
                                      !!                      3 -> anu on nuclide.
  !
  integer                            :: i
  integer                            :: stat
  integer                            :: len
  integer                            :: err
  character(5), dimension(2)         :: nam
  integer, dimension(2)              :: parts,ind
  real(r_kind), dimension(nt_nugrid) :: cs
  type(nurate_type), dimension(:), allocatable  :: nuc_tmp


!----- determine approximate number of entries in file for typ=2 or 3
  if (typ.ne.1) then
     len=1
     do
        read(nufile,'(2/)',iostat=stat)
        if (stat.ne.0) exit
        len = len+1
     end do
     allocate(nuc_tmp(len))
     rewind(nufile)
  end if
!----- read nufile for the different cases
  select case(typ)
  case(1)  ! only nu on n,p
     allocate(nurate(2))
     do i=1,2
        nurate(i)%parts(:) = 0
        read(nufile,111)nam(1:2),nurate(i)%source,nurate(i)%cs(1:nt_nugrid),nurate(i)%avE(1:nt_nugrid)
        nurate(i)%parts(1) = benam(nam(1))
        nurate(i)%parts(2) = benam(nam(2))
        nurate(i)%ch_amount(:) = 0
        nurate(i)%ch_amount(1) = -1
        nurate(i)%ch_amount(2) =  1
     end do
  case(2)  ! nu on nuclides
     read(nufile,'(3/)')
     i=0
     do
        read(nufile,112,iostat=stat)parts(1:2),cs(1:nt_nugrid)
        if(stat.ne.0) exit
        call ident(parts(2),parts(1),parts(2)+1,parts(1)-1 &
             ,nam(1),nam(2),ind(1),ind(2),err)
        if (err.eq.0) then
            i = i+1
           nuc_tmp(i)%parts(:) = 0
           nuc_tmp(i)%parts(1:2) = ind
           nuc_tmp(i)%ch_amount(:) = 0
           nuc_tmp(i)%ch_amount(1) = -1
           nuc_tmp(i)%ch_amount(2) =  1
           nuc_tmp(i)%source = 'lznu'
           nuc_tmp(i)%cs = cs
           nuc_tmp(i)%avE = 0
        end if
     end do
     allocate(nunuc(i))
     nunuc = nuc_tmp(1:i)
     deallocate(nuc_tmp)
     call write_data_to_std_out('Neutrino reactions on (A,Z)',int_to_str(i))
  case(3)  ! anu on nuclides
     read(nufile,*)
     i=0
     do
        read(nufile,112,iostat=stat)parts(1:2),cs(1:nt_nugrid)
        if(stat.ne.0) exit
        call ident(parts(2),parts(1),parts(2)-1,parts(1)+1 &
             ,nam(1),nam(2),ind(1),ind(2),err)
        if (err.eq.0) then
            i = i+1
           nuc_tmp(i)%parts(:) = 0
           nuc_tmp(i)%parts(1:2) = ind
           nuc_tmp(i)%ch_amount(:) = 0
           nuc_tmp(i)%ch_amount(1) = -1
           nuc_tmp(i)%ch_amount(2) =  1
           nuc_tmp(i)%source = 'lzan'
           nuc_tmp(i)%cs = cs
           nuc_tmp(i)%avE = 0
        end if
     end do
     allocate(anunuc(i))
     anunuc = nuc_tmp(1:i)
     deallocate(nuc_tmp)
     call write_data_to_std_out('Anti-neutrino reactions on (A,Z)',int_to_str(i))
  end select
  return

111 format(5x,2a5,28x,a4/7(f6.2,4x)/7(f6.2,4x))
112 format(/4x,i3,5x,i3/13x,7(f13.2))

end subroutine readnucs




!> Routine to merge neutrino rates into rrate array
!!
!! This subroutine merges the neutrino rates into a larger rate
!! array via a mergesort (see \ref mergesort_module::rrate_ms)
!!
!! @remark This used to be done after reading chapter 2 in the reaclib
!!         while the reaclib file is still open. Now it is done
!!         afterwards which may be not as efficient, but more clear.
!!
!! @author M. Reichert
!! @date 27.01.21
subroutine merge_neutrino_rates(rrate_array,rrate_length)
   use error_msg_class,  only: raise_exception
   use parameter_class,  only: nuflag
   use mergesort_module, only: rrate_ms,rrate_sort
   use parameter_class,  only: use_prepared_network
   implicit none
   type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate_array  !< Large rate array, containing all reactions
   type(reactionrate_type),dimension(:),allocatable               :: rrate_tmp    !< Temporary rate array
   integer,intent(inout)                                          :: rrate_length !< length of rrate_array
   integer                                                        :: alloc_stat   !< Allocation state
   integer                                                        :: new_length   !< New length of rrate_array

   if (.not. use_prepared_network) then
    if (nuflag .ge. 1) then
        new_length = rrate_length+nnu
        if (nnu .ne. 0) then
            if (.not. allocated(rrate_array)) then
            !-- Allocate the reaclib rate array
            allocate(rrate_array(nnu),stat=alloc_stat)
            if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                        "merge_neutrino_rates",330001)
            rrate_array(1:nnu) = rrate_nu(1:nnu)
            rrate_length = new_length
            else
            call rrate_sort(rrate_length,rrate_array(1:rrate_length))
            call rrate_ms(rrate_array(1:rrate_length),rrate_length,&
                            rrate_nu,size(rrate_nu),0,new_length,rrate_tmp)
            rrate_length = new_length
            ! Reallocate rrate_array with new length
            if (allocated(rrate_array)) deallocate(rrate_array)
            allocate(rrate_array(rrate_length),stat=alloc_stat)
            if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                        "merge_neutrino_rates",330001)
            rrate_array(1:rrate_length) = rrate_tmp(1:rrate_length)
            end if
            !-- Deallocate the reaclib rate array
            deallocate(rrate_nu)
        end if
    end if
   end if

end subroutine merge_neutrino_rates




!> Set the type of a neutrino reaction
!!
!! This subroutine modifies the "kind" entry in global_class::nurate
!! based on the neutrino labels. It is only called for
!! reactions in the format of Froehlich et al. 2006 and
!! for neutrino reactions on nucleons.
!! The different kinds of neutrino reactions are:
!!
!! <table>
!! <caption id="multi_row">Neutrino reaction types</caption>
!! <tr><th> kind <th> Explanation
!! <tr><td> 1    <td> Electron neutrino reaction
!! <tr><td> 2    <td> Electron anti-neutrino reaction
!! <tr><td> 3    <td> Muon/tau neutrino reaction
!! <tr><td> 4    <td> Muon/tau anti-neutrino reaction
!! </table>
!! Kind 3 and 4 are not implemented at the moment.
!!
!! @warning This subroutine uses the neutrino reaction labels
!!          nen, lznu, nebp, and lzan hard coded. It is only used for
!!          neutrino reactions on nucleons.
subroutine set_nutype(nnu_in)

  use global_class, only: nurate

  implicit none
  integer, intent(in)            :: nnu_in !< Length of global_class::nurate
  character(4),dimension(4)      :: nutypes
  integer                        :: i,j

  nutypes = (/' nen' , &
       'lznu' , &
       'nebp' , &
       'lzan'   /)
  do i=1,nnu_in
     do j=1,4
        if(nurate(i)%source == nutypes(j))then
           if (j.le.2) then
              nurate(i)%kind = 1
           else if (j.gt.2) then
              nurate(i)%kind = 2
           end if
           cycle
        end if
     end do
  end do

  return

end subroutine set_nutype

end module nuflux_class
