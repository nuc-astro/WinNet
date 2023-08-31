!> @file timestep.f90
!!
!! The error file code for this file is ***W43***.
!!
!! Contains the module \ref timestep_module
!!
!! @brief Timestepping subroutines: \ref timestep, \ref exp_timestep etc.


!>
!! @brief The timestep_module contains all timestep subroutines and newton-raphsons
!! to solve the ODE system.
!!
!! This module contains several different treatments for the timestep and for
!! solving the ODE with different solvers (so far implicit euler and gear).
!! The gear solver develops its own timestep in \ref gear_module. It is only
!! affected by the subroutine \ref restrict_timestep.
!!
#include "macros.h"
Module timestep_module
    use parameter_class, only: timestep_max,timestep_factor,timestep_Ymin,timestep_traj_limit
    use error_msg_class, only: raise_exception, num_to_str
    implicit none

    !
    ! Public and private fields and methods of the module.
    !
    public:: &
         select_optimistic_timestep, advance_implicit_euler, advance_gear
    private:: &
         restrict_timestep, ye_timestep, abchange_timestep, hydro_timestep, &
         update_hydro

contains

!>
!! @brief Network timestepping along the thermodynamic trajectory,
!! specified by {ztime, ztemp, zdens}[zsteps]
!!
!! This subroutine calls more specific subroutines to calculate the
!! timestep. This also depends on the conditions. For example, in NSE it
!! call \ref ye_timestep, while it calls \ref abchange_timestep otherwise.
!!
!! @note Gear solver will only be affected by the last call of this subroutine
!! of \ref restrict_timestep
!!
!! \b Edited:
!!    - OK: 13.01.14
!!    - OK:  7.11.16
!!    - OK: 27.08.17 - merged in the rest of timestepping subroutines
!!    - MR: 27.06.18 - introduced additional subroutine "restrict_timestep"
!! .
subroutine timestep (ctime,temp,dens,entropy,Y,dYdt,Ye,delYe,h,evmode)
   use global_class,     only: net_size
   use parameter_class,  only: solver, freeze_rate_temp
   use gear_module,      only: get_timestep
   implicit none

   real(r_kind), intent(in)                   :: ctime !< current time
   real(r_kind), intent(in)                   :: temp  !< temperature [GK]
   real(r_kind), intent(in)                   :: dens  !< density [g/cm^3]
   real(r_kind), intent(in)                   :: entropy !< [kB/baryon]
   real(r_kind),dimension(net_size),intent(in):: Y     !< abundances
   real(r_kind),dimension(net_size),intent(in):: dYdt  !< change of abundances
   real(r_kind), intent(in)                   :: Ye    !< electron fraction
   real(r_kind), intent(in)                   :: delYe !< change of elec.fraction (Ye-Ye_p)
   real(r_kind), intent(inout)                :: h     !< resulting timestep
   integer, intent(in)                        :: evmode!< zone evolution mode
   !
   real(r_kind) :: h_old
   real(r_kind) :: rad
   real(r_kind) :: tdens =1e24  ! interpolated density: initialize to huge value
   real(r_kind) :: ttemp =1e24  ! interpolated temperature: the same

   INFO_ENTRY("timestep")

   !-- save old timestep
   h_old = h

   ! MR: even without the if statement, the gear solver does not care
   ! about the timestep, so I implemented it here
   if (solver .ne. 1) then
      !-- attempt to increase the timestep
      h = h*timestep_max

      if (evmode.eq.EM_NSE) then
      !-- timescale of Ye changes
         h= min(h, ye_timestep(h, Ye, delYe))
         if (h .ne. h) call raise_exception("Timestep got NaN after ye_timestep.",&
                                            "timestep",430003)

      else
      !-- timescale of abundance change (h_abchange)
         h= min(h, abchange_timestep(h,Y,dYdt))
         if (h .ne. h) call raise_exception("Timestep got NaN after abchange_timestep.",&
                                            "timestep",430004)
      endif

      !-- hydrodynamic timescale (h_hydro)
      !-- Only do this in case temperature is not at lower boundary
      if (temp .gt. freeze_rate_temp*1.000001) then
        h= min(h, hydro_timestep(h,ctime,temp,dens,entropy,ttemp,tdens,Ye))
        if (h .ne. h) call raise_exception("Timestep got NaN after hydro_timestep.",&
                                           "timestep",430005)
      end if
   elseif (solver .eq. 1) then
      ! Also gear should be regulated by the restriction
      h = get_timestep()
      ! Get hydro quantities
      call update_hydro(ctime+h, ttemp, tdens, rad, entropy, Ye)
   end if


   ! Restrict the timestep to match specific values
   call restrict_timestep(ctime,h,temp,dens,ttemp,tdens)
   if (h .ne. h) call raise_exception("Timestep got NaN after restrict_timestep.",&
                                      "timestep",430006)

   INFO_EXIT("timestep")

end subroutine timestep



!>
!! @brief Restricting the previously selected timestep
!!
!! Possible restrictions are given by the termination criterion and
!! the trajectory datapoints as well as custom selected snapshots.
!! Furthermore, it is restricted to hit the NSE transition temperature
!! (\ref parameter_class::nsetemp_cold).
!!
!! @note The timestep of the gear solver, which is evolved independently in
!! \ref gear_module is also effected by this subroutine.
!!
!! \b Edited:
!!    - MR: 27.06.18
!!    - MR: 22.12.20  - Removed "const" parameter
!! .
subroutine restrict_timestep(ctime,h,temp,dens,ttemp,tdens)
   use parameter_class, only: &
       termination_criterion,final_dens,final_temp,final_time,&
       custom_snapshots,trajectory_mode,heating_mode, &
       nsetemp_cold,T9_analytic,rho_analytic,solver,h_custom_snapshots, &
       heating_density
   use gear_module,      only: get_timestep, set_timestep
   use hydro_trajectory, only: ztime,zsteps,ztime,ztemp,zdens
   use analysis,         only: snapshot_time,snapshot_amount
   use single_zone_vars, only: evolution_mode
   use global_class,     only: heating_switch
   use inter_module,     only: inverse_interp1d, interp1d
   use parser_module
   implicit none
   real(r_kind), intent(in)                   :: ctime     !< current time
   real(r_kind), intent(in)                   :: temp      !< temperature [GK]
   real(r_kind), intent(in)                   :: dens      !< density [g/cm^3]
   real(r_kind), intent(in)                   :: ttemp     !< temperature [GK] of the next step
   real(r_kind), intent(in)                   :: tdens     !< density [g/cm^3] of the next step
   real(r_kind), intent(inout)                :: h         !< resulting timestep
   real(r_kind)                               :: h_in      !< ingoing timestep
   real(r_kind)                               :: tnext     !< Next time of the trajectory
   real(r_kind)                               :: temp_tmp  !< Temperature on hypothetical timestep
   real(r_kind)                               :: dens_tmp  !< Density on hypothetical timestep
   integer                                    :: i         !< Loop variable
   logical                                    :: converged !< Convergence indicator for analytic mode

   ! Store the ingoing timestep for later comparison
   h_in = h

   !-- cut by termination criterion
   select case (termination_criterion)
   case (0) ! by the end of trajectory file (regardless of timestep_traj_limit)
      if (ctime+h.gt.ztime(zsteps)) h = ztime(zsteps) - ctime
   case (1) ! by final_time
      if (ctime+h.gt.final_time)    h = final_time - ctime
   case (2) ! by final_temp
      if (ttemp .lt. final_temp*0.999d0) then
         if ((trim(adjustl(trajectory_mode)) .eq. 'from_file')) then
             ! Try to interpolate if still within the trajectory
            if ((ctime.lt.ztime(zsteps)) .and. (.not. (heating_switch))) then
               ! Interpolate the NSE temperature
               tnext = ctime+h
               call inverse_interp1d(zsteps,ztemp,final_temp*0.999d0,ztime,tnext)
               ! Check that it worked, sometimes for very flat trajectories it has
               ! rounding errors
               call interp1d (zsteps,ztime,tnext,ztemp,temp_tmp)
               if (temp_tmp .lt. final_temp) h = max(tnext-ctime,1e-15) ! it worked
               if (temp_tmp .ge. final_temp) h = min(1.1*max(tnext-ctime,1e-15),h) ! make the timestep larger
            else
              ! otherwise, basic interpolation
              h = h * (temp-final_temp)/(temp-ttemp)
            end if
         elseif ((trim(adjustl(trajectory_mode)) .eq. 'analytic')) then
            ! Find the time of T=final_temp
            call find_start_value(T9_analytic,final_temp*0.999d0,ctime+h,converged,tnext)
            if (converged .and. (tnext .gt. ctime)) h = max(tnext-ctime,1e-15)
         end if
      end if
   case (3) ! by final_dens
      if (tdens.le.final_dens*0.999d0) then
         if ((trim(adjustl(trajectory_mode)) .eq. 'from_file')) then
             ! Try to interpolate if still within the trajectory
            if ((ctime.lt.ztime(zsteps))) then
               ! Interpolate the NSE temperature
               tnext = ctime+h
               call inverse_interp1d(zsteps,zdens,final_dens*0.999d0,ztime,tnext)
               ! Check that it worked, sometimes for very flat trajectories it has
               ! rounding errors
               call interp1d (zsteps,ztime,tnext,zdens,temp_tmp)
               if (temp_tmp .lt. final_dens) h = max(tnext-ctime,1e-15) ! it worked
               if (temp_tmp .ge. final_dens) h = min(1.1*max(tnext-ctime,1e-15),h) ! make the timestep larger
            else
              ! otherwise, basic interpolation
              h = h * (dens-final_dens)/(dens-tdens)
            end if
         elseif ((trim(adjustl(trajectory_mode)) .eq. 'analytic')) then
            ! Find the time of rho=final_dens
            call find_start_value(rho_analytic,final_dens*0.999d0,ctime+h,converged,tnext)
            if (converged .and. (tnext .gt. ctime)) h = max(tnext-ctime,1e-15)
         end if
      end if
   endselect


   !-- Prevent NSE overshooting
   if ((evolution_mode .eq. EM_NSE) .and. (ttemp .lt. nsetemp_cold)) then
      if (trajectory_mode .eq. "from_file") then
         ! No heating
         if ((.not. (heating_switch)) .and. (.not. (ctime+h .gt.ztime(zsteps)))) then
            if (ttemp .lt. nsetemp_cold*(1d0-1d-9)) then
               ! Interpolate the NSE temperature
               tnext = ctime+h
               call inverse_interp1d(zsteps,ztemp,nsetemp_cold*(1d0-1d-9),ztime,tnext)
               ! Check that it worked, sometimes for very flat trajectories it has
               ! rounding errors
               call interp1d (zsteps,ztime,tnext,ztemp,temp_tmp)
               if (temp_tmp .lt. nsetemp_cold) h = min(max(tnext-ctime,1e-15),h) ! it worked
               if (temp_tmp .ge. nsetemp_cold) h = min(1.1*max(tnext-ctime,1e-15),h) ! make the timestep larger
            end if
         end if
      else if (trajectory_mode .eq. "analytic") then
         ! No heating
         if (.not. (heating_switch)) then
            if (ttemp .lt. nsetemp_cold*(1d0-1d-9)) then
               ! Find the time of T=T_nsecold
               call find_start_value(T9_analytic,nsetemp_cold*(1d0-1d-9),ctime+h,converged,tnext)
               if (converged .and. (tnext .gt. ctime)) h = max(tnext-ctime,1e-15)
            end if
         end if
      end if
   end if


   !-- Prevent density overshooting for heating switch
   if ((heating_mode .gt. 0) .and. (.not. heating_switch) .and. (tdens .lt. heating_density)) then
    if (trajectory_mode .eq. "from_file") then
        if (tdens .lt. heating_density*(1d0-1d-8)) then
            ! Interpolate the NSE temperature
            tnext = ctime+h
            call inverse_interp1d(zsteps,zdens,heating_density*(1d0-1d-8),ztime,tnext)
             ! Check that it worked, sometimes for very flat trajectories it has
             ! rounding errors
             call interp1d (zsteps,ztime,tnext,zdens,dens_tmp)
             if (dens_tmp .lt. heating_density) h = min(max(tnext-ctime,1e-15),h) ! it worked
             if (dens_tmp .ge. heating_density) h = min(1.1*max(tnext-ctime,1e-15),h) ! make the timestep larger
          end if
    else if (trajectory_mode .eq. "analytic") then
          if (tdens .lt. heating_density*(1d0-1d-8)) then
             ! Find the time of T=T_nsecold
             call find_start_value(rho_analytic,heating_density*(1d0-1d-8),ctime+h,converged,tnext)
             if (converged .and. (tnext .gt. ctime)) h = max(tnext-ctime,1e-15)
          end if
    end if
   end if

   !-- if timestep is limited by trajectory steps, cut to the next grid point
   if (timestep_traj_limit .and. &
       & (trim(adjustl(trajectory_mode)) .eq. 'from_file')) then
         if (ctime.lt.ztime(zsteps)) then
            find_next_time: do i=1,zsteps-1
               if (ztime(i).gt.ctime) exit find_next_time
            enddo find_next_time
            tnext= ztime(i)
            if ((ctime.lt.tnext).and.(ctime+h.gt.tnext)) then
               ! To ensure that the timestep is not to small
               h = max(tnext-ctime,1e-15)
            end if
         endif
   endif

   !-- Jump at least at the last trajectory point
   if ((trim(adjustl(trajectory_mode)) .eq. 'from_file')) then
         if ((ctime.lt.ztime(zsteps)) .and. (ctime+h .gt.ztime(zsteps))) then
               ! To ensure that the timestep is not to small
               h = max(ztime(zsteps)-ctime,1e-15)
         endif
   endif

   !-- Restrict timestep to match custom snapshots
   if ((custom_snapshots) .or. (h_custom_snapshots)) then
      ! Find the next time, should be cheap for short arrays
      next_time: do i=1,snapshot_amount
        if (snapshot_time(i) .gt. ctime) exit next_time
      enddo next_time
      ! Be careful with the array bounds
      if (i .gt. snapshot_amount) i = snapshot_amount
      ! Restrict the timestep if necessary
      tnext=snapshot_time(i)
      if ((ctime.lt.tnext).and.(ctime+h.gt.tnext)) then
         ! To ensure that the timestep is not to small
         h = max(tnext-ctime,1e-15)
      end if
   end if

   ! Also gear should be regulated by the restriction so change the timestep
   if ((solver == 1) .and. (h_in .ne. h)) then
      call set_timestep(h)
   end if
end subroutine restrict_timestep



!>
!! @brief Estimate the timestep from the change of electron fraction
!!
!! In NSE, only (slow) weak reactions are taken into account when calculating the
!! timestep. Therefore, only the change in \f$ y_e \f$ is considered here and
!! the timestep is usually larger than in the full network mode.
!! The timestep is calculated by:
!! \f[ h = h_\mathrm{old} \cdot \eta \frac{y_e}{|\Delta y_e|} \f]
!! with \f$ \eta \f$ being the parameter \ref parameter_class::timestep_factor.
!!
!! \b Created: OK 28.08.2017
!!
function ye_timestep(h,Ye,delYe)
use parameter_class, only: timestep_factor
implicit none
real(r_kind):: ye_timestep
real(r_kind),intent(in):: Ye    !< electron fraction
real(r_kind),intent(in):: delYe !< change of the electron fraction
real(r_kind),intent(in):: h     !< timescale on which this change happens
!
   if (abs(delYe) .gt. 1d-15) then
      ye_timestep  =  h * timestep_factor * Ye/abs(delYe)
   else
      ye_timestep = h
   end if

end function ye_timestep


!>
!! @brief Estimate the timestep from the abundances change
!!
!! In the network evolution mode, the timestep is calculated based on the
!! most changing isotope. The timestep is the maximum of :
!! \f[ h = \eta \cdot \left|\frac{Y}{\Delta Y} \right|  \f]
!! for all nuclei above a certain threshold abundance
!! \ref parameter_class::timestep_Ymin. The factor \f$ \eta \f$ is determined
!! by the parameter \ref parameter_class::timestep_factor.
!!
!! \b Created: OK 28.08.2017
!!
function abchange_timestep(h,Y,dYdt)
use parameter_class, only: timestep_factor,timestep_Ymin
use global_class, only: net_size
implicit none
real(r_kind):: abchange_timestep
real(r_kind),dimension(net_size),intent(in):: Y     !< abundances
real(r_kind),dimension(net_size),intent(in):: dYdt  !< change of abundances
real(r_kind),intent(in):: h     !< timescale on which this change happens
!
integer:: i
real(r_kind):: h1

   h1= h
   do i = 1, net_size
      if (Y(i) .le. timestep_Ymin) cycle
      if (dYdt(i) .eq. 0.d0) cycle
      if (dYdt(i)*h1 .gt. abs(Y(i))*timestep_factor) then
         h1= timestep_factor * abs(Y(i)/dYdt(i))
      endif
   end do
   abchange_timestep= h1

end function abchange_timestep


!>
!! @brief Estimate the hydro timestep
!!
!! The timestep is restricted by the change of the hydrodynamic conditions
!! (i.e., temperature and density). The maximum change of them is determined
!! by the parameter \ref parameter_class::timestep_hydro_factor.
!! If the temperature or density change to rapidly the timestep will be decreased
!! and ultimately, if the timestep is lower than \f$ 10^{-24}s \f$ an exception
!! will be raised and the code terminates. The subroutine is also responsible
!! to update the next temperature (ttemp) and density (tdens) that are used
!! in \ref restrict_timestep.
!!
!! @note This subroutine restricts only the timestep for the implicit euler solver,
!! not for the gear solver.
!!
!! \b Created: OK 28.08.2017
!!
function hydro_timestep(h1,ctime,temp,dens,entropy,ttemp,tdens,Ye)
use parameter_class, only: heating_mode
use parameter_class, only: timestep_hydro_factor
use gear_module,     only: get_timestep
use global_class,    only: heating_switch
use expansion_module,only: rad
implicit none
real(r_kind):: hydro_timestep
real(r_kind),intent(in):: h1    !< timescale on which this change happens
real(r_kind),intent(in):: ctime !< current time
real(r_kind),intent(in):: temp  !< temperature at ctime
real(r_kind),intent(in):: dens  !< density at ctime
real(r_kind),intent(in):: entropy !< [kB/baryon]
real(r_kind),intent(inout):: ttemp !< computed temperature at ctime + h
real(r_kind),intent(inout):: tdens !< computed density at ctime + h
real(r_kind),intent(in):: Ye    !< electron fraction
!
real(r_kind)   :: ttime,h
character(200) :: e_message2  !< Error message


h= h1
hydro_loop: do
   if (h .ne. h) call raise_exception("Timestep got NaN.","hydro_timestep",430007)

   ! Suggested next time
   ttime = ctime + h
   ! Initialize old values
   ttemp = temp
   tdens = dens
   ! Get temperature and density in the next step
   call update_hydro(ttime, ttemp, tdens, rad, entropy, Ye)

   ! check density (and temperature if no heating) increment
   ! is small enough
   if (heating_switch) then
      if (dabs(tdens-dens)/dens.lt.timestep_hydro_factor)   exit hydro_loop
      h = h * timestep_hydro_factor * dens / abs(tdens-dens) &
            * .999 ! black magic
   else
      if ((dabs(tdens-dens)/dens.lt.timestep_hydro_factor) .and. &
          (dabs(ttemp-temp)/temp.lt.timestep_hydro_factor)) exit hydro_loop
      h = h * timestep_hydro_factor  &
            / max (abs(tdens-dens)/dens,abs(ttemp-temp)/temp) &
            * .999 ! black magic
   endif

   ! if timestep too small, complain and exit
   if (h.lt.1e-24) then
      ! Prepare a meaningful message
      ! Say something about the time where it occured
      e_message2 = 'time: t   ='//trim(adjustl(num_to_str(ctime)))//', h     ='//trim(adjustl(num_to_str(h)))//NEW_LINE('A')
      e_message2 = trim(adjustl(e_message2))//'dens: d(t)='//trim(adjustl(num_to_str(dens)))//', d(t+h)='//trim(adjustl(num_to_str(tdens)))//NEW_LINE('A')
      e_message2 = trim(adjustl(e_message2))//'temp: d(t)='//trim(adjustl(num_to_str(temp)))//', d(t+h)='//trim(adjustl(num_to_str(ttemp)))//NEW_LINE('A')

      ! Raise the exception
      call raise_exception(&
         "Timestep too small (<1e-24). This happened because the hydrodynamic"//NEW_LINE('A')//&
         'conditions were changing to fast. Check the hydrodynamic conditions'//NEW_LINE('A')//&
         'or choose the "timestep_hydro_factor" parameter larger ('//trim(adjustl(num_to_str(timestep_hydro_factor)))//').'//&
         NEW_LINE('A')//&
         trim(adjustl(e_message2))&
         ,'hydro_timestep',430008)
   endif
enddo hydro_loop


hydro_timestep= h

end function hydro_timestep


!>
!! Selects the timestep
!!
!! This subroutine is the main interface to the \ref driver.f90 .
!! It calls \ref timestep and initializes the timestep in case
!! of the first call.
!!
!! \b Edited:
!!      - 17.12.2016
!! .
subroutine select_optimistic_timestep(init)
use parameter_class, only: initial_stepsize, solver
use single_zone_vars,only: time,Y,Ye,Ye_p,T9,rhob,Rkm,dYdt,ent
use single_zone_vars,only: stepsize,evolution_mode
use nucstuff_class,  only: el_ab
use jacobian_class,  only: abchange
use gear_module, only: init_dYdt
implicit none
logical, intent(inout):: init

   if (init) then

      ! calculate dYdt from reaction rates, needed for timestep
      call el_ab (Y, Ye)
      call abchange (time, T9, rhob, Ye, Rkm, Y, dYdt, evolution_mode)

      if(solver==1) call init_dYdt(dYdt)

      ! adjust the stepsize
      stepsize= initial_stepsize
   endif

   call el_ab (Y, Ye)
   call timestep (time,T9,rhob,ent,Y,dYdt,Ye,Ye-Ye_p,stepsize,evolution_mode)
   init = .false.

end subroutine select_optimistic_timestep



!>
!! @brief Returns temperature, density, and radius for a given time "time_i".
!!
!! The time, entropy, and electron fraction are only inputs and are not changed.
!! The update depends on the trajectory mode of the run and can be either
!! interpolated from a trajectory, got from the expansion, or from an analytic
!! expression.
!!
!! \b Edited:
!!       - MR: 15.01.2020 - Moved it from advance_implicit_euler/gear to this new subroutine
!!       - MR: 10.02.2023 - Gave an optional output for the temperature
!! .
subroutine update_hydro(time_i, temperature, density, radius, entropy, efraction, T_tr)
   use single_zone_vars, only: Y, T9_p, time_p
   use parameter_class,  only: trajectory_mode, heating_mode,&
                               rho_analytic, Rkm_analytic, T9_analytic,&
                               freeze_rate_temp
   use inter_module,     only: interp1d
   use hydro_trajectory, only: zsteps,ztime,ztemp,zdens,zrad
   use expansion_module, only: expansion,vel
   use global_class,     only: net_size, isotope, heating_switch
   use parser_module,    only: parse_string
   implicit none
   real(r_kind),intent(in)    :: time_i      !< Time [s]
   real(r_kind),intent(in)    :: entropy     !< Entropy [kB/nuc]
   real(r_kind),intent(in)    :: efraction   !< Electron fraction
   real(r_kind),intent(inout) :: temperature !< temperature at time_i
   real(r_kind),intent(inout) :: density     !< density at time_i
   real(r_kind),intent(inout) :: radius      !< radius at time_i
   real(r_kind),intent(out),optional :: T_tr !< temperature at time_i, ignoring heating
   real(r_kind)               :: ysum, yasum
   real(r_kind)               :: htmp        !< Temporary stepsize


   ! Calculate necessary sums for expansion
   ysum    = sum(Y(1:net_size))
   yasum   = sum(Y(1:net_size)*isotope(1:net_size)%mass)
   ! yzsum   = sum(Y(1:net_size)*isotope(1:net_size)%p_nr)

   ! interpolate temperature, density and radius
   if (trim(adjustl(trajectory_mode)) .eq. "analytic") then
      if (.not. (heating_switch)) then
         ! Update temperature from analytic expression
         temperature = parse_string(T9_analytic,time_i)
      else
         temperature = T9_p
         if (present(T_tr)) then
            T_tr = parse_string(T9_analytic,time_i)
         end if
      end if
      ! Update density and radius from analytic expression
      density = parse_string(rho_analytic,time_i)
      radius  = parse_string(Rkm_analytic,time_i)
      ! Avoid underflows/overflows
      if (temperature .le. freeze_rate_temp)  temperature = freeze_rate_temp
      if (density     .le. 1d-99)             density     = 1d-99
      if (radius      .ge. 1d99)              radius      = 1d99

   else if (trim(adjustl(trajectory_mode)) .eq. "from_file") then
      if (time_i.gt.ztime(zsteps)) then
         htmp = time_i - time_p
         call expansion(time_p,htmp,efraction,density,radius,temperature,entropy,vel,yasum/ysum)
         if (present(T_tr)) then
            T_tr = temperature
        end if
        if (heating_switch) then
            temperature = T9_p
        end if
      else
         if (.not. (heating_switch)) then
            call interp1d (zsteps,ztime,time_i,ztemp,temperature)
         else
            if (present(T_tr)) then
                call interp1d (zsteps,ztime,time_i,ztemp,temperature)
                T_tr = temperature
            end if
            temperature = T9_p
         end if

         call interp1d (zsteps,ztime,time_i,zdens,density)
         call interp1d (zsteps,ztime,time_i,zrad,radius,itype=itype_LINEAR)
      end if
   end if

   ! Generally only allow temperatures above 1e-6 GK
   temperature = max(temperature,freeze_rate_temp)

end subroutine update_hydro


!>
!! Advance system to the next step for the implicit Euler method
!!
!! This subroutine solves one timestep for the implicit Euler method
!! via a newton-raphson (NR) scheme. The maximum amount of iterations for
!! the NR is controlled by the parameter \ref parameter_class::nr_maxcount.
!! The NR is successful if the mass conservation deviates less than
!! parameter_class::nr_tol. If the maximum iterations have been reached without
!! convergence, it repeats the calculation with half of the timestep.
!! This is done parameter_class::adapt_stepsize_maxcount many times. If after
!! this amount of iterations still no convergence could be achieved it will
!! raise an error.
!!
!! @see [Winteler 2013](https://edoc.unibas.ch/29895/)
!!
!! \b Edited:
!!     - OK: 15.06.2017 - Made into a separate subroutine
!!     - MR: 06.07.2022 - Let it crash if the result does not converge properly
!! .
subroutine advance_implicit_euler(cnt)
use parameter_class, only: adapt_stepsize_maxcount,&
    heating_mode,nr_maxcount,nr_tol,nse_calc_every, &
    trajectory_mode,initial_stepsize,freeze_rate_temp,&
    timestep_hydro_factor,nr_mincount
use global_class,    only: net_size, isotope,ipro,ineu,&
    heating_switch
use hydro_trajectory,only: zsteps,ztime
use analysis,        only: output_nr_diagnostic
use winnse_module,   only: pf,winnse_guess
use nucstuff_class,  only: el_ab,masscalc
use jacobian_class,  only: jacobi_init
use pardiso_class,   only: netsolve
use nuclear_heating, only: nuclear_heating_update
use ls_timmes_eos_module, only: timmes_eos_state,timmes_eos,ink
use error_msg_class, only: num_to_str
use single_zone_vars
implicit none
integer, intent(in)    :: cnt          !< global iteration counter
!
integer                :: k            !< Current iteration of adapt stepsize
integer                :: nr_count     !< Current iteration of the newton-raphson
real(r_kind)           :: m_tot        !< Mass conservation
type(timmes_eos_state) :: state        !< State variable for Helmholtz EOS
integer                :: eos_status   !< Status of the eos, 0=working, 0!=failing
logical                :: exit_condition!< Helper variable that keeps track of the exit condition of the NR
real(r_kind)           :: T9_nrlast    !< Last temperature of the NR
   ! Check the mass conservation and kill the calculation if something goes odd.
   ! Mass conservation should not be within 0.001% to nr_tol
   ! simulataneously with a small timestep.
   ! Note: Other networks such as SkyNet rescale the abundances at this point.
   call masscalc(Y_p, m_tot)
   if ((1.00001*dabs(m_tot-1d0) .ge. nr_tol) .and. (stepsize .le. 1e-15)) then
     call raise_exception("The result did not converge properly. "//NEW_LINE("A")//&
                          "Mass conservation was "//num_to_str(dabs(m_tot-1d0))//"."//NEW_LINE("A")//&
                          "Stepsize was "//num_to_str(stepsize)//"s."//NEW_LINE("A")//&
                          "This error could be a result of inconsistent reaction rates, "//NEW_LINE("A")//&
                          "if this is not the case,"//&
                          "try to use different values of nr_tol, nr_maxcount, or timestep_factor.",&
                          'advance_implicit_euler',430011)
   end if

   adapt_stepsize: do k=0,adapt_stepsize_maxcount

      ! compute some handy reductions and auxiliary variables
      Y = Y_p
      Rkm     = Rkm_p
      time    = time_p + stepsize
      ent     = ent_p
      P       = P_p
      T9      = T9_p

      call el_ab (Y, Ye)
      ! Get temperature, density, and radius
      call update_hydro(time, T9, rhob, Rkm, ent, Ye, T9h)

      ! Set an initial value to the last temperature of the NR
      T9_nrlast = T9

      !-- Newton-Rhapson iterative solve ---------
      newton_raphson: do nr_count=1,nr_maxcount
         exit_condition = .True.

         ! calculate jacobian and rhs of system of DE (f())
         call jacobi_init(time,T9,rhob,Rkm,Y,Y_p,dYdt,f,stepsize,evolution_mode)

         ! solve the linearized system
         call netsolve(f, Y)

         ! thermodynamics update
         call el_ab(Y, Ye)

         ! Check if thinks are reasonable or
         ! went terribly wrong
         if ((Ye .lt. 0) .or. ( Ye .gt. 1)) then
            if (VERBOSE_LEVEL .gt. 1) then
               write(*,*) "Ye went wrong, obtained Ye = "//num_to_str(Ye)
            end if
            exit_condition = .False.
            exit newton_raphson
         end if

         if (heating_switch) then
            ! Update NSE if heating is turned on
            if (evolution_mode .eq. EM_NSE) then
                call winnse_guess( T9, rhob, Ye, Y(ineu), Y(ipro), Y)
            end if
            ! update entropy from nuclear heating, then update temperature
            call nuclear_heating_update(nr_count,rhob,Ye,pf,Y_p,Y,ent_p,ent,T9_p,T9,T9h_p,T9h,stepsize)
         end if

         ! mass conservation check
         call masscalc(Y, m_tot)

         ! periodic output
         call output_nr_diagnostic(cnt,k,nr_count,time,T9,stepsize,m_tot,Y_p,Y)

         ! exit if mass is conserved sufficiently well
         if (((.not. heating_switch) .and. (nr_count.ge.nr_mincount) .and. (dabs(m_tot-1d0).lt.nr_tol)) .or. &
             ((heating_switch) .and. (nr_count.ge.nr_mincount) .and. (dabs(m_tot-1d0).lt.nr_tol)  .and. &
              ((max(T9_nrlast,T9)/min(T9_nrlast,T9)-1 .lt. timestep_hydro_factor))) .or. &
              ((heating_switch) .and. (nr_count.ge.nr_mincount) .and. (dabs(m_tot-1d0).lt.nr_tol) .and. &
               k .eq. adapt_stepsize_maxcount)) &
              exit newton_raphson


         T9_nrlast = T9
      end do newton_raphson


      ! See if one has to repeat the timestep
      exit_condition = exit_condition .and. (nr_count.le.nr_maxcount)
      if (heating_switch) then
        ! Dont be too strict in the last adapt stepsize iteration
        if (k .ne. adapt_stepsize_maxcount) then
         exit_condition = exit_condition .and. (max(T9,T9_p)/min(T9,T9_p)-1 .lt. timestep_hydro_factor)
        end if
      end if
      if (exit_condition) exit adapt_stepsize

      ! Check if the maxcount is reached and try one last time
      ! to get it to convergence
      if (k .eq. adapt_stepsize_maxcount) then
        stepsize = max(initial_stepsize,1e-3*stepsize)
      else
        stepsize = 0.5*stepsize
      end if
   end do adapt_stepsize

   ! Update entropy
   if (.not. (heating_switch)) then
      if (trim(adjustl(trajectory_mode)) .eq. "from_file") then
         if (time.le.ztime(zsteps) .and. (T9 .gt. freeze_rate_temp*1.0000001)) then
            ! update entropy using the temperature
            state%abar = sum(Y(1:net_size)*isotope(1:net_size)%mass) &
                       / sum(Y(1:net_size))
            call timmes_eos(ink,T9*1.d9,rhob,Ye,state,eos_status)
            if(eos_status.ne.0) call raise_exception("An error occured in the EOS.",&
                                                     'advance_implicit_euler',430009)
            ent = state%s
            P   = state%p
         else
            ! adiabatic expansion
            ent = ent_p
            ! update pressure
            if (T9 .gt. freeze_rate_temp*1.0000001) then
                ! update entropy using the temperature
                state%abar = sum(Y(1:net_size)*isotope(1:net_size)%mass) &
                           / sum(Y(1:net_size))
                call timmes_eos(ink,T9*1.d9,rhob,Ye,state,eos_status)
                if(eos_status.ne.0) call raise_exception("An error occured in the EOS.",&
                                                          'advance_gear',430009)
                P   = state%p
              end if

         end if
      else if (trim(adjustl(trajectory_mode)) .eq. "analytic") then
         if (T9 .gt. freeze_rate_temp*1.0000001) then
           ! update entropy using the temperature
           state%abar = sum(Y(1:net_size)*isotope(1:net_size)%mass) &
                      / sum(Y(1:net_size))
           call timmes_eos(ink,T9*1.d9,rhob,Ye,state,eos_status)
           if(eos_status.ne.0) call raise_exception("An error occured in the EOS.",&
                                                     'advance_implicit_euler',430009)
           ent = state%s
           P   = state%p
         end if
      end if
   end if


   ! if still not converged, complain and exit
   if (k.gt.adapt_stepsize_maxcount) then
      call raise_exception("Max. number of steps reached. Probably try to change "//NEW_LINE('A')//&
                           'the "adapt_stepsize_maxcount", "nr_maxcount", or "nr_tol" parameter.',&
                           'advance_implicit_euler',430010)
   end if

   ! Update NSE abundances, only in NSE evolution mode
   if (evolution_mode .eq. EM_NSE) then
      ! Update at all?
      if (nse_calc_every .ne. 0) then
         ! Update every nse_calc_step
         if (modulo(cnt,nse_calc_every) .eq. 0) then
            call winnse_guess( T9, rhob, Ye, Y(ineu), Y(ipro), Y)
         end if
      end if
   end if

end subroutine advance_implicit_euler


!>
!! Advance system to the next step using the Gear method
!! (using \ref gear_module)
!!
!! This subroutine is similar to \ref advance_implicit_euler. However, no
!! additional check for convergence is done after the newton-raphson. This is
!! not necessary as the Gear solver has its estimate of a possible error
!! in the next step, given by evolving higher orders and comparing to lower ones.
!!
!! @see \ref gear_module, [Martin, D. 2017](https://tuprints.ulb.tu-darmstadt.de/6301/)
!!
!! \b Edited:
!!       -  DM 15.06.2017
!!       -  OK 27.08.2017
!!       -  MR 22.12.2020
!! .
subroutine advance_gear(cnt)
use parameter_class, only: heating_mode,nse_calc_every, &
                           trajectory_mode, gear_nr_eps, gear_escale, &
                           gear_nr_maxcount, freeze_rate_temp, &
                           initial_stepsize, adapt_stepsize_maxcount,&
                           timestep_hydro_factor, gear_nr_mincount, &
                           gear_ignore_adapt_stepsize
use global_class,    only: net_size, isotope,ipro,ineu, heating_switch
use hydro_trajectory,only: zsteps,ztime
use analysis,        only: output_nr_diagnostic
use winnse_module,   only: pf,winnse_guess
use nucstuff_class,  only: el_ab,masscalc
use jacobian_class,  only: jacobi_init
use pardiso_class,   only: netsolve
use nuclear_heating, only: nuclear_heating_update
use ls_timmes_eos_module, only: timmes_eos_state,timmes_eos,ink
use single_zone_vars
use gear_module, only: get_time, get_timestep, get_l1, &
                       get_solution, get_predictor_Y, get_predictor_dYdt, &
                       shift_tau, nordsieck_product, set_xi, set_l, &
                       set_new_result, prepare_next_step, gear_nr_count, &
                       get_solution_at, set_timestep, revert_timestep
implicit none
integer, intent(in)    :: cnt !< global iteration counter
!
integer        :: k           !< Dummy for nr_diagnostic
real(r_kind)   :: m_tot       !< Mass conservation
!
real(r_kind),dimension(net_size) :: predictor
real(r_kind),dimension(net_size) :: delta
real(r_kind),dimension(net_size) :: ydiff
real(r_kind),dimension(net_size) :: olol

type(timmes_eos_state) :: state
integer                :: eos_status
logical                :: exit_condition!< Helper variable that keeps
                                        !  track of the exit condition of the NR
real(r_kind)           :: T9_nrlast     !< Last temperature in NR

   adapt_stepsize: do k=0,adapt_stepsize_maxcount
        ! compute some handy reductions and auxiliary variables
        stepsize = get_timestep()

        olol = get_solution()
        Y = Y_p
        Rkm     = Rkm_p
        time    = time_p + stepsize
        ent     = ent_p
        P       = P_p

        call el_ab (Y, Ye)

        ! Get temperature, density, and radius
        call update_hydro(time, T9, rhob, Rkm, ent, Ye, T9h)

        T9_nrlast = T9

        ! shift tau list and insert h, set time
        call shift_tau
        time = get_time()

        ! predictor
        call nordsieck_product

        ! calculate xi vector
        call set_xi

        ! calculate l vector
        call set_l

        predictor = get_predictor_Y()
        Y = predictor

        !-- Newton-Rhapson iterative solve ---------
        newton_raphson: do gear_nr_count=1,gear_nr_maxcount
            ! Set a variable that takes track of problems
            exit_condition = .True.

            ! calculate jacobian and rhs of system of DE (f())
            call jacobi_init(time,T9,rhob,Rkm,Y,Y_p,dYdt,f,stepsize,evolution_mode)

            ! solve the system
            call netsolve(f, delta)

            ! add correction to the abundance array
            Y= Y + delta

            ! thermodynamics update
            call el_ab(Y, Ye)

            ! Check if thinks are reasonable or
            ! went terribly wrong
            if (Ye .lt. 0 .or. Ye .gt. 1) then
                exit_condition = .False.
                exit newton_raphson
            end if


            if (heating_switch) then
                ! Update NSE if heating is turned on
                if (evolution_mode .eq. EM_NSE) then
                    call winnse_guess( T9, rhob, Ye, Y(ineu), Y(ipro), Y)
                end if
                call nuclear_heating_update(gear_nr_count,rhob,Ye,pf,Y_p,Y,ent_p,ent,T9_p,T9,T9h_p,T9h,stepsize)
            end if

            ! mass conservation check
            call masscalc(Y, m_tot)

            ! periodic output
            call output_nr_diagnostic(cnt,k,gear_nr_count,time,T9,stepsize,m_tot,Y_p,Y)

            ! exit if mass is conserved sufficiently well
            ! Also check the temperature for the heating case
            if (((.not. heating_switch) .and. maxval(dabs(delta)/max(Y,gear_escale)).lt.gear_nr_eps .and. &
                dabs(m_tot-1.0d0).lt.gear_nr_eps .and. gear_nr_count .ge. gear_nr_mincount) .or. &
                (heating_switch .and. maxval(dabs(delta)/max(Y,gear_escale)).lt.gear_nr_eps .and. &
                dabs(m_tot-1.0d0).lt.gear_nr_eps .and. gear_nr_count .ge. gear_nr_mincount .and.&
                ((max(T9_nrlast,T9)/min(T9_nrlast,T9)-1 .lt. timestep_hydro_factor))) .or. &
                 ((heating_switch) .and. maxval(dabs(delta)/max(Y,gear_escale)).lt.gear_nr_eps .and. &
                 dabs(m_tot-1.0d0).lt.gear_nr_eps .and. &
                 (k .eq. adapt_stepsize_maxcount))) then
                exit newton_raphson
            end if

            T9_nrlast = T9

        end do newton_raphson

        ! Gear usually converges quite well without adapting the time step
        if (.not. gear_ignore_adapt_stepsize) then
            ! if previous loop didn't have problems converging ...
            ! This basically never happens, but it's good to have it
            ! Dont do it in NSE, because it will always fail
            exit_condition = exit_condition .and. ((gear_nr_count.le.gear_nr_maxcount) .or. (evolution_mode .eq. EM_NSE))
            if (heating_switch) then
            ! Dont be too strict in the last adapt stepsize iteration
                if (k .ne. adapt_stepsize_maxcount) then
                    exit_condition = exit_condition .and. (max(T9,T9_p)/min(T9,T9_p)-1 .lt. timestep_hydro_factor)
                end if
            end if
            if (exit_condition) exit adapt_stepsize
        else
            exit adapt_stepsize
        end if

        ! Check if the maxcount is reached and try one last time
        ! to get it to convergence
        if (k .eq. adapt_stepsize_maxcount) then
          stepsize = max(initial_stepsize,1e-3*stepsize)
        else
          stepsize = 0.5*stepsize
        end if
        call revert_timestep(stepsize)

    end do adapt_stepsize

    ! if still not converged, complain and exit
    if (k.gt.adapt_stepsize_maxcount) then
        call raise_exception("Max. number of steps reached. Probably try to change "//NEW_LINE('A')//&
                             'the "adapt_stepsize_maxcount", "gear_nr_maxcount", or "gear_nr_eps" parameter.',&
                             'advance_gear',430010)
     end if


   ! Update entropy
   if (.not. (heating_switch)) then
      if (trim(adjustl(trajectory_mode)) .eq. "from_file") then
         if (time.le.ztime(zsteps) .and. (T9 .gt. freeze_rate_temp*1.0000001)) then
            ! update entropy using the temperature
            state%abar = sum(Y(1:net_size)*isotope(1:net_size)%mass) &
                       / sum(Y(1:net_size))
            call timmes_eos(ink,T9*1.d9,rhob,Ye,state,eos_status)
            if(eos_status.ne.0) call raise_exception("An error occured in the EOS.",&
                                                     'advance_gear',430009)
            ent = state%s
            P   = state%p
         else
            ! adiabatic expansion
            ent = ent_p

            if (T9 .gt. freeze_rate_temp*1.0000001) then
                ! update entropy using the temperature
                state%abar = sum(Y(1:net_size)*isotope(1:net_size)%mass) &
                           / sum(Y(1:net_size))
                call timmes_eos(ink,T9*1.d9,rhob,Ye,state,eos_status)
                if(eos_status.ne.0) call raise_exception("An error occured in the EOS.",&
                                                          'advance_gear',430009)
                ! update pressure
                P   = state%p
              end if
         end if
      else if (trim(adjustl(trajectory_mode)) .eq. "analytic") then
         if (T9 .gt. freeze_rate_temp*1.0000001) then
           ! update entropy using the temperature
           state%abar = sum(Y(1:net_size)*isotope(1:net_size)%mass) &
                      / sum(Y(1:net_size))
           call timmes_eos(ink,T9*1.d9,rhob,Ye,state,eos_status)
           if(eos_status.ne.0) call raise_exception("An error occured in the EOS.",&
                                                     'advance_gear',430009)
           ent = state%s
           P   = state%p
         end if
      end if
   end if

   ! Update NSE abundances, only in NSE evolution mode
   if (evolution_mode .eq. EM_NSE) then
      ! Update at all?
      if (nse_calc_every .ne. 0) then
         ! Update every nse_calc_step
         if (modulo(cnt,nse_calc_every) .eq. 0) then
            call winnse_guess( T9, rhob, Ye, Y(ineu), Y(ipro), Y)
         end if
      end if
   end if

   ydiff = Y-predictor

   ! set the new result based on ydiff obtained in NR scheme
   call set_new_result(ydiff)

   ! calculate errors and possible stepsizes and take q(+-1) with largest stepsize
   call prepare_next_step

end subroutine advance_gear

end module timestep_module
