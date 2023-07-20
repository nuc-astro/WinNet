!> @file expansion_module.f90
!!
!! The error file code for this file is ***W17***.
!! @brief Module \ref expansion_module for parametric expansions
!!

!>
!! Subroutines to handle parametric evolution of hydrodynamic quantities after the
!! final timestep of the trajectory
!!
!! \b Edited:
!!         - 01.11.16
!!         - M.R. 26.12.20  - Got rid of expansion parameter
!! .
#include "macros.h"
module expansion_module
  use parameter_class,only:heating_mode
  use ls_timmes_eos_module
  use error_msg_class, only: raise_exception,int_to_str,&
                             num_to_str,write_data_to_std_out
  implicit none

save
  real(r_kind)           :: mdot     !< 1/4pi*dM/dr = rho*r^2 = const
  real(r_kind)           :: rho_0    !< Density [g/ccm] at the end of the trajectory
  real(r_kind)           :: rad_0    !< Radius [km] at the end of the trajectory
  real(r_kind)           :: t9_0     !< Temperature [GK] at the end of the trajectory
  real(r_kind)           :: t_0      !< Time at the end of the trajectory
  real(r_kind)           :: tau
  real(r_kind)           :: S_0       !< Entropy at the end of the trajectory (deprecated)
  type(timmes_eos_state) :: state     !< EOS state
  real(r_kind)           :: vel       !< velocity [km/s]
  real(r_kind)           :: rad       !< radial distance [km]
  logical                :: expand    !< flag if expansion will be present
  integer                :: expand_count !< 1 for the first timestep 2 for following steps, 0 for no expansion yet
  !-- Debugging
  integer, private       :: debug_exp !< file ID for debugging expansion



contains



!>
!! Initialize the expansion module and open files for debugging.
!!
!! This subroutine also determines if the simulation
!! has to use the expansion routine or not (depending
!! on parameter_class::termination_criterion).
!!
!! @author: Moritz Reichert
!!
!! \b Edited:
!!         - 28.02.17
!!         - MR 22.12.20
!!         - MR 22.01.21
!! .
  subroutine expansion_init()
    use file_handling_class
    use parameter_class, only: trajectory_mode,termination_criterion,final_time,&
                               final_temp,final_dens,expansiontype
    use hydro_trajectory
    implicit none


    expand = .False.

    if ((trim(adjustl(trajectory_mode)) .eq. 'from_file')) then
      ! Check if it has to expand
      if          (((termination_criterion .eq. 1) .and. (ztime(zsteps) .lt. final_time)) &
            & .or. ((termination_criterion .eq. 2) .and. (minval(dexp(ztemp(init_index:))) .gt. final_temp)) &
            & .or. ((termination_criterion .eq. 3) .and. (minval(dexp(zdens(init_index:))) .gt. final_dens)) &
            & .or. ((termination_criterion .eq. 4) )) then
         expand = .True.
      end if

      if (expand) then
         if (VERBOSE_LEVEL .gt. 2) then
            debug_exp = open_outfile('debug_expansion.dat')
            write(debug_exp,'(A)') '# File for debugging subroutine expansion in expansion_module.f90'
            write(debug_exp,'(*(A,4x))') '# Time [s]','Density [g/ccm]','Radius [km]','Old Radius [km]', &
            'Velocity [km/s]', 'Temperature [GK]'
         end if

         ! Consistency check
         if ((expansiontype .lt. 1) .or. (expansiontype .gt. 8)) then
            call raise_exception("Expansiontype ("//trim(adjustl(int_to_str(expansiontype)))&
                              //") not implemented yet. "//&
                              "Set it to a supported value!",&
                              "expansion_init",&
                              170003)
         end if

         ! Get velocity, final density and final radius
         call residual_expansion
         expand_count = 0
       end if
    end if
  end subroutine expansion_init



!>
!! Calculates the velocity for trajectory extrapolation by linear regression
!! on a few last points; the number of points is controlled by the parameter
!! parameter_class::extrapolation_width.
!!
!! @remark At verbose Level >2 this creates a file called "debug_vel.dat"
!!         which lists extrapolated points and the computed value of velocity.
!!
!! \b Edited:
!!       - 09.03.17
!!       - OK 13.01.2014
!!       - OK 17.06.2017 - cleanup, renamed: variance -> residual_expansion
!!       - MR 22.01.2021 - Changed entropy calculation
!! .
subroutine residual_expansion
   use parameter_class, only: extrapolation_width
   use file_handling_class, only: open_outfile
   use hydro_trajectory
   implicit none
   ! integer :: eos_status

   real(r_kind) :: v  ! velocity
   !
   integer :: i, debug_vel
   real(r_kind) :: er,ert,et,ets,det

   ! sanity check
   if (extrapolation_width.gt.zsteps) then
      call raise_exception('The parameter "extrapolation_width" ('//trim(adjustl(int_to_str(extrapolation_width)))//&
                           ') should be smaller'//NEW_LINE("A")//&
                           'than the amount of trajectory points (zsteps = '//trim(adjustl(int_to_str(zsteps)))//").",&
                           'residual_expansion',&
                           170004)
   endif
   if (extrapolation_width.lt.2) then
      call raise_exception('The parameter "extrapolation_width" ('//trim(adjustl(int_to_str(extrapolation_width)))//&
                           ') should be larger than 1.',&
                           'residual_expansion',&
                           170005)
   endif

   if (extrapolation_width.eq.2) then
      ! handle case of two endpoints
      v= (exp(zrad(zsteps)) - exp(zrad(zsteps-1))) &
           / (ztime(zsteps) - ztime(zsteps-1))
   else
      ! more than two points
      et  = 0.0
      ets = 0.0
      er  = 0.0
      ert = 0.0
      do i=zsteps-extrapolation_width+1,zsteps
         et  = et  + ztime(i)
         ets = ets + ztime(i)**2
         er  = er  + exp(zrad(i))
         ert = ert + ztime(i)*exp(zrad(i))
      enddo
      det = et*et - extrapolation_width*ets
      v= 0.0
      if (abs(det).gt.1e-20) then
         v = (er*et - extrapolation_width*ert)/det
      endif
   endif

   if (VERBOSE_LEVEL.ge.1) then
      call write_data_to_std_out("Expansion velocity:",num_to_str(v),"[km/s]")
      ! print '(A,ES12.5,A)', &
      !    " Using expansion with residual velocity: ", v," [km/s]"
   endif

   ! initialize module parameters
   vel = v
   rad_0 = exp(zrad(zsteps))
   rho_0 = exp(zdens(zsteps))
   T9_0 = exp(ztemp(zsteps))
   t_0 = ztime(zsteps)
   mdot = rho_0 * rad_0**2

   ! MR: The difference between the EOS entropy and
   !     T9^3/rho (taking some dummy ye and abar)
   !     differs by 5 magnitudes. This seems wrong.
   !     In my understanding the entropy is only proportional
   !     to T9^3/rho and not equal
   !     (see Python script of expansion test).
   !     The entropy S_0 is, however, also not used so far.
   !     Therefore, I commented it out.

   ! calculate entropy
   ! if (.not. include_nuclear_heating) then
      ! state%abar =20
      ! call timmes_eos(ink,T9_0*1.d9,rho_0,0.5,state,eos_status) ! Careful, dummy abar and ye used!
      ! if(eos_status.ne.0) stop 'error in EOS'
      ! S_0 = state%s
      ! write(*,*) "EOS: S=",S_0
      ! S_0 = T9_0**3/rho_0
      ! write(*,*) "T9**3/rho",S_0
   ! endif

   ! Debugging
   if (VERBOSE_LEVEL .gt. 2) then
      ! open a file 'debug_vel.dat' and record extrapolated value
      debug_vel = open_outfile('debug_vel.dat')
      write(debug_vel,'(A)') '# Velocity extrapolation info'
      write(debug_vel,'(A)') '# 1:step 2:time[s] 3:R[km]'
      do i=zsteps-extrapolation_width+1,zsteps
         write (debug_vel,'(i4,2(1x,es15.7))') i, ztime(i), exp(zrad(i))
      enddo
      write (debug_vel,'("vel = ",es15.7, " km/s")') v
      close(debug_vel)
   end if

  if (vel .lt. 0.) call raise_exception("Expansion velocity was negative!",&
                                        'residual_expansion',&
                                        170006)

end subroutine residual_expansion

!>
!! Returns temperature, radius, density, and entropy after the trajectory has ended.
!!
!! This depends on the \ref parameter_class::expansiontype parameter.
!! Possible values are:
!! <table>
!! <caption id="multi_row">Expansiontype parameter</caption>
!! <tr><th> Value     <th> Expansion                          <th> Literature
!! <tr><td> 1         <td> Adiabatic expansion                <td> [Korobkin et al. 2012](https://ui.adsabs.harvard.edu/abs/2012MNRAS.426.1940K/abstract)
!! <tr><td> 2         <td> Exponential expansion              <td> -
!! <tr><td> 3         <td> Expansion with const. velocity     <td> -
!! <tr><td> 4         <td> Expansion with const. velocity     <td> [Fujimoto et al. 2008](https://ui.adsabs.harvard.edu/abs/2008ApJ...680.1350F/abstract)
!! </table>
!!
!! \b Edited:
!!     - 26.01.21 - MR
!!     - 17.02.21 - MR: Cleaned up the routine and removed some expansion types
!! .
subroutine expansion(t,dt,ye,dens,rad_in,temp,entropy,vel_in,abar)
    use ls_timmes_eos_module
    use parameter_class,  only: expansiontype, freeze_rate_temp
    use hydro_trajectory, only: ztemp,zsteps
    implicit none

    real(r_kind),intent(in)       :: t
    real(r_kind),intent(in)       :: dt
    real(r_kind),intent(in)       :: abar  ! Average mass number (TODO: add this as an use statement after restructure of the code)
    real(r_kind),intent(inout)    :: dens
    real(r_kind),intent(inout)    :: rad_in
    real(r_kind),intent(inout)    :: temp
    real(r_kind),intent(in)       :: entropy
    real(r_kind),intent(inout)    :: vel_in
    real(r_kind),intent(in)       :: ye
    real(r_kind)                  :: time
    real(r_kind)                  :: rold,dold
    logical                       :: calc_temperature
    integer                       :: eos_status

    ! It can happen that the function is called even when there should be
    ! no expansion (e.g., because the timestep is restricted after the
    ! hydro update call). Therefore return without doing something.
    if (.not. expand) return

    rold = rad_in
    dold = dens
    time = t+dt-t_0 ! time that has passed since the last step in the trajectory

    calc_temperature = ((heating_mode .eq. 0) .and. (temp .gt. freeze_rate_temp))

    ! calc_temperature = (((heating_mode .eq. 2) .and. (temp .gt. freeze_rate_temp)) &
    !                    .or. (temp .gt. min(max(dexp(ztemp(zsteps)),1d-6),freeze_rate_temp) &
    !                    .and. (.not. (heating_mode .eq. 1))))

    select case (expansiontype)
    case(1) !adiabatic expansion (from Korobkin et al. (2012))
       rad_in = rad_0 + time*vel_in
       dens = rho_0*(t_0/(time+t_0))**3
       ! Ensure to call the EOS only at valid points, i.e., above some temperature
       if (calc_temperature) then
          state%abar = abar
          call timmes_eos(ins,entropy,dens,ye,state,eos_status)
          if(eos_status.ne.0) call raise_exception('Error when calling the EOS.','expansion',&
                                                   170007)
          temp = 1.d-9*state%t
       end if
    case(2) !exponential expansion
       rad_in = rad_0 + time*vel_in
       tau  = t_0 ! TODO: move to parameter_class
       dens = rho_0*dexp(-time/tau)
       ! Prevent underflows
       dens = max(dens,1d-99)
       ! Ensure to call the EOS only at valid points, i.e., above some temperature
       if (calc_temperature) then
          state%abar = abar
          call timmes_eos(ins,entropy,dens,ye,state,eos_status)
          if(eos_status.ne.0) call raise_exception('Error when calling the EOS.','expansion',&
                                                   170007)
          temp = 1.d-9*state%t
       end if
    case(3) !expansion with constant velocity
       rad_in  = rad_0 + time*vel_in
       dens = mdot/(rad_in**2)! rho_0 *(r_0/r)**2
       ! Ensure to call the EOS only at valid points, i.e., above some temperature
       if (calc_temperature) then
          state%abar = abar
          call timmes_eos(ins,entropy,dens,ye,state,eos_status)
          if(eos_status.ne.0) call raise_exception('Error when calling the EOS.','expansion',&
                                                   170007)
          temp = 1.d-9*state%t
       end if
    case(4) !expansion with constant velocity as in fujimoto apj 680
       rad_in  = rad_0 + time*vel_in
       ! Free spherical symmetric expansion
       dens = rho_0*(rad_0/rad_in)**3
       ! Ensure to call the EOS only at valid points, i.e., above some temperature
       if (calc_temperature) then
         state%abar = abar
         call timmes_eos(ins,entropy,dens,ye,state,eos_status)
         if(eos_status.ne.0) call raise_exception('Error when calling the EOS.','expansion',&
                                                  170007)
         temp = 1.d-9*state%t
       end if
       ! Alternatively:
       ! The temperature refers to an adiabatic expansion as S \propto T^3/rho = const
       ! temp = t9_0*rad_0/rad
    case default
       call raise_exception('Expansion type not (yet) supported.',&
                            'expansion',170008)
    end select

    ! Set the temperature above some threshold value,
    ! that is either the minimum of the trajectory,
    ! 1d-2GK or 1d-6GK.
    ! The rates are frozen in any case at 1d-2GK in the Jacobian class
    if (calc_temperature) then
        temp = max(temp,freeze_rate_temp)
    end if

    ! Debug file
    if (VERBOSE_LEVEL .gt. 2) then
        write(debug_exp,'(*(10es23.14,4x))')time,dens,rad_in,rold,vel_in,temp
    end if

    ! Check which expansion step this was
    if (expand_count .lt. 2) expand_count = expand_count+1

    return

end subroutine expansion


!>
!! Close debug files of the module
!!
!! @author: Moritz Reichert
!!
!! \b Edited:
!!    - 28.02.17
!! .
subroutine expansion_finalize()
 use file_handling_class
 implicit none

 if ((VERBOSE_LEVEL .gt. 2) .and. (expand)) then
   call close_io_file(debug_exp,'debug_expansion.dat')
 end if

end subroutine expansion_finalize

end module expansion_module
