!> @file analysis.f90
!!
!! The error file code for this file is ***W10***.
!! @brief Module \ref analysis with various analysis subroutines
!!

!> Contains various routines for analysis and diagnostic output
#include "macros.h"
module analysis

  use parameter_class, only: nrdiag_every,timescales_every
  use global_class,    only: net_size, ihe4, ineu, ipro
  use file_handling_class

  implicit none

  real(r_kind),dimension(:),allocatable :: Sn                                    !< array of n-separation energies
  real(r_kind)                          :: Sn_ave                                !< average n-separation energy
  integer                               :: mainout_unit                          !< file descriptor for the main analysis output file
  integer                               :: nrdiag_unit                           !< diagnostic output inside the Newton-Rhapson loop
  integer                               :: track_unit                            !< file id for tracked nuclei
  integer                               :: nu_loss_gain_unit                     !< file id for neutrino loss/gain
  integer                               :: cl_start, cl_end, cl_rate, cl_cmax    !< system clock variables
  integer,private                       :: tsfile                                !< timescale file unit
  integer,private                       :: nuc_heat_file                         !< nuclear heating file unit
  real(r_kind)                          :: tsnp, T9snp, rosnp                    !< for snapshot printing triggers
  integer                               :: snapcnt                               !< snapshot count
  integer                               :: flowcnt                               !< flow printing counter
  real,dimension(:),allocatable         :: snapshot_time                         !< point in time for custom snapshot
  integer                               :: snapshot_amount                       !< Amount of custom snapshots

  !
  ! Public and private fields and methods of the module
  !
  public:: &
     analysis_init, output_initial_step, output_final_step, output_iteration, &
     analysis_finalize, output_nr_diagnostic, output_final_stats, finab

  private:: &
     output_track_nuclei, output_mainout, output_snapshot, printsnap, &
     calc_nseparation_energy, calc_nuclear_heating, calc_av_timescales, &
     finab_sort, output_nu_loss

contains

!>
!! Open files, write headers, allocate space etc.
!!
subroutine analysis_init()
  use parameter_class
  use global_class,         only: net_names,track_nuclei_nr,track_nuclei_indices
  use benam_class,          only: findaz
  use ls_timmes_eos_module
#ifdef USE_HDF5
  use hdf5_module,          only: init_hdf5_module
#endif
  implicit none

  ! neutron separation energies
  real(r_kind)  :: val                        !< aux variable for reading Sn file
  integer       :: ai, ni, zi, snfile, fstat  !< temporary variables needed for reading the Sn file
  character     :: nl, buf(30)                !< newline character, buffer
  integer       :: i,k
  character*10,dimension(track_nuclei_nr) :: tmp_nucleinames !< Name of the nuclei that will be tracked

  INFO_ENTRY("analysis_init")

#ifdef USE_HDF5
   ! Only initialize if it is really possible to execute it
   ! (i.e., correct compiler h5fc)
   call init_hdf5_module()
#endif

  cl_cmax= 100000
  snapcnt= 1
  flowcnt= 1
  nl = NEW_LINE('A')

  ! main output file: create, write header
  if (mainout_every .gt. 0) then
     mainout_unit= open_outfile("mainout.dat")
     write(mainout_unit,'(4A)') &
            '# 1:iteration 2:time[s], 3:T[GK], 4:rho[g/cm3], 5:Ye '//nl &
          , '# 6:R[km], 7:Y_n, 8:Y_p, 9:Y_alpha, 10:Y_lights '//nl  &
          , '# 11:Y_heavies 12:<Z> 13:<A> 14:entropy [kB/baryon]'//nl &
          , '# 15: Pressure [dyn/cm^2] (16:Sn [MeV])'
  endif

  ! read neutron separation energies
  if (calc_nsep_energy) then
     snfile= open_outfile(nsep_energies_file)

     allocate(Sn(net_size))
     Sn(1:net_size) = 0d0
     read(snfile, *) buf ! skip the first line
     do
        read(snfile, *, iostat=fstat) Zi, Ni, Ai, val
        if (fstat /= 0) exit
        k = findaz(ai, zi)
        if (k>0) then
           Sn(k) = val
        endif
     enddo
     call close_io_file(snfile,nsep_energies_file)
  endif

  ! Open the timescale file and write a header
  if (timescales_every .gt. 0) then
    tsfile = open_outfile('timescales.dat')
    write(tsfile, '((A1,(A22,3x),*(A23,3x)))') &
       "#","time[s]","Temperature [GK]","tau_ga [s]","tau_ag [s]","tau_ng [s]","tau_gn [s]",   &
       "tau_pg [s]", "tau_gp [s]", "tau_np [s]", "tau_pn [s]", &
       "tau_an [s]", "tau_na [s]", "tau_ap [s]", "tau_pa [s]", "tau_beta [s]", "tau_alpha [s]", &
       "tau_nfiss [s]", "tau_sfiss [s]", "tau_bfiss [s]"
  end if

  ! Open the track_nuclei file and write a header
  if (track_nuclei_every .gt. 0) then
     ! Convert indices to names
    do i=1, track_nuclei_nr
      tmp_nucleinames(i) = "Y( "//net_names(track_nuclei_indices(i))//" )"
    end do
    ! Write the header
    track_unit = open_outfile('tracked_nuclei.dat')
    write(track_unit, '((A,20X,*(A,18X)))') "# time[s]",tmp_nucleinames(:)
  end if


  if (nrdiag_every.gt.0) then
     nrdiag_unit= open_outfile("nrdiag.dat")
     write(nrdiag_unit,'(11A)') &
            '# Newton-Raphson diagnostic output'//nl               &
          , '# 1 : global iteration count (cnt)'//nl               &
          , '# 2 : global time [s]'//nl                            &
          , '# 3 : temperature T9 [GK]'//nl                        &
          , '# 4 : adapt stepsize loop counter (k)'//nl            &
          , '# 5 : NR loop counter (nr_count)'//nl                 &
          , '# 6 : global timestep [s]'//nl                        &
          , '# 7 : mass conservation error'//nl                    &
          , '# 8 : maximal abundance change (eps)'//nl             &
          , '# 9 : most rapidly evolving isotope (epsl)'//nl &
          , '# 10: abundance of the isotope epsl, y(epsl)'
  end if

  ! Open the energy generation file and write a header
  if (engen_every .gt. 0) then
    nuc_heat_file = open_outfile('generated_energy.dat')
    ! Header
    write(nuc_heat_file, '((A,5X))') "# Generated Energy for each nuclear reaction type given in erg/g/s"
    write(nuc_heat_file, '(17(A23,3X))') &
       "#               time[s]", "Engen(Total)", "S_src","Engen(weak)","Engen(n,g)",&
       "Engen(p,g)","Engen(a,g)","Engen(n,p)","Engen(n,a)", &
       "Engen(p,a)", "Engen(fiss)"
  end if



  if ((nu_loss_every .gt. 0)) then
    nu_loss_gain_unit = open_outfile('nu_loss_gain.dat')
     write(nu_loss_gain_unit, '(A)') "# File containing the neutrino loss/gain for the nuclear heating."
     write(nu_loss_gain_unit, '(A)') "# Negative values mean that energy is added to the system."
     write(nu_loss_gain_unit, '(A)') "# Neutrino energies are given in Erg/g/s."
     write(nu_loss_gain_unit, '((A1,(A18,3x),*(A19,3x)))') &
        "#","time[s]","T[GK]", "rho[g/cm3]", "R[km]", "nu_total", "nu_beta",&
        "nu_nuheat", "nu_thermal"
  end if


  INFO_EXIT("analysis_init")

end subroutine analysis_init


!>
!! Output first step
!!
!! This routine calls output_iteration at the first step
!!
!! @see output_final_step
subroutine output_initial_step(ctime,T9,rho_b,entropy,pressure,Rkm,Y,pf)
  implicit none
  real(r_kind),intent(in) :: ctime    !< actual time [s]
  real(r_kind),intent(in) :: T9       !< initial temperature [GK]
  real(r_kind),intent(in) :: rho_b    !< initial density [gcc]
  real(r_kind),intent(in) :: entropy  !< entropy [kB/nucleon]
  real(r_kind),intent(in) :: pressure !< Pressure [dyn/cm**2]
  real(r_kind),intent(in) :: Rkm      !< Radius [km]
  real(r_kind),dimension(net_size),intent(in)   :: Y !< array of abundances (Y_i)
  real(r_kind),dimension(0:net_size),intent(inout) :: pf !< partition functions

  INFO_ENTRY('output_initial_step')
  call output_iteration(0,ctime,T9,rho_b,entropy,pressure,Rkm,Y,pf)
  INFO_EXIT('output_initial_step')

end subroutine output_initial_step


!>
!! Output final step
!!
!! This routine calls output_iteration at the last step
!!
!! @see output_final_step
!! \b Edited: OK 26.08.2017
subroutine output_final_step(cnt,ctime,T9,rho_b,entropy,pressure,Rkm,Y,pf)
  use parameter_class, only: out_every,snapshot_every,flow_every,mainout_every, &
                             track_nuclei_every,termination_criterion
  implicit none
  integer,intent(in)      :: cnt      !< current iteration counter
  real(r_kind),intent(in) :: ctime    !< actual time [s]
  real(r_kind),intent(in) :: T9       !< initial temperature [GK]
  real(r_kind),intent(in) :: rho_b    !< initial density [gcc]
  real(r_kind),intent(in) :: entropy  !< entropy [kB/nucleon]
  real(r_kind),intent(in) :: pressure !< Pressure [dyn/cm**2]
  real(r_kind),intent(in) :: Rkm      !< Radius [km]
  real(r_kind),dimension(net_size),intent(in)   :: Y  !< array of abundances (Y_i)
  real(r_kind),dimension(0:net_size),intent(inout) :: pf !< partition functions

  INFO_ENTRY('output_final_step')

  ! hack output frequency parameters to output the final step out-of-order
  if (out_every.gt.0)          out_every= 1
  if (mainout_every.gt.0)      mainout_every= 1
  if (snapshot_every.gt.0)     snapshot_every= 1
  if (flow_every.gt.0)         flow_every= 1
  if (timescales_every.gt.0)   timescales_every= 1
  if (track_nuclei_every.gt.0) track_nuclei_every= 1
  print '("---------------")'
  call output_iteration(cnt,ctime,T9,rho_b,entropy,pressure,Rkm,Y,pf)
  print '("============================")'
  select case(termination_criterion)
  case(0)
     print'(A)',"End of trajectory file reached."
  case(1)
     print'(A)',"Final time reached."
  case(2)
     print'(A)',"Final temperature reached."
  case(3)
     print'(A)',"Final density reached."
  case(4)
     print'(A)',"Freeze-out reached."
  endselect

  INFO_EXIT('output_final_step')

end subroutine output_final_step


!>
!! Periodic output of analysis data for current iteration
!!
!! Prints various outputs depending on the input parameter.
!!
!! \b Edited:
!!       - MR : 19.01.21
!! .
subroutine output_iteration(cnt,ctime,T9,rho_b,entropy,pressure,Rkm,Y,pf)
  use parameter_class,  only: out_every,snapshot_every,mainout_every, &
                              track_nuclei_every,heating_mode, &
                              custom_snapshots,flow_every,engen_every,&
                              top_engen_every,h_snapshot_every,h_mainout_every,&
                              h_custom_snapshots,h_track_nuclei_every,&
                              h_timescales_every,h_flow_every,h_engen_every,&
                              nu_loss_every, h_nu_loss_every, unit, use_thermal_nu_loss
  use nucstuff_class,   only: inter_partf
  use nuclear_heating,  only: output_nuclear_heating,output_top_contributor,&
                              qdot_nu, qdot_bet, qdot_th
  use thermal_neutrino_module, only: thermal_neutrinos
  use single_zone_vars, only: evolution_mode, stepsize
  use global_class,     only: isotope
  use flow_module,      only: flowcalc, flows
#ifdef USE_HDF5
  use hdf5_module
#endif
  implicit none

  integer,intent(in)      :: cnt      !< current iteration counter
  real(r_kind),intent(in) :: ctime    !< actual time [s]
  real(r_kind),intent(in) :: T9       !< temperature [GK]
  real(r_kind),intent(in) :: rho_b    !< density [gcc]
  real(r_kind),intent(in) :: entropy  !< entropy [kB/nucleon]
  real(r_kind),intent(in) :: pressure !< Pressure [dyn/cm**2]
  real(r_kind),intent(in) :: Rkm      !< length scale [km]
  real(r_kind),dimension(net_size),intent(in)      :: Y  !< array of abundances (Y_i)
  real(r_kind),dimension(0:net_size),intent(inout) :: pf !< partition functions
  integer                 :: i              !< Loop variable
  real(r_kind)            :: ysum,yasum,abar,ye

  INFO_ENTRY("output_iteration")

  ! screen output
  if (out_every.gt.0) then
     if (mod(cnt,50*out_every)==0) then
        print *
        print '(A)',"#-- 1:it 2:time[s] 3:T9[GK] 4:rho[g/cm3] &
                     5:entropy[kB/baryon]"
     endif
     if (mod(cnt,out_every) .eq. 0) then
        print '(I8,X,ES14.7,X,F10.7,X,ES14.7,X,F12.7)', &
              cnt,ctime,T9,rho_b,entropy
     endif

     if(mod(cnt,out_every).eq.0 .and. (heating_mode .eq. 1)) then
        if (VERBOSE_LEVEL.ge.2) then
            call output_nuclear_heating(cnt,ctime)
        end if
     endif
  endif

  ! This is only done if the correct compiler (h5fc) was used
#ifdef USE_HDF5
  ! Store snapshots in hdf5
  if (h_snapshot_every.gt.0) then
     if (mod(cnt,h_snapshot_every).eq. 0) then
        ! Ye = sum(isotope(1:net_size)%p_nr*Y(1:net_size))
        call extend_snaps(ctime,Y)
     end if
  end if

  ! Store mainout in hdf5
  if (h_mainout_every.gt.0) then
     if (mod(cnt,h_mainout_every).eq. 0) then
        call extend_mainout(cnt,ctime,T9,rho_b,entropy,pressure,Rkm,Y)
     end if
  end if

  ! Snapshot output for specific times
  if (h_custom_snapshots) then
     ! Test if a snapshot should be done
     h_make_snap: do i=1,snapshot_amount
        ! Time is within certain tolerance? Do the snapshot!
        if (abs(ctime-snapshot_time(i)) .le. 1e-20) then
           call extend_snaps(ctime,Y)
           exit h_make_snap
        end if
     enddo h_make_snap
  end if

  ! Output the tracked nuclei
  if (h_track_nuclei_every .gt. 0) then
      if (mod(cnt,h_track_nuclei_every).eq.0) then
         call extend_track_nuclei(ctime,Y)
      end if
  end if

  ! Output average timescales
  if (h_timescales_every .gt. 0) then
     if ((mod(cnt,h_timescales_every).eq.0) .and. (evolution_mode .ne. EM_NSE)) then
        ! calculate averaged timescales of various types of reactions
        ! calculate partition function first
        call inter_partf (T9, pf)
        call calc_av_timescales(ctime,T9,Y,.True.)
     end if
  end if

  ! Output the flow
  if ((h_flow_every .gt. 0) .and. (cnt .ne. 0)) then
     if ((mod(cnt,h_flow_every).eq.0)) then
        call flowcalc(Y)
        call extend_flow(ctime,T9,rho_b,flows,size(flows),Y)
     end if
  end if

  ! Calculate and output generated energy
  if (h_engen_every .gt. 0) then
     if (mod(cnt,h_engen_every).eq.0) then
        call calc_nuclear_heating(.True.)
     end if
  end if

  ! Calculate and output generated energy
  if ((h_nu_loss_every .gt. 0)) then
     if (mod(cnt,h_nu_loss_every).eq.0) then

        if ((heating_mode .eq. 0) .and. (use_thermal_nu_loss)) then
            ysum= sum(Y(1:net_size))
            yasum= sum(Y(1:net_size)*isotope(1:net_size)%mass)
            Ye = sum(isotope(1:net_size)%p_nr*Y(1:net_size))
            abar= yasum / ysum
            call thermal_neutrinos(abar, Ye, T9, rho_b, stepsize, qdot_th)
            qdot_th = qdot_th * unit%hix ! erg/g -> MeV/baryon
        else
            qdot_th = 0d0
        end if

        call extend_nu_loss_gain(ctime,T9,rho_b,rkm,qdot_th/unit%hix/stepsize,&
                                                    qdot_nu/unit%hix/stepsize,&
                                                    qdot_bet/unit%hix/stepsize)
     end if
  end if

#endif

  ! Output list of isotopes that dominantly contribute to the energy generation
  if ((top_engen_every .ne. 0) .and. (heating_mode .gt. 0)) then
    if(mod(cnt,top_engen_every).eq.0) call output_top_contributor(cnt,ctime,T9,rho_b)
  endif

  ! main data log
  if (mainout_every .gt.0) then
     if (mod(cnt,mainout_every).eq.0) then
        call output_mainout(cnt,ctime,T9,rho_b,entropy,pressure,Rkm,Y)
     end if
  end if

  ! Iterative Snapshot and flow output
  if(snapshot_every.gt.0) then
     if(mod(cnt,snapshot_every).eq.0) then
        call output_snapshot(ctime,T9,rho_b,Y)
     end if
  end if

  if(flow_every.gt.0) then
     if(mod(cnt,flow_every).eq.0) then
        call output_flow(ctime,T9,rho_b,Y)
     end if
  end if

  if ((nu_loss_every .gt. 0)) then
     if (mod(cnt,nu_loss_every).eq.0) then

        if ((heating_mode .eq. 0) .and. (use_thermal_nu_loss)) then
            ysum= sum(Y(1:net_size))
            yasum= sum(Y(1:net_size)*isotope(1:net_size)%mass)
            Ye = sum(isotope(1:net_size)%p_nr*Y(1:net_size))
            abar= yasum / ysum
            call thermal_neutrinos(abar, Ye, T9, rho_b, stepsize, qdot_th)
            qdot_th = qdot_th * unit%hix ! erg/g -> MeV/baryon
        else
            qdot_th = 0d0
        end if

        call output_nu_loss(ctime,T9,rho_b,Rkm)
     end if
  end if


  ! Snapshot output for specific times
  if (custom_snapshots) then
     ! Test if a snapshot should be done
     make_snap: do i=1,snapshot_amount
        ! Time is within certain tolerance? Do the snapshot!
        if (abs(ctime-snapshot_time(i)) .le. 1e-20) then
           call output_snapshot(ctime,T9,rho_b,Y)
           exit make_snap
        end if
     enddo make_snap
  end if

  if (timescales_every .gt. 0) then
     if ((mod(cnt,timescales_every).eq.0) .and. (evolution_mode .ne. EM_NSE)) then
        ! calculate averaged timescales of various types of reactions
        ! calculate partition function first
        call inter_partf (T9, pf)
        call calc_av_timescales(ctime,T9,Y)
     end if
  end if

 ! Output the tracked nuclei
  if (track_nuclei_every .gt. 0) then
     if (mod(cnt,track_nuclei_every).eq.0) then
        call output_track_nuclei(ctime,Y)
     end if
  end if

  ! Calculate generated energy
  if (engen_every .gt. 0) then
     if (mod(cnt,engen_every).eq.0) then
        call calc_nuclear_heating(.False.)
     end if
  end if

  INFO_EXIT("output_iteration")
end subroutine output_iteration


!> Output the abundances of specific nuclei
!!
!! First column is always the time, followed by the abundances of
!! nuclei given in "track_nuclei_file". An example of the file looks like:
!! \l_file{
!! # time[s]                    Y(     n )                  Y(   o24 )                  Y(   f24 )
!!  0.0000000000000000E+00      0.0000000000000000E+00      2.0833333333333332E-02      0.0000000000000000E+00
!!  1.9999999999999999E-20      8.5150206007805894E-19      2.0833333333333332E-02      1.8661570495808601E-21
!!  5.9999999999999994E-20      2.5545061802341771E-18      2.0833333333333332E-02      5.5984711487425806E-21
!!  1.3999999999999998E-19      5.9605144205464133E-18      2.0833333333333332E-02      1.3063099347066021E-20
!! ...}
!!
!! @author: Moritz Reichert
!!
!! \b Edited:
!!      - 09.06.17
!! .
subroutine output_track_nuclei(ctime,Y)
   use global_class,       only: track_nuclei_indices,track_nuclei_nr
   implicit none
   real(r_kind),intent(in)                          :: ctime    !< Actual time
   real(r_kind),dimension(net_size),intent(in)      :: Y        !< array of abundances (Y_i)
   real(r_kind),dimension(track_nuclei_nr)          :: tmp_Y    !< Stores the abundance of the tracked nuclei
   integer                                          :: i        !< Loop variable

   INFO_ENTRY("output_track_nuclei")
   ! Store the abundance temporary
   do i=1, track_nuclei_nr
      tmp_Y(i) = Y(track_nuclei_indices(i))
      if (tmp_Y(i) .lt. 1d-99) tmp_Y(i) = 0.
   end do
   ! Output
   write(track_unit,'(*(es23.16,5x))') ctime,tmp_Y(:)
   INFO_EXIT("output_track_nuclei")

end subroutine output_track_nuclei


!> Output the energy that enters and leaves the system
!!
!! Outputs the energy that enters and leaves the system when
!! nuclear_heating is turned on.
!! This output is controlled by the parameter \ref nu_loss_every .
!! The output is written to the file "nu_loss_gain_file".
!! An example of the file could look like:
!! \file{
!! \# File containing the neutrino loss/gain for the nuclear heating.
!! \# Negative values mean that energy is added to the system.
!! \# Neutrino energies are given in MeV/baryon/s.
!! \#           time[s]                 T[GK]            rho[g/cm3]                 R[km]              nu_total  ...
!!         1.10385E+00           1.00000E+01           4.89334E+07           3.55626E+02           0.00000E+00  ...
!!         1.10385E+00           1.00000E+01           4.89334E+07           3.55626E+02          -5.05524E+00  ...
!!         1.10385E+00           1.00000E+01           4.89334E+07           3.55626E+02          -5.05524E+00  ...
!!         1.10385E+00           1.00000E+01           4.89334E+07           3.55626E+02          -5.05524E+00  ...
!! ...}
!! @author Moritz Reichert
!! @date   12.04.2023
subroutine output_nu_loss(ctime,temp,rho,Rkm)
    use single_zone_vars,only: stepsize
    use nuclear_heating, only: qdot_bet, qdot_nu, qdot_th
    use parameter_class, only: unit
    implicit none
    real(r_kind),intent(in) :: ctime !< actual time (s)
    real(r_kind),intent(in) :: temp  !< temperature (GK)
    real(r_kind),intent(in) :: rho   !< density (g/cm3)
    real(r_kind),intent(in) :: Rkm   !< radius (km)
    ! local variables
    real(r_kind) :: qdot_tot

    qdot_tot = qdot_bet + qdot_nu + qdot_th

    write(nu_loss_gain_unit, '(*(es19.5,3x))') ctime, temp, rho, Rkm, &
                                               qdot_tot/unit%hix/stepsize, &
                                               qdot_bet/unit%hix/stepsize, &
                                               qdot_nu/unit%hix/stepsize,  &
                                               qdot_th/unit%hix/stepsize

end subroutine output_nu_loss


!>
!! Output Mainout file
!!
!! This file contains all basic quantities such as neutron abundance Yn,
!! time, temperature, density, Abar, and much more. An example file
!! looks like:
!! \l_file{
!! # 1:iteration 2:time[s], 3:T[GK], 4:rho[g/cm3], 5:Ye
!! # 6:R[km], 7:Y_n, 8:Y_p, 9:Y_alpha, 10:Y_lights
!! # 11:Y_heavies 12:<Z> 13:<A> 14:entropy [kB/baryon] (15:Sn [MeV])
!!    0  0.0000000000000000E+00   1.0964999999999998E+01   8.7095999999999795E+12   3.4880000000001680E-02   4.9272999999999975E+01   8.7430471311642388E-01   1.8154340816658842E-15   1.8313910002116333E-13   8.4293473893924100E-09   1.6154109228777289E-03   3.9820982195819934E-02   1.1416565996507526E+00   9.4638587929828290E-03
!!   10  4.0940000000000001E-09   1.0965038071994945E+01   8.7093701186135752E+12   3.4880006806212963E-02   4.9273176232580205E+01   8.7430472419630612E-01   1.8157922458135375E-15   1.8317495847108614E-13   8.4302150496944692E-09   1.6154086944791908E-03   3.9820989563731306E-02   1.1416565881127740E+00   9.4639764072160116E-03
!! ...}
!!
!!
subroutine output_mainout(cnt,ctime,T9,rho_b,entropy,pressure,Rkm,Y)
  use parameter_class, only: calc_nsep_energy
  use global_class,    only: isotope
  implicit none

  integer,intent(in)      :: cnt                         !< current iteration counter
  real(r_kind),intent(in) :: ctime                       !< current time
  real(r_kind),intent(in) :: T9                          !< temperature [GK]
  real(r_kind),intent(in) :: rho_b                       !< density [gcc]
  real(r_kind),intent(in) :: entropy                     !< entropy [kB/nucleon]
  real(r_kind),intent(in) :: pressure                    !< pressure [dyn/cm**2]
  real(r_kind),intent(in) :: Rkm                         !< length scale [km]
  real(r_kind),dimension(net_size),intent(in)      :: Y  !< array of abundances (Y_i)
  !
  real(r_kind) :: Ye,Yneu,Ypro,Yhe4,ylight,yheavies,ysum,yasum,yzsum,abar,zbar

  ! calculate neutron separation energy
  if (calc_nsep_energy) then
     call calc_nseparation_energy(Y)
  endif

  ! output diagnostics
  Ye = sum(isotope(1:net_size)%p_nr*Y(1:net_size))
  ylight=   max(1.0d-99,sum(Y(ipro+1:ihe4-1)))
  yheavies= 0.0
  if (ihe4.ge.1) yheavies= max(1.0d-99,sum(Y(ihe4+1:net_size)))

  ysum  = sum(Y(1:net_size))
  yasum = sum(Y(1:net_size)*isotope(1:net_size)%mass)
  yzsum = sum(Y(1:net_size)*isotope(1:net_size)%p_nr)
  abar  = yasum/ysum
  zbar  = yzsum/ysum
  !zbar = yzsum/sum(Y(2:net_size)) !as neutrons have Z=0

  Yneu= 0.0
  if (ineu.ge.1) Yneu= max(1.0d-99,Y(ineu))

  Ypro= 0.0
  if (ipro.ge.1) Ypro= max(1.0d-99,Y(ipro))

  Yhe4= 0.0
  if (ihe4.ge.1) Yhe4= max(1.0d-99,Y(ihe4))

  if (calc_nsep_energy) then
    write (mainout_unit,'(i8,X,100(es23.16,2X))') cnt,ctime,T9,rho_b,Ye,Rkm, &
       Yneu,Ypro,Yhe4,ylight,yheavies,zbar,abar,entropy,pressure,Sn_ave
  else
    write (mainout_unit,'(i8,X,100(es23.16,2X))') cnt,ctime,T9,rho_b,Ye,Rkm, &
       Yneu,Ypro,Yhe4,ylight,yheavies,zbar,abar,entropy,pressure
  endif

end subroutine output_mainout


!>
!! Output Snapshot files
!!
!! Snapshots store the abundances and mass fractions of all nuclei
!! at certain times. An example file looks like:
!! \file{
!! time    temp    dens
!! 9.77198412010000E-02   9.99927332947922E+00   3.88380373937259E+06
!!  nin     zin       y       x
!! 1    0     5.72194135572689E-01   5.72194135572689E-01
!! 0    1     3.92194484636783E-01   3.92194484636783E-01
!! 1    1     7.83033243870819E-05   1.56606648774164E-04
!! ...}
!!
!! \b Edited:
!!           - 03.04.18, M. Jacobi
!! .
subroutine output_snapshot(ctime,T9,rho_b,Y)
  implicit none

  real(r_kind),intent(in) :: ctime                 !< current time
  real(r_kind),intent(in) :: T9                    !< temperature [GK]
  real(r_kind),intent(in) :: rho_b                 !< density [gcc]
  real(r_kind),dimension(net_size),intent(in):: Y  !< array of abundances Y_i


  ! Output Snapshots
  call printsnap (ctime,T9,rho_b,Y)

end subroutine output_snapshot


!>
!! Output flow files
!!
!! \b Edited:
!!           - 19.02.20, M. Jacobi
subroutine output_flow(ctime,T9,rho_b,Y)
  use flow_module,     only: flowprint, flowcalc
  implicit none

  real(r_kind),intent(in) :: ctime                 !< current time
  real(r_kind),intent(in) :: T9                    !< temperature [GK]
  real(r_kind),intent(in) :: rho_b                 !< density [gcc]
  real(r_kind),dimension(net_size),intent(in):: Y  !< array of abundances Y_i

  call flowcalc(Y)
  call flowprint(ctime,T9,rho_b,Y,flowcnt)
  flowcnt = flowcnt + 1
end subroutine output_flow


!>
!! Output of the Newton-Raphson loop diagnostics
!!
!! The diagnostic output contains NR iterations, timesteps,
!! the nucleus that restricts the timestep and more.
!! An example file looks like:
!! \l_file{
!! # Newton-Raphson diagnostic output
!! # 1 : global iteration count (cnt)
!! # 2 : global time [s]
!! # 3 : temperature T9 [GK]
!! # 4 : adapt stepsize loop counter (k)
!! # 5 : NR loop counter (nr_count)
!! # 6 : global timestep [s]
!! # 7 : mass conservation error
!! # 8 : maximal abundance change (eps)
!! # 9 : most rapidly evolving isotope (epsl)
!! # 10: abundance of the isotope epsl, y(epsl)
!!       0   5.20740000000E-22   7.42450000000E+00    0    1   5.20740000000E-22   9.32587340685E-15   2.22044604925E-16 cl42   2.43899441583E-08
!!       0   5.20740000000E-22   7.42450000000E+00    0    2   5.20740000000E-22   9.32587340685E-15   2.22044604925E-16  s40   2.73320077046E-08
!!       1   1.56220000000E-21   7.42450000000E+00    0    1   1.04146000000E-21   9.32587340685E-15   2.22044604925E-16  he4   3.13811866442E-03
!!       1   1.56220000000E-21   7.42450000000E+00    0    2   1.04146000000E-21   9.32587340685E-15   2.22044604925E-16 mg30   1.57736344914E-10
!! ...}
subroutine output_nr_diagnostic(cnt,k,nr_c,ctime,T9,stepsize,mtot,Y_p,Y)
  use parameter_class, only: nrdiag_every
  use global_class, only: isotope
  implicit none
  integer,intent(in)      :: cnt               !< global iteration counter
  integer,intent(in)      :: k                 !< counter in adapt_stepsize loop
  integer,intent(in)      :: nr_c              !< counter in the NR loop
  real(r_kind),intent(in) :: ctime             !< global current time
  real(r_kind),intent(in) :: T9                !< current temperature [GK]
  real(r_kind),intent(in) :: stepsize          !< current step size
  real(r_kind),intent(in) :: mtot              !< total mass
  real(r_kind),dimension(net_size) :: Y_p,Y    !< old/new abundances
  !
  integer      :: i
  real(r_kind) :: eps,epst
  integer      :: epsl

  INFO_ENTRY("output_nr_diagnostic")

     ! check periodicity condition
     if (nrdiag_every.gt.0) then
        if (mod(cnt,nrdiag_every) .ne. 0) return
     else
        return
     endif

     ! which isotope has changed the most
     eps = 0.d0
     epsl= 1
     do i=1,net_size
        if(Y_p(i).lt.1.d-10) cycle
        if(dabs(Y_p(i)-Y(i)).eq.0.d0) cycle
        epst = dabs(Y(i)/Y_p(i)-1.d0)
        if (epst.gt.eps) then
           eps = epst
           epsl = i
        end if
     end do

     ! output convergence diagnostic
     write(nrdiag_unit, &
        '(i7,2(1x,es19.11),2(1x,i4),3(1x,es19.11),a5,x,es19.11)')           &
        cnt,ctime,T9,k,nr_c,stepsize,dabs(mtot-1d0),eps,isotope(epsl)%name, &
        Y(epsl)


  INFO_EXIT("output_nr_diagnostic")

end subroutine output_nr_diagnostic


!>
!! Print final statistics
!!
!! An example looks like:
!! \file{
!! Final time reached.
!! Final time:  2.000000000000E+17
!! Total number of iterations:       5334
!! Elapsed time [s]:  5.41572E+02
!! }
!!
!! \b Edited:
!!  - 2023-02-06, M.R. took care of wrapping of the simulation time
!!                     that could lead to negative simulation times.
!! .
subroutine output_final_stats(cnt,ctime)
  implicit none
  integer,intent(in)      :: cnt
  real(r_kind),intent(in) :: ctime
  !
  real(r_kind) :: simtime

  INFO_ENTRY('output_final_stats')

  ! calculate elapsed time, check wrapping at cl_cmax
  if (cl_end .gt. cl_start) then
    simtime= float(cl_end-cl_start)/float(cl_rate)
  else
    simtime= (float(cl_end)+float(cl_cmax)-float(cl_start))/float(cl_rate)
  end if

  ! output final statistics
  print '(A,ES19.12)', 'Final time: ', ctime
  print '(A,I10)',     'Total number of iterations: ', cnt
  print '(A,ES12.5)',  'Elapsed time [s]: ', simtime

  INFO_EXIT('output_final_stats')

end subroutine output_final_stats


!>
!! Print current snapshot
!!
!! Snapshots store the abundances and mass fractions of all nuclei
!! at certain times. An example file looks like:
!! \file{
!! time    temp    dens
!! 9.77198412010000E-02   9.99927332947922E+00   3.88380373937259E+06
!!  nin     zin       y       x
!! 1    0     5.72194135572689E-01   5.72194135572689E-01
!! 0    1     3.92194484636783E-01   3.92194484636783E-01
!! 1    1     7.83033243870819E-05   1.56606648774164E-04
!! ...}
subroutine printsnap (t,T9,rho_b,Y)
  use global_class, only: isotope
  implicit none

  real(r_kind),intent(in)                :: t,T9,rho_b
  real(r_kind),dimension(net_size),intent(in) :: Y
  character (len=21) :: snapfile
  integer :: snapunit, z
  real(r_kind) :: Xz

  INFO_ENTRY("printsnap")

  write(snapfile,'("snaps/snapsh_",i4.4,".dat")') snapcnt
  snapunit= open_outfile(snapfile)

  write(snapunit,'(3a8)')'time','temp','dens'
  write(snapunit,'(3es23.14)')t,T9,rho_b
  write(snapunit,'(4(a8))')'nin','zin','y','x'
  do z=1,net_size
     Xz= Y(z)*isotope(z)%mass
     write(snapunit,'(2(i3,2x),2(es23.14))') &
          isotope(z)%n_nr, isotope(z)%p_nr, max(1d-99,Y(z)), max(1d-99,Xz)
  end do
  call close_io_file(snapunit,snapfile)
  snapcnt= snapcnt + 1
  tsnp= t
  T9snp= T9
  rosnp=rho_b

  INFO_EXIT("printsnap")

end subroutine printsnap


!>
!! Calculates average neutron separation energy
!!
!! This neutron separation energy is printed in the mainout when
!! parameter_class::calc_nsep_energy was set to "yes" and
!! a valid file for parameter_class::nsep_energies_file.
subroutine calc_nseparation_energy(Y)
  use global_class, only: isotope
  implicit none

  real(r_kind),dimension(net_size),intent(in) :: Y  !< abundances
  !
  real(r_kind) :: ysum,yasum,yzsum,abar_ions,zbar_ions
  integer :: i

  INFO_ENTRY("calc_nseparation_energy")

  abar_ions = 0d0
  zbar_ions = 0d0
  do i=max(1,ihe4),net_size
     abar_ions = abar_ions + isotope(i)%mass * Y(i)
     zbar_ions = zbar_ions + isotope(i)%p_nr * Y(i)
     !Sn_ave = Sn_ave + Sn(i) * Y(i)
  enddo
  ysum = sum(Y(max(1,ihe4):net_size))
  abar_ions = abar_ions / ysum
  zbar_ions = zbar_ions / ysum
  ysum      = 0d0
  Sn_ave    = 0d0
  do i=max(1,ihe4),net_size
     if (Sn(i) /= 0) then
        if ((isotope(i)%mass.ge.220).and.(isotope(i)%mass.le.260)) then
           Sn_ave = Sn_ave + Sn(i) * Y(i)
           ysum = ysum + Y(i)
        end if
     endif
  enddo
  Sn_ave = Sn_ave / ysum
  ysum = sum(Y)
  yasum = sum(Y*isotope%mass)
  yzsum = sum(Y*isotope%p_nr)

  INFO_EXIT("calc_nseparation_energy")

end subroutine calc_nseparation_energy



!>
!! @brief: Calculate the generated energy for each class of reactions separately.
!!
!! This subroutine calculates the generated energy of the nuclear reactions
!! with the help of the Q-values. This energy is written to an output only.
!! It is not fed back to the temperature and has therefore
!! no impact on the final result.
!!
!! @see \ref nuclear_heating
!!
!! @author: Moritz Reichert
!! @date   29.02.20
subroutine calc_nuclear_heating(hdf5_mode)
   use global_class,        only: reactionrate_type,rrate,nreac,ihe4,ineu,ipro
   use fission_rate_module, only: fissionrate_type, fissrate, nfiss
   use nucstuff_class,      only: isotope,inter_partf
   use single_zone_vars,    only: Y,Y_p,T9,rhob,time,stepsize,evolution_mode,Ye
   use parameter_class,     only: unit,fissflag,h_engen_detailed
   use hdf5_module
   implicit none
   logical, intent(in)     :: hdf5_mode  !< whether to store as hdf5 or ascii
   type(reactionrate_type) :: rr_tmp     !< auxiliary shorthand for reaction
   type(fissionrate_type)  :: fr_tmp     !< auxiliary shorthand for fission reaction
   character*4  :: src                   !< reaction source
   integer      :: i                     !< Loop variable
   integer      :: help_parts            !< Helper variable
   real(r_kind) :: rat                   !< storage of the cached reaction rate
   ! Heating variables
   real(r_kind) :: engen_q_value         !< All reactions
   real(r_kind) :: engen_q_weak          !< weak reactions
   real(r_kind) :: engen_q_ng,engen_q_gn !< gamma-n reactions
   real(r_kind) :: engen_q_pg,engen_q_gp !< gamma-p reactions
   real(r_kind) :: engen_q_ga,engen_q_ag !< gamma-alpha reactions
   real(r_kind) :: engen_q_pa,engen_q_ap !< p-alpha reactions
   real(r_kind) :: engen_q_pn,engen_q_np !< p-n reactions
   real(r_kind) :: engen_q_na,engen_q_an !< n-alpha reactions
   real(r_kind) :: engen_fiss            !< fission reactions
   real(r_kind) :: engen_tot_ncap, engen_tot_pcap, engen_tot_acap, &
                   engen_tot_npcap,engen_tot_nacap, engen_tot_pacap !< net energy
   real(r_kind) :: engen_tot_exc         !< total energy with mass excess formula
   real(r_kind) :: entropy_src_term      !< Source term of entropy * kBT
   real(r_kind) :: mic2,gi,twopi_h2c2,mui!< Helper variables
   real(r_kind) :: kBT                   !< Helper variables
   real(r_kind),dimension(0:net_size)::pf!< partition functions
   real(r_kind) :: etaele,etapos         !< electron and positron chemical potental
   real(r_kind) :: mue                   !< electron chemical potental [MeV]
   ! Helpers
   real(r_kind) :: tmp_energy                         !< Temporary energy per reaction
   real(r_kind),dimension(net_size)  ::det_bet_engen  !< Energy of decays per parent nucleus
   real(r_kind),dimension(net_size)  ::det_ag_engen   !< Energy of (a,g)+(g,a) per parent nucleus
   real(r_kind),dimension(net_size)  ::det_pg_engen   !< Energy of (p,g)+(g,p) per parent nucleus
   real(r_kind),dimension(net_size)  ::det_ng_engen   !< Energy of (n,g)+(g,n) per parent nucleus
   real(r_kind),dimension(net_size)  ::det_pn_engen   !< Energy of (p,n)+(n,p) per parent nucleus
   real(r_kind),dimension(net_size)  ::det_ap_engen   !< Energy of (a,p)+(p,a) per parent nucleus
   real(r_kind),dimension(net_size)  ::det_an_engen   !< Energy of (a,n)+(n,a) per parent nucleus
   real(r_kind),dimension(net_size)  ::det_other_engen!< Energy of other reactions per parent nucleus
   real(r_kind),dimension(net_size)  ::det_fiss_engen !< Energy of fission per parent nucleus
   real(r_kind),dimension(net_size)  ::det_tot_engen  !< Total energy per parent nucleus


   INFO_ENTRY("calc_nuclear_heating")

   ! Calculate the heating with help of the q-value
   engen_q_value = 0 ! Total energy [erg/g/s]
   engen_q_weak  = 0 ! Energy generated by weak reactions [erg/g/s]
   engen_q_ng    = 0 ! Energy generated by n-gamma reactions [erg/g/s]
   engen_q_gn    = 0 ! Energy generated by gamma-n reactions [erg/g/s]
   engen_q_pg    = 0 ! Energy generated by p-gamma reactions [erg/g/s]
   engen_q_gp    = 0 ! Energy generated by gamma-p reactions [erg/g/s]
   engen_q_ga    = 0 ! Energy generated by gamma-alpha reactions [erg/g/s]
   engen_q_ag    = 0 ! Energy generated by alpha-gamma reactions [erg/g/s]
   engen_q_np    = 0 ! Energy generated by n-p reactions [erg/g/s]
   engen_q_pn    = 0 ! Energy generated by p-n reactions [erg/g/s]
   engen_q_na    = 0 ! Energy generated by n-alpha reactions [erg/g/s]
   engen_q_an    = 0 ! Energy generated by alpha-n reactions [erg/g/s]
   engen_q_pa    = 0 ! Energy generated by p-alpha reactions [erg/g/s]
   engen_q_ap    = 0 ! Energy generated by alpha-p reactions [erg/g/s]
   engen_fiss    = 0 ! Energy generated by fission [erg/g/s]

   ! Initialize decay energy
   if (h_engen_detailed) then
       det_bet_engen(:)   = 0.0
       det_ag_engen(:)    = 0.0
       det_pg_engen(:)    = 0.0
       det_ng_engen(:)    = 0.0
       det_pn_engen(:)    = 0.0
       det_ap_engen(:)    = 0.0
       det_an_engen(:)    = 0.0
       det_other_engen(:) = 0.0
       det_fiss_engen(:)  = 0.0
       det_tot_engen(:)   = 0.0
   end if

   ! For NSE the reaction rates are not cached and the generated energy would be therefore weird
   ! Set it to zero for this case
   if (evolution_mode .ne. EM_NSE) then
       ! Loop through all reactions
       do i=1,nreac
          rr_tmp = rrate(i)
          src = rr_tmp%source
          rat = rr_tmp%cached

          ! Do nothing if the rate is zero
          if (rat .eq. 0) cycle

          ! dYdt depends on the reaclib chapter and we therefore have to separate them
          select case (rr_tmp%group)
          case(1:3,11)! Reaclib chapter 1-3 and 11
             ! Energy per reaction
             tmp_energy = rat*Y(rr_tmp%parts(1)) * rr_tmp%q_value
             ! Total energy
             engen_q_value = engen_q_value + tmp_energy
             ! Seperate the reaction types, weak reactions and photodissociations should be in chapter 1-3

             ! weak rates (only contained in chapter 1-3)
             if (rr_tmp%is_weak) then
                engen_q_weak = engen_q_weak + tmp_energy
                if (h_engen_detailed) then
                  det_bet_engen(rr_tmp%parts(1)) = det_bet_engen(rr_tmp%parts(1))+tmp_energy
                end if
             end if
             ! (gamma,n) reactions
             if(rr_tmp%reac_type.eq.rrt_gn) then
                engen_q_gn = engen_q_gn + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(2) .eq. ineu) then
                    help_parts = rr_tmp%parts(3)
                  else
                    help_parts = rr_tmp%parts(2)
                  end if
                  det_ng_engen(help_parts) = det_ng_engen(help_parts)+tmp_energy
                end if

             endif
             ! (gamma,p) reactions
             if(rr_tmp%reac_type.eq.rrt_gp) then
                engen_q_gp = engen_q_gp + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(2) .eq. ipro) then
                    help_parts = rr_tmp%parts(3)
                  else
                    help_parts = rr_tmp%parts(2)
                  end if
                  det_pg_engen(help_parts) = det_pg_engen(help_parts)+tmp_energy
                end if
             endif
             ! (gamma,alpha) reactions
             if(rr_tmp%reac_type.eq.rrt_ga) then
                engen_q_ga = engen_q_ga + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(2) .eq. ihe4) then
                    help_parts = rr_tmp%parts(3)
                  else
                    help_parts = rr_tmp%parts(2)
                  end if
                  det_ag_engen(help_parts) = det_ag_engen(help_parts)+tmp_energy
                end if
             endif


            if(rr_tmp%reac_type.eq.rrt_o) then
              if (h_engen_detailed) then
                det_other_engen(rr_tmp%parts(1)) = det_other_engen(rr_tmp%parts(1))+tmp_energy
              end if
            end if

          case(4:7)! Reaclib chapter 4-7
             ! Energy per reaction
             tmp_energy = rat*Y(rr_tmp%parts(1)) * Y(rr_tmp%parts(2)) &
                          * rr_tmp%q_value
             ! Total energy
             engen_q_value = engen_q_value + tmp_energy

             ! Seperate the reaction types
             ! (n-gamma) reactions
             if(rr_tmp%reac_type.eq.rrt_ng) then
                engen_q_ng = engen_q_ng + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(2) .eq. ineu) then
                    help_parts = rr_tmp%parts(1)
                  else
                    help_parts = rr_tmp%parts(2)
                  end if
                  det_ng_engen(help_parts) = det_ng_engen(help_parts)+tmp_energy
                end if
             endif
             ! (p-gamma) reactions
             if(rr_tmp%reac_type.eq.rrt_pg) then
                engen_q_pg = engen_q_pg + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(2) .eq. ipro) then
                    help_parts = rr_tmp%parts(1)
                  else
                    help_parts = rr_tmp%parts(2)
                  end if
                  det_pg_engen(help_parts) = det_pg_engen(help_parts)+tmp_energy
                end if
             endif
             ! (alpha gamma) reactions
             if(rr_tmp%reac_type.eq.rrt_ag) then
                engen_q_ag = engen_q_ag + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(2) .eq. ihe4) then
                    help_parts = rr_tmp%parts(1)
                  else
                    help_parts = rr_tmp%parts(2)
                  end if
                  det_ag_engen(help_parts) = det_ag_engen(help_parts)+tmp_energy
                end if
             endif
             ! (p,alpha) reactions
             if(rr_tmp%reac_type.eq.rrt_pa) then
                engen_q_pa = engen_q_pa + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(3) .eq. ihe4) then
                    help_parts = rr_tmp%parts(4)
                  else
                    help_parts = rr_tmp%parts(3)
                  end if
                  det_ap_engen(help_parts) = det_ap_engen(help_parts)+tmp_energy
                end if
             endif
             ! (alpha,p) reactions
             if(rr_tmp%reac_type.eq.rrt_ap) then
                engen_q_ap = engen_q_ap + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(2) .eq. ihe4) then
                    help_parts = rr_tmp%parts(1)
                  else
                    help_parts = rr_tmp%parts(2)
                  end if
                  det_ap_engen(help_parts) = det_ap_engen(help_parts)+tmp_energy
                end if
             endif
             ! (n,p) reactions
             if(rr_tmp%reac_type.eq.rrt_np) then
                engen_q_np = engen_q_np + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(3) .eq. ipro) then
                    help_parts = rr_tmp%parts(4)
                  else
                    help_parts = rr_tmp%parts(3)
                  end if
                  det_pn_engen(help_parts) = det_pn_engen(help_parts)+tmp_energy
                end if
             endif
             ! (p,n) reactions
             if(rr_tmp%reac_type.eq.rrt_pn) then
                engen_q_pn = engen_q_pn + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(1) .eq. ipro) then
                    help_parts = rr_tmp%parts(2)
                  else
                    help_parts = rr_tmp%parts(1)
                  end if
                  det_pn_engen(help_parts) = det_pn_engen(help_parts)+tmp_energy
                end if
             endif
             ! (alpha,n) reactions
             if(rr_tmp%reac_type.eq.rrt_an) then
                engen_q_an = engen_q_an + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(1) .eq. ihe4) then
                    help_parts = rr_tmp%parts(2)
                  else
                    help_parts = rr_tmp%parts(1)
                  end if
                  det_an_engen(help_parts) = det_an_engen(help_parts)+tmp_energy
                end if
             endif
             ! (n,alpha) reactions
             if(rr_tmp%reac_type.eq.rrt_na) then
                engen_q_na = engen_q_na + tmp_energy
                ! Detailed output
                if (h_engen_detailed) then
                  if (rr_tmp%parts(3) .eq. ihe4) then
                    help_parts = rr_tmp%parts(4)
                  else
                    help_parts = rr_tmp%parts(3)
                  end if
                  det_an_engen(help_parts) = det_an_engen(help_parts)+tmp_energy
                end if
             endif

             if(rr_tmp%reac_type.eq.rrt_o) then
               if (h_engen_detailed) then
                 det_other_engen(rr_tmp%parts(1)) = det_other_engen(rr_tmp%parts(1))+tmp_energy
               end if
             end if

          case(8:9)! Reaclib chapter 8-9
             ! Total energy
             engen_q_value = engen_q_value + rat*Y(rr_tmp%parts(1)) * Y(rr_tmp%parts(2)) &
                  * Y(rr_tmp%parts(3)) * rr_tmp%q_value

             if(rr_tmp%reac_type.eq.rrt_o) then
               if (h_engen_detailed) then
                 det_other_engen(rr_tmp%parts(1)) = det_other_engen(rr_tmp%parts(1))+tmp_energy
               end if
             end if

          case(10)! Reaclib chapter 8
            ! Total energy
            engen_q_value = engen_q_value + rat*Y(rr_tmp%parts(1)) * Y(rr_tmp%parts(2)) &
                          * Y(rr_tmp%parts(3)) * Y(rr_tmp%parts(4))* rr_tmp%q_value

            if(rr_tmp%reac_type.eq.rrt_o) then
              if (h_engen_detailed) then
                det_other_engen(rr_tmp%parts(1)) = det_other_engen(rr_tmp%parts(1))+tmp_energy
              end if
            end if


          case default ! This should not happen
             cycle
          end select
       end do

       ! Do the same for the fission reactions
       if (fissflag.ne.0) then
          fission:do i=1,nfiss
             fr_tmp = fissrate(i)
             rat = fr_tmp%cached
             if (rat .eq. 0) cycle
             if ((fr_tmp%mode.eq.2).or.(fr_tmp%mode.eq.3)) then  ! spont. and b-delayed fission
                engen_q_value = engen_q_value + rat*Y(fissrate(i)%fissparts(1)) * fissrate(i)%averageQ
                engen_fiss    = engen_fiss + rat*Y(fissrate(i)%fissparts(1)) * fissrate(i)%averageQ
                tmp_energy    = rat*Y(fissrate(i)%fissparts(1)) * fissrate(i)%averageQ
             else if (fr_tmp%mode.eq.1) then ! n-induced fission
                engen_q_value = engen_q_value + rat*Y(fissrate(i)%fissparts(1))*Y(fissrate(i)%fissparts(2)) * fissrate(i)%averageQ
                engen_fiss    = engen_fiss + rat*Y(fissrate(i)%fissparts(1))*Y(fissrate(i)%fissparts(2)) * fissrate(i)%averageQ
                tmp_energy    = rat*Y(fissrate(i)%fissparts(1))*Y(fissrate(i)%fissparts(2)) * fissrate(i)%averageQ
             end if

             if (h_engen_detailed) then
               det_fiss_engen(fissrate(i)%fissparts(1)) = det_fiss_engen(fissrate(i)%fissparts(1))+tmp_energy
             end if

          end do fission
       end if
   end if

   engen_tot_exc = 0
   do i=1,net_size
     if (Y(i).gt.1d-25) then
         engen_tot_exc = engen_tot_exc -(isotope(i)%mass_exc)*(Y(i)-Y_p(i))  ! energy generated by species i [MeV/baryon]
         if (h_engen_detailed) then
           det_tot_engen(i) = -(isotope(i)%mass_exc)*(Y(i)-Y_p(i))
         end if
     endif
   end do
   engen_tot_exc = engen_tot_exc/stepsize / unit%hix

   entropy_src_term = 0.d0
   call inter_partf (T9, pf)

   ! Take care that one does not go off-grid
   if (rhob .gt. 1d-10) then
    ! Get the electron chemical potential
    call chempot(T9*1d9,rhob,Ye,etaele,etapos)
    mue = etaele * unit%k_mev * T9*1d9
   else
    mue = 0.d0
   end if
   mue = mue + unit%mass_e

   kBT = T9*1d9 * unit%k_mev ! 8.617332d-11 * T [MeV]
   do i=1,net_size
     if (Y(i).gt.1d-25) then
         mic2 = isotope(i)%mass*unit%mass_u - isotope(i)%p_nr * unit%mass_e + isotope(i)%mass_exc    ! [MeV]
         gi = 2d0*isotope(i)%spin + 1d0                               ! statistical weight
         twopi_h2c2 = 2d0*unit%pi/(unit%h_mevs*unit%clight)**2        ! 2\pi/(hc)^2 [MeV^-2 cm^-2]
         mui = -kBT*log( gi*pf(i) / (rhob*Y(i)*unit%n_a) * &
             (twopi_h2c2 * mic2 * kBT)**1.5d0 )
         ! entropy generated by species i [MeV/baryon]
         entropy_src_term = entropy_src_term -(mic2 + mui + isotope(i)%p_nr*mue)*(Y(i)-Y_p(i))
     endif
   end do

   entropy_src_term = entropy_src_term/ unit%hix/stepsize
   ! Convert from MeV/Baryon/s - > Erg/g/s
   engen_q_value = engen_q_value / unit%hix
   engen_q_weak  = engen_q_weak  / unit%hix
   engen_q_ng    = engen_q_ng    / unit%hix
   engen_q_gn    = engen_q_gn    / unit%hix
   engen_q_pg    = engen_q_pg    / unit%hix
   engen_q_gp    = engen_q_gp    / unit%hix
   engen_q_ga    = engen_q_ga    / unit%hix
   engen_q_ag    = engen_q_ag    / unit%hix
   engen_q_np    = engen_q_np    / unit%hix
   engen_q_pn    = engen_q_pn    / unit%hix
   engen_q_na    = engen_q_na    / unit%hix
   engen_q_an    = engen_q_an    / unit%hix
   engen_q_pa    = engen_q_pa    / unit%hix
   engen_q_ap    = engen_q_ap    / unit%hix
   engen_fiss    = engen_fiss    / unit%hix
   ! Detailed output
   if (h_engen_detailed) then
     det_tot_engen   = det_tot_engen/stepsize / unit%hix
     det_bet_engen   = det_bet_engen   / unit%hix
     det_ag_engen    = det_ag_engen    / unit%hix
     det_pg_engen    = det_pg_engen    / unit%hix
     det_ng_engen    = det_ng_engen    / unit%hix
     det_pn_engen    = det_pn_engen    / unit%hix
     det_ap_engen    = det_ap_engen    / unit%hix
     det_an_engen    = det_an_engen    / unit%hix
     det_other_engen = det_other_engen / unit%hix
     det_fiss_engen  = det_fiss_engen  / unit%hix
   end if

   !Calculate net energy
   engen_tot_ncap = engen_q_ng + engen_q_gn
   engen_tot_pcap = engen_q_pg + engen_q_gp
   engen_tot_acap = engen_q_ag + engen_q_ga
   engen_tot_npcap= engen_q_np + engen_q_pn
   engen_tot_nacap= engen_q_na + engen_q_an
   engen_tot_pacap= engen_q_pa + engen_q_ap

   if (.not. hdf5_mode) then
      ! Write the output (detailed)
      !write(nuc_heat_file,'(*(es23.14,3x))') time,engen_q_value,engen_q_weak,&
      !                                       engen_q_ng,engen_q_gn,engen_tot_ncap,&
      !                                       engen_q_pg,engen_q_gp,engen_tot_pcap,&
      !                                       engen_q_ag,engen_q_ga,engen_tot_acap,&
      !                                       engen_q_np,engen_q_pn,engen_tot_npcap,&
      !                                       engen_q_na,engen_q_an,engen_tot_nacap,&
      !                                       engen_q_pa,engen_q_ap,engen_tot_pacap,&
      !                                       engen_fiss
      ! Write the output (more relaxed)
      write(nuc_heat_file,'(*(es23.14,3x))') time,engen_tot_exc,entropy_src_term,&
                                             engen_q_weak,&
                                             engen_tot_ncap,&
                                             engen_tot_pcap,&
                                             engen_tot_acap,&
                                             engen_tot_npcap,&
                                             engen_tot_nacap,&
                                             engen_tot_pacap,&
                                             engen_fiss

#ifdef USE_HDF5
   else
      call extend_engen(time,engen_tot_exc, entropy_src_term,engen_tot_ncap, &
                        engen_tot_pcap, engen_tot_acap,engen_tot_npcap, &
                        engen_tot_nacap, engen_tot_pacap, engen_q_weak, &
                        engen_fiss, &
                        det_bet_engen, det_ag_engen, det_pg_engen, &
                        det_ng_engen, det_pn_engen, det_ap_engen, &
                        det_an_engen, det_other_engen, det_fiss_engen, &
                        det_tot_engen, &
                        h_engen_detailed)
#endif

   end if
   INFO_EXIT("calc_nuclear_heating")

end subroutine calc_nuclear_heating






!>
!! Calculate average timescales of different classes of reactions
!!
!! These timescales are calculated as in
!! [Arcones et al. 2012](https://ui.adsabs.harvard.edu/abs/2012ApJ...750...18A/abstract)
!! (Eq. 1-6). The (p,gamma) timescale is given by e.g.,
!! \f[
!! \tau_{p\gamma} = \left[ \frac{\rho Y_p}{Y_h}N_A \sum \langle \sigma \nu \rangle_{p,\gamma}(Z,A)Y(Z,A) \right]^{-1}
!! \f]
!!
!! An example file looks like:
!! \l_file{
!! # time[s]     Temperature [GK]     tau_ng [s]     tau_gn [s]     tau_pg [s]     tau_gp [s]     tau_np [s]     tau_pn [s]     tau_an [s]     tau_na [s]     tau_beta [s]     tau_alpha [s]     tau_nfiss [s]     tau_sfiss [s]     tau_bfiss [s]
!!  1.1130000000000000E+00    9.1821750000000009E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    1.2642767370327558E+01    2.9487624276191804E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00
!!  1.1220000000000001E+00    8.2638489999999951E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00    1.0557999180036711E+02    2.5923655028366777E+01    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00
!!  1.1247618491767180E+00    7.9904974735555152E+00    9.5690472538317728E-13    5.3150487750907188E-11    7.4419112944811684E-10    9.8479480155439862E-08    1.9792689873245490E-12    1.9444248037144314E-12    4.8686627517806918E-10    3.9657685448682009E-12    1.4600540103519526E+01    5.2634156341152767E+01    0.0000000000000000E+00    0.0000000000000000E+00    0.0000000000000000E+00
!! ...}
subroutine calc_av_timescales(ctime,T9,Y,hdf5)
  use global_class,        only: reactionrate_type,rrate,nreac,ihe4,ineu,ipro
  use fission_rate_module, only: fissionrate_type, fissrate, nfiss
  use file_handling_class
#ifdef USE_HDF5
  use hdf5_module,         only: extend_av_timescales
#endif
  implicit none
  real(r_kind),intent(in)       :: ctime !< current time [s]
  real(r_kind),intent(in)       :: T9    !< Temperature [GK]
  logical,optional,intent(in)   :: hdf5  !< Whether the output is written to hdf5 or not
  real(r_kind),dimension(net_size),intent(in)   :: Y  !< abundances

  ! average timescales
  real(r_kind) :: tau_ng           !< (n,gamma)
  real(r_kind) :: tau_gn           !< (gamma,n)
  real(r_kind) :: tau_pg           !< (p,gamma)
  real(r_kind) :: tau_gp           !< (gamma,p)
  real(r_kind) :: tau_np           !< (n,p)
  real(r_kind) :: tau_pn           !< (p,n)
  real(r_kind) :: tau_an           !< (alpha,n)
  real(r_kind) :: tau_na           !< (n,alpha)
  real(r_kind) :: tau_ag           !< (alpha,gamma)
  real(r_kind) :: tau_ga           !< (gamma,alpha)
  real(r_kind) :: tau_pa           !< (p,alpha)
  real(r_kind) :: tau_ap           !< (alpha,p)
  !real(r_kind) :: tau_npro        !< neutron production timescale
  !real(r_kind) :: tau_fdis        !< photodissociation timescale
  real(r_kind) :: tau_beta         !< beta-decay timescale
  real(r_kind) :: tau_alpha        !< alpha-decay timescale
  real(r_kind) :: tau_nfiss        !< n-induced fission timescale
  real(r_kind) :: tau_sfiss        !< spontaneous fission timescale
  real(r_kind) :: tau_bfiss        !< beta-delayed fission timescale
  integer      :: nuc              !< auxiliary var
  real(r_kind) :: rat              !< Value of reaction rate
  real(r_kind) :: ng_sum,gn_sum,pg_sum,gp_sum,beta_sum
  real(r_kind) :: np_sum,pn_sum,an_sum,na_sum,ag_sum,ga_sum
  real(r_kind) :: ap_sum,pa_sum
  real(r_kind) :: alpha_sum, nfiss_sum, sfiss_sum, bfiss_sum
  real(r_kind) :: ysum               !< Sum of abundances heavier than helium
  integer      :: i                  !< Loop variable
  type(reactionrate_type) :: rr_tmp  !< auxiliary shorthand for reaction
  type(fissionrate_type)  :: fr_tmp  !< auxiliary shorthand for reaction
  character*4  :: src                !< reaction source
  real(r_kind) :: q_sm = 1.e-20      !< small quantity to add to denominator

  INFO_ENTRY("calc_av_timescales")


  ! to avoid division by zero
  ng_sum    = q_sm
  gn_sum    = q_sm
  pg_sum    = q_sm
  gp_sum    = q_sm
  np_sum    = q_sm
  pn_sum    = q_sm
  an_sum    = q_sm
  na_sum    = q_sm
  ag_sum    = q_sm
  ga_sum    = q_sm
  ap_sum    = q_sm
  pa_sum    = q_sm
  beta_sum  = q_sm
  alpha_sum = q_sm
  nfiss_sum = q_sm
  sfiss_sum = q_sm
  bfiss_sum = q_sm


  ysum = sum(Y(max(1,ihe4):net_size))

  do i=1,nreac
     rr_tmp = rrate(i)
     src = rr_tmp%source
     rat = rr_tmp%cached

     ! beta-decays (not including beta-delayed fission)
     if ((rr_tmp%group.eq.1).or.(rr_tmp%group.eq.2).or.(rr_tmp%group.eq.3) .or.(rr_tmp%group.eq.11)) then
        if ((rr_tmp%reac_type .eq.rrt_betm).or.(rr_tmp%reac_type .eq.rrt_betp).and.(rr_tmp%source.ne.'ms99')) then
        ! only include significant reactions
           nuc = rr_tmp%parts(1)
           if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
              beta_sum = beta_sum + rat*Y(nuc)
           end if
        endif
     endif

     ! alpha-decays
     if ((rr_tmp%reac_type.eq.rrt_alpd)) then
        nuc = rr_tmp%parts(1)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           alpha_sum = alpha_sum + rat*Y(nuc)
        end if
     end if

     ! (a,p)
     if(rr_tmp%reac_type.eq.rrt_ap) then
        nuc = rr_tmp%parts(2)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           ap_sum  = ap_sum  + rat*Y(nuc)*Y(ihe4)
        end if
     endif

     ! (p,a)
     if(rr_tmp%reac_type.eq.rrt_pa) then
        nuc = rr_tmp%parts(2)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           pa_sum  = pa_sum  + rat*Y(nuc)*Y(ipro)
        end if
     endif

     ! (n,gamma)
     if(rr_tmp%reac_type.eq.rrt_ng) then
        nuc = rr_tmp%parts(2)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           ng_sum  = ng_sum  + rat*Y(nuc)*Y(ineu)
        end if
     endif

     ! (gamma,n)
     if(rr_tmp%reac_type.eq.rrt_gn) then
        nuc = rr_tmp%parts(1)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           gn_sum  = gn_sum  + rat*Y(nuc)
        end if
     endif

     ! (p,gamma)
     if(rr_tmp%reac_type.eq.rrt_pg) then
        nuc = rr_tmp%parts(2)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           pg_sum  = pg_sum  + rat*Y(nuc)*Y(ipro)
        end if
     endif

     ! (gamma,p)
     if(rr_tmp%reac_type.eq.rrt_gp) then
        nuc = rr_tmp%parts(1)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           gp_sum  = gp_sum  + rat*Y(nuc)
        end if
     endif

     ! (a,g)
     if(rr_tmp%reac_type.eq.rrt_ag) then
        nuc = rr_tmp%parts(2)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           ag_sum  = ag_sum  + rat*Y(nuc)*Y(ihe4)
        end if
     endif

     ! (g,a)
     if(rr_tmp%reac_type.eq.rrt_ga) then
        nuc = rr_tmp%parts(1)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           ga_sum  = ga_sum  + rat*Y(nuc)
        end if
     endif

     ! (n,p)
     if(rr_tmp%reac_type.eq.rrt_np) then
        nuc = rr_tmp%parts(2)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           np_sum  = np_sum  + rat*Y(nuc)*Y(ineu)
        end if
     endif

     ! (p,n)
     if(rr_tmp%reac_type.eq.rrt_pn) then
        nuc = rr_tmp%parts(2)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           pn_sum  = pn_sum  + rat*Y(nuc)*Y(ipro)
        end if
     endif

     ! (a,n)
     if(rr_tmp%reac_type.eq.rrt_an) then
        nuc = rr_tmp%parts(2)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           an_sum  = an_sum  + rat*Y(nuc)*Y(ihe4)
        end if
     endif

     ! (n,a)
     if(rr_tmp%reac_type.eq.rrt_na) then
        nuc = rr_tmp%parts(2)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           na_sum  = na_sum  + rat*Y(nuc)*Y(ineu)
        end if
     endif

  enddo

  fission:do i=1,nfiss
     fr_tmp = fissrate(i)
     rat = fr_tmp%cached
     ! neutron-induced fission
     if(fr_tmp%reac_type.eq.rrt_nf) then
        nuc = fr_tmp%fissparts(2)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           nfiss_sum = nfiss_sum + rat*Y(nuc)*Y(ineu)
        end if
     end if

     ! spontaneous fission
     if(fr_tmp%reac_type.eq.rrt_sf) then
        nuc = fr_tmp%fissparts(1)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           sfiss_sum = sfiss_sum + rat*Y(nuc)
        end if
     end if

     ! beta-delayed fission
     if(fr_tmp%reac_type.eq.rrt_bf) then
        nuc = fr_tmp%fissparts(1)
        if(rat.gt.1.d-20 .and. Y(nuc).gt.1.d-20) then
           bfiss_sum = bfiss_sum + rat*Y(nuc)
        end if
     end if

  end do fission

  tau_ng    = ysum / ng_sum
  tau_gn    = ysum / gn_sum
  tau_pg    = ysum / pg_sum
  tau_gp    = ysum / gp_sum
  tau_ag    = ysum / ag_sum
  tau_ga    = ysum / ga_sum
  tau_np    = ysum / np_sum
  tau_pn    = ysum / pn_sum
  tau_an    = ysum / an_sum
  tau_na    = ysum / na_sum
  tau_pa    = ysum / pa_sum
  tau_ap    = ysum / ap_sum
  tau_beta  = ysum / beta_sum
  tau_alpha = ysum / alpha_sum
  tau_nfiss = ysum / nfiss_sum
  tau_sfiss = ysum / sfiss_sum
  tau_bfiss = ysum / bfiss_sum

  if (.not. present(hdf5)) then
     ! write it to an ascii file
     write(tsfile, '(*(es23.16,3x))') &
        ctime, T9, tau_ga,tau_ag, tau_ng, tau_gn, tau_pg, tau_gp, tau_np, &
                   tau_pn, tau_an, tau_na,tau_pa,tau_ap, tau_beta, tau_alpha, &
                   tau_nfiss, tau_sfiss, tau_bfiss
  else
#ifdef USE_HDF5
     ! Write it to the hdf5
     call extend_av_timescales(ctime,tau_ga, tau_ag, tau_ng, tau_gn, tau_pg,   &
                                     tau_gp, tau_np, tau_pn, tau_an, tau_na,   &
                                     tau_ap, tau_pa, tau_beta, tau_alpha,  &
                                     tau_nfiss, tau_sfiss, tau_bfiss)

#else
     continue
#endif
  endif

  INFO_EXIT('calc_av_timescales')

end subroutine calc_av_timescales


!>
!! Compute final abundances and write final output
!!
!! There are three files written _finabsum.dat_
!! containing abundances summed over equal A,
!! _finab.dat_ containing all final abundances,
!! _finabelem.dat_ containing all abundances
!! summed over equal Z.
subroutine finab(Y)
  use global_class, only: isotope, net_size
  use file_handling_class
  implicit none

  real(r_kind),dimension(net_size),intent(in)  :: Y
  integer,dimension(net_size)                  :: list
  real(r_kind),dimension(:),allocatable        :: XA,YZ
  integer            :: finfile,finsfile,finzfile
  integer            :: mval,Zmval
  integer            :: i,k
  integer            :: mass,protonnum

  INFO_ENTRY('finab')

  finfile= open_outfile ('finab.dat')
  write(finfile,'(A)') "# 1:A 2:Z 3:N 4:Yi 5:Xi"

  call finab_sort(net_size,list)

  do k=1,net_size
     i=list(k)
     if (Y(i).le.1.d-20) cycle
     write(finfile,'(3(i5,2x),2es23.14)') isotope(i)%mass,   &
          isotope(i)%p_nr,isotope(i)%n_nr,Y(i),              &
          Y(i)*isotope(i)%mass
  end do
  call close_io_file(finfile,'finab.dat')

  mval = maxval(isotope%mass)
  allocate(XA(mval))
  Zmval = maxval(isotope%p_nr)
  allocate(YZ(Zmval))

  XA = 0.d0
  YZ = 0.d0

  do i=1,net_size
     mass=isotope(i)%mass
     protonnum = isotope(i)%p_nr
     if(Y(i).gt.1.d-20) then
        XA(mass) = XA(mass) + Y(i)*mass
        if (protonnum.gt.0) YZ(protonnum) = YZ(protonnum) + Y(i) ! omit neutrons
     end if
  end do

  finsfile= open_outfile ('finabsum.dat')
  write(finsfile,'(A)') "# 1:A 2:Y(A) 3:X(A)"
  do i=1,mval
     !i:Mass, XA(i)/i: Abundance, XA(i): Mass fraction
     write(finsfile,'(i5,2x,es23.14,2x,es23.14)')i, XA(i)/i ,XA(i)
  end do
  finzfile= open_outfile ('finabelem.dat')
  write(finzfile,'(A)') "# 1:Z 2:Y(Z)"
  do i=1,Zmval
     write(finzfile,'(i5,2x,es23.14)')i, YZ(i)
  end do

  call close_io_file(finsfile,'finabsum.dat')

  INFO_EXIT('finab')
end subroutine finab


!>
!! auxiliary sorting subroutine
!!
subroutine finab_sort(nsize,list)
  use global_class, only:isotope
  use file_handling_class

  implicit none

  type :: sorttype
     integer :: a
     integer :: z
     logical :: chk
  end type sorttype

  integer,intent(in)                        :: nsize
  integer,dimension(nsize),intent(out)      :: list
  type(sorttype),dimension(nsize)           :: tmp
  integer                                   :: i,j
  integer                                   :: amin,zmin
  integer                                   :: nxt

  INFO_ENTRY('finab_sort')

  amin = 10000
  zmin = 10000

  do i=1,nsize
     tmp(i)%a = isotope(i)%mass
     tmp(i)%z = isotope(i)%p_nr
     tmp(i)%chk = .false.
  end do

  do i=1,nsize
     do j=1,nsize
        if(tmp(j)%chk) cycle
        select case (tmp(j)%a-amin)
        case(:-1)
           nxt = j
           amin = tmp(j)%a
           zmin = tmp(j)%z
        case(0)
           if (tmp(j)%z.lt.zmin) then
              nxt = j
              amin = tmp(j)%a
              zmin = tmp(j)%z
           end if
        case default
           cycle
        end select
     end do
     list(i) = nxt
     tmp(nxt)%chk = .true.
     amin=10000
     zmin=10000
  end do

  INFO_EXIT('finab_sort')

  return

end subroutine finab_sort


!>
!! Close files at the end of the calculation
!!
!! @author: Moritz Reichert
!!
!! \b Edited:
!!     - 01.03.17
!!     - 25.01.21
!! .
subroutine analysis_finalize()
   use parameter_class, only: track_nuclei_every
   implicit none

   INFO_ENTRY('analysis_finalize')

   ! Close the timescale file
   if (timescales_every .gt. 0) then
      call close_io_file(tsfile,'timescales.dat')
   end if

   ! Close the Newton-Raphson diagnostic file
   if (nrdiag_every.gt.0) then
      call close_io_file(nrdiag_unit,"nrdiag.dat")
   end if

   ! Close track nuclei file
   if (track_nuclei_every .gt. 0) then
      call close_io_file(track_unit,"tracked_nuclei.dat")
   end if

   INFO_EXIT('analysis_finalize')

end subroutine analysis_finalize


end module analysis
