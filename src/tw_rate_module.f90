!> @file tw_rate_module.f90
!!
!! The error file code for this file is ***W44***.
!!
!! Contains the module \ref tw_rate_module.


!>
!! @brief This module contains everything for the theoretical weak rates that
!!        are added into the rate array
!!
!! The theoretical weak rates in the current form were implemented by M. Ugliano.
!! They contain electron captures and beta decays.
!! Rates that are contained in the weak rate file and in reaclib
!! get replaced by the ones contained in weak rates for small
!! temperatures (see parameter_class::temp_reload_exp_weak_rates).
!!
!! @see [Langanke & Martinez-Pinedo 2001](https://ui.adsabs.harvard.edu/abs/2001ADNDT..79....1L/abstract),
!!      [Fuller et al. 1985](https://ui.adsabs.harvard.edu/abs/1985ApJ...293....1F/abstract),
!!      [Oda et al. 1994](https://ui.adsabs.harvard.edu/abs/1994ADNDT..56..231O/abstract),
!!      [Weak rate library](https://groups.nscl.msu.edu/charge_exchange/weakrates.html)
#include "macros.h"
module tw_rate_module
   use global_class, only: reactionrate_type
   implicit none


   !> data fields for weak rates given in Langanke&Martinez-Pinedo 2001
   type,private                              :: weakrate_type
      integer,dimension(2)                   :: parts             !< isotope index of participants [1:2]
      logical                                :: is_ec             !< true if reaction is electron capture reaction
      character(4)                           :: source            !< source of reaction rate
      integer                                :: n_points          !< total amount of grid_points
      integer                                :: n_temp_grid       !< Length of temperature grid
      integer                                :: n_rho_grid        !< Length of density grid
      real(r_kind)                           :: q_value           !< reaction Q-value [MeV]
      real(r_kind)                           :: min_rho,max_rho   !< min and max rho*ye
      real(r_kind)                           :: min_temp,max_temp !< min and max temp
      real(r_kind),dimension(:),allocatable  :: temp_grid         !< temperature grid for tabulation
      real(r_kind),dimension(:),allocatable  :: rho_grid          !< rho*ye grid for tabulation
      real(r_kind),dimension(:,:),allocatable:: mue_kin           !< Electron chemical potential (only kinetic, no rest mass)
      real(r_kind),dimension(:,:),allocatable:: dmue_kin_dT       !< Temperature derivative of the chemical potential
      real(r_kind),dimension(:,:),allocatable:: dmue_kin_dR       !< Density derivative of the chemical potential
      real(r_kind),dimension(:,:),allocatable:: dmue_kin_dT_dR    !< Mixed partial derivative of the chemical potential
      real(r_kind),dimension(:,:),allocatable:: beta_rate         !< tabulated beta decay rates
      real(r_kind),dimension(:,:),allocatable:: dbeta_dT          !< Temperature derivative of the beta rate
      real(r_kind),dimension(:,:),allocatable:: dbeta_dR          !< Density derivative of the beta rate
      real(r_kind),dimension(:,:),allocatable:: dbeta_dT_dR       !< Mixed partial derivative of the beta rate
      real(r_kind),dimension(:,:),allocatable:: ft_rate           !< tabulated log<ft> electron/positron capture rates
      real(r_kind),dimension(:,:),allocatable:: dft_dT            !< Temperature derivative of the log<ft> rates
      real(r_kind),dimension(:,:),allocatable:: dft_dR            !< Density derivative of the log<ft>
      real(r_kind),dimension(:,:),allocatable:: dft_dT_dR         !< Mixed partial derivative of the log<ft>
      real(r_kind),dimension(:,:),allocatable:: nu_loss           !< tabulated energy loss due to (anti-)neutrino emission
   end type weakrate_type

   type(weakrate_type),dimension(:),allocatable,public         :: weak_rate !< array containing the weak reactions

   logical,public                                              :: weak         !< switch for weak rates
   integer,private                                             :: nweak        !< number of weak reactions taken from Langanke&Martinez-pinedo
   type(reactionrate_type), dimension(:), allocatable,private  :: rrate_weak   !< Reaction rate representative in global_class::rrate
   integer,dimension(2,2)   , private                          :: wk_index     !< Multi-index for the weak rates
   real(r_kind),public                                         :: mue          !< chemical potential of the electron
   integer,parameter,private                                   :: mue_ident=1  !< Identifier for chemical potential
   integer,parameter,private                                   :: beta_ident=2 !< Identifier for beta rate
   integer,parameter,private                                   :: ft_ident=3   !< Identifier for forward rate
   integer,parameter,private                                   :: nu_loss_ident=4   !< Identifier for the neutrino loss

   integer,private :: n_ec, n_o !< individual reaction types

    !
    ! Public and private fields and methods of the module
    !
    public:: &
        init_theoretical_weak_rates, weak_index, weak_inter, reload_exp_weak_rates, &
        calculate_twr_rate
    private:: &
        readweak_logft, read_theoretical_weak_rates,&
        sort, write_reac_verbose_out

contains


   !> Initialize theoretical weak rates
   !!
   !! This subroutine set the flag weak and calls
   !! the necessary subroutines to read the rates.
   !!
   !! @author Moritz Reichert
   !! @date 24.01.21
   subroutine init_theoretical_weak_rates()
      use parameter_class, only: iwformat
      implicit none

      ! Set the amount to zero, later it is set in readweak_logft
      nweak = 0
      ! Set default for individual reactions
      n_ec=0; n_o=0

      weak= (iwformat.ne.0)
      if (weak) call read_theoretical_weak_rates()
      if (weak) call initialize_cubic_interp()
      if (weak) call write_reac_verbose_out()
      if (weak) call output_n_p()
   end subroutine init_theoretical_weak_rates


   !> Calculate derivatives of weak reactions
   !!
   !! These derivatives are necessary for the bicubic interpolation.
   !!
   !! @author: M. Reichert
   subroutine initialize_cubic_interp
     use inter_module, only: calc_derivative_2D
     use parameter_class, only: iwinterp,use_timmes_mue
     implicit none
     integer :: nt,nr !< Number of temperature and rho points
     integer :: i !< Loop variable

     ! (Bi-)Cubic interpolation
     if (iwinterp .eq. 1) then
       do i=1, nweak
         nt = weak_rate(i)%n_temp_grid
         nr = weak_rate(i)%n_rho_grid
         ! Allocate arrays
         allocate(weak_rate(i)%dbeta_dT(nt,nr),&
                  weak_rate(i)%dbeta_dR(nt,nr),&
                  weak_rate(i)%dbeta_dT_dR(nt,nr),&
                  weak_rate(i)%dft_dT(nt,nr),&
                  weak_rate(i)%dft_dR(nt,nr),&
                  weak_rate(i)%dft_dT_dR(nt,nr) )

        ! Calculate derivatives
        call calc_derivative_2D(weak_rate(i)%beta_rate,weak_rate(i)%temp_grid,&
                               weak_rate(i)%rho_grid,weak_rate(i)%dbeta_dT,&
                               weak_rate(i)%dbeta_dR,weak_rate(i)%dbeta_dT_dR,&
                               nt,nr)
        call calc_derivative_2D(weak_rate(i)%ft_rate,weak_rate(i)%temp_grid,&
                               weak_rate(i)%rho_grid,weak_rate(i)%dft_dT,&
                               weak_rate(i)%dft_dR,weak_rate(i)%dft_dT_dR,&
                               nt,nr)
       end do
     end if

     ! The chemical potential should be interpolated bicubic
     if (.not. use_timmes_mue) then
       do i=1, nweak
         nt = weak_rate(i)%n_temp_grid
         nr = weak_rate(i)%n_rho_grid
         ! Allocate arrays
         allocate(weak_rate(i)%dmue_kin_dT(nt,nr),&
                  weak_rate(i)%dmue_kin_dR(nt,nr),&
                  weak_rate(i)%dmue_kin_dT_dR(nt,nr))
        ! Calculate derivatives
        call calc_derivative_2D(weak_rate(i)%mue_kin,weak_rate(i)%temp_grid,&
                               weak_rate(i)%rho_grid,weak_rate(i)%dmue_kin_dT,&
                               weak_rate(i)%dmue_kin_dR,weak_rate(i)%dmue_kin_dT_dR,&
                               nt,nr)
       end do
     end if

   end subroutine initialize_cubic_interp


   !> Write the amount of individual reactions to the out
   !!
   !! The rates are always counted, for a certain verbose level they
   !! are also printed to the OUT file
   !!
   !! @author M. Reichert
   !! @date 27.01.21
   subroutine write_reac_verbose_out()
      use error_msg_class, only: int_to_str, write_data_to_std_out
      implicit none
      character(len=7) :: tmp !< temporary character for pretty output

      if (VERBOSE_LEVEL .ge. 1) then
         call write_data_to_std_out("Amount weak (ffn) rates",int_to_str(nweak))
      elseif (VERBOSE_LEVEL .ge. 2) then
         if (nweak .gt. 0) write(*,"(A)") ""
         if (nweak .gt. 0) write(*,"(A)")"    FFN rates:  "
         if (nweak .gt. 0) write(*,"(A)")"   |---------------------------|"
         tmp = int_to_str(nweak)
         if (nweak .gt. 0) write(*,"(A)")"   | Total            :"//adjustr(tmp)//" |"
         tmp = int_to_str(n_ec)
         if (n_ec .gt. 0) write(*,"(A)") "   | Electron-capture :"//adjustr(tmp)//" |"
         tmp = int_to_str(n_o)
         if (n_o .gt. 0)  write(*,"(A)") "   | Beta-decay       :"//adjustr(tmp)//" |"
         if (nweak .gt. 0)write(*,"(A)") "   |---------------------------|"
         if (nweak .gt. 0)write(*,"(A)") ""
      end if

   end subroutine write_reac_verbose_out


   !> Merge theoretical weak rates into larger rate array.
   !!
   !! This merge is partly done instantaneously and partly when
   !! the temperature reaches parameter_class::reload_exp_weak_rates.
   subroutine merge_theoretical_weak_rates(rrate_array,rrate_length)
      use error_msg_class,  only: raise_exception
      use mergesort_module, only: rrate_ms,rrate_sort
      use global_class,     only: common_weak_rates, only_theo_weak_rates
      implicit none
      type(reactionrate_type),dimension(:),allocatable,intent(inout) :: rrate_array  !< Large rate array, containing all reactions
      integer,intent(inout)                                          :: rrate_length !< length of rrate_array
      type(reactionrate_type),dimension(:),allocatable               :: rrate_tmp    !< Temporary rate array
      integer                                                        :: alloc_stat   !< Allocation state
      integer                                                        :: new_length   !< New length of the array

      if (weak) then
         if (.not. allocated(rrate_array)) then
            !-- Allocate the reaclib rate array
            allocate(rrate_array(nweak),stat=alloc_stat)
            if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                       "merge_theoretical_weak_rates",440001)
            rrate_array(1:nweak) = rrate_weak(1:nweak)
            rrate_length = nweak
         else
            call rrate_sort(rrate_length,rrate_array(1:rrate_length))
            common_weak_rates     = 0
            only_theo_weak_rates  = 0
            !----- merge weak rates into rrate
            call rrate_ms(rrate_array(1:rrate_length),rrate_length,rrate_weak,&
                         size(rrate_weak),1,new_length,rrate_tmp)
            rrate_length = new_length
            ! Reallocate rrate_array with new length
            if (allocated(rrate_array)) deallocate(rrate_array)
            allocate(rrate_array(rrate_length),stat=alloc_stat)
            if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                       "merge_theoretical_weak_rates",440001)
            rrate_array(1:rrate_length) = rrate_tmp(1:rrate_length)
         end if
         !-- Deallocate the reaclib rate array
         deallocate(rrate_weak)
      end if

   end subroutine merge_theoretical_weak_rates


   !> Calculate the theoretical weak rate.
   !!
   !! For \ref parameter_class::iwformat equal to 1, this routine interpolates
   !! the rate on the temperature - density grid. Otherwise, for parameter_class::iwformat
   !! equal 2, this subroutine calculates the effective phase space integral to translate
   !! the tabulated log <ft> values to actual rates. For this also the chemical
   !! potential has to be known.
   !!
   !! \b Edited:
   !!    - 26.07.22, MR: Created this subroutine from previously existing code in Jacobi_init
   !! .
   !!
   !! @warning In case of \ref parameter_class::iwformat equal to two
   !!          and using timmes chemical potential (\ref parameter_class::use_timmes_mue)
   !!          \ref mue has to be calculated prior to calling this subroutine.
   subroutine calculate_twr_rate(rrate, temp, rho, Ye, rat_calc, nuloss)
     use parameter_class, only: iwformat, use_timmes_mue, unit, heating_frac
     use effphase_class,  only: effphase
     implicit none
     type(reactionrate_type),intent(inout):: rrate    !< rate instance
     real(r_kind),intent(in)              :: temp     !< Temperature [GK]
     real(r_kind),intent(in)              :: rho      !< Density [g/ccm]
     real(r_kind),intent(in)              :: Ye       !< Electron fraction
     real(r_kind),intent(out)             :: rat_calc !< rate value
     logical,optional                     :: nuloss   !< calculate neutrino loss?
     ! Internal variables
     integer                              :: wind     !< Weak index
     real(r_kind)                         :: rbeta,rft,repsi,rnu
     real(r_kind)                         :: epsi     !< effective phase-space integral
     logical                              :: is_ec    !< flag for electron captures

     ! Set default value of neutrino loss
     rrate%nu_frac = heating_frac

     ! weak_index determines the cornerpoints for the temp and dens interpolation
     wind = int(rrate%param(1))
     call weak_index (temp, rho, Ye, weak_rate(wind))
     is_ec = weak_rate(wind)%is_ec
     ! weak_inter interpolates the weak rate using the result from weak_index
     ! for weak rates rrate(i)%param gives the index of the rate in weak_rate array
     if (iwformat .eq. 1) then
        rft   = weak_inter(ft_ident,temp,dlog10(rho*Ye),weak_rate(wind))
        rbeta = weak_inter(beta_ident,temp,dlog10(rho*Ye),weak_rate(wind))
        if (present(nuloss)) then
           if (nuloss) then
              rnu   = weak_inter(nu_loss_ident,temp,dlog10(rho*Ye),weak_rate(wind))
              rnu   = rnu/(rft+rbeta)
           end if
        end if
        rat_calc = rft+rbeta
     else if (iwformat .eq. 2) then
        ! Use the electron chemical potential directly from the rate file
        if (.not. use_timmes_mue) then
          mue = weak_inter(mue_ident,temp,dlog10(rho*Ye),weak_rate(wind),.True.)
          mue = dlog10(mue)+ unit%mass_e
        end if
        rft   = weak_inter(ft_ident,temp,dlog10(rho*Ye),weak_rate(wind))
        rbeta = weak_inter(beta_ident,temp,dlog10(rho*Ye),weak_rate(wind))
        if (present(nuloss)) then
           if (nuloss) then
              rnu   = weak_inter(nu_loss_ident,temp,dlog10(rho*Ye),weak_rate(wind))
           end if
        end if
        if (is_ec) then
           call effphase(temp,mue,weak_rate(wind)%q_value,epsi)
        else
           call effphase(temp,-mue,weak_rate(wind)%q_value,epsi)
        end if
        if ((epsi.lt.1.d-100).or.(rft.le.9.d-99)) then
           repsi = 1.d-100
        else
           repsi = dlog(2.d0)*epsi/rft
        end if


        ! total rate = beta rate + ft-rate
        rat_calc = repsi + rbeta
      end if

      ! Calculate neutrino loss fractions
      ! Note that the fraction can be larger than one
      if (present(nuloss)) then
        if (nuloss) then
            ! This can also be negative, but later on it is multiplied by the
            ! q-value again
            rrate%nu_frac = rnu/rrate%q_value
        end if
      end if


   end subroutine calculate_twr_rate





   !> Reads and counts the amount of theoretical weak rates
   !!
   !! Furthermore, this routine creates the rrate_weak array that is
   !! merged into the total rate array rrate.
   !!
   !! \b Edited:
   !!   - 24.01.21, MR - copied it into this separate subroutine
   !! .
   subroutine read_theoretical_weak_rates()
      use parameter_class,  only: iwformat, weak_rates_file, unit, heating_frac
      use mergesort_module, only: rrate_sort
      use benam_class,      only: getcoefficients
      use file_handling_class
      implicit none

      integer  :: twr_unit  !< File ID of theoretical weak rate file
      integer  :: cnt       !< Count variable for the amount of rates
      integer  :: i         !< Loop variable


      if ((iwformat.eq.1) .or. (iwformat.eq.2)) then
        twr_unit= open_infile(weak_rates_file)

   !----- readweak(_logft) returns weak_rate(:) and number of weakrates
        if ((iwformat.eq.1) .or. (iwformat.eq.2)) call readweak_logft(twr_unit,cnt)

        ! Set the correct number of rates
        nweak = cnt

   !----- create rrate-type array of weak-rates, needed to merge them later
        allocate(rrate_weak(cnt))
        do i = 1,cnt
            rrate_weak(i)%group       = 1
            rrate_weak(i)%parts       = 0
            rrate_weak(i)%parts(1:2)  = weak_rate(i)%parts(1:2)
            rrate_weak(i)%source      = weak_rate(i)%source
            rrate_weak(i)%is_reverse  = .false.
            rrate_weak(i)%is_resonant = .false.
            rrate_weak(i)%is_weak     = .true.
            rrate_weak(i)%reac_src    = rrs_twr
            rrate_weak(i)%cached      = -1
            rrate_weak(i)%nu_frac     = heating_frac
   !----- Q-value of weakrates does not contain the electron rest mass
   !      here it is added to get full Q-values in rrate. The weak_rate
   !      Q-values are needed for the phase-space integral!
            if (weak_rate(i)%is_ec) then
               rrate_weak(i)%q_value = weak_rate(i)%q_value + unit%mass_e
               rrate_weak(i)%reac_type   = rrt_betp
               n_ec=n_ec+1
            else
               n_o=n_o+1
               rrate_weak(i)%q_value = weak_rate(i)%q_value - unit%mass_e
               rrate_weak(i)%reac_type   = rrt_betm
            end if
            rrate_weak(i)%param(1) = dble(i)
        end do
   !----- sort them(just to be sure)
        call rrate_sort(cnt,rrate_weak)

        ! Close the file again
        close(twr_unit)
        ! get the correct coefficients
        call getcoefficients(rrate_weak,nweak)
      end if

   end subroutine



   !>
   !! Reads weak reaction rates in log format
   !!
   !! This subroutine is called when parameter_class::iwformat was
   !! set to 1. It assumes weak reaction rates that are
   !! logarithmic. It will fill the global_class::weak_rate array
   !! with rates that are later indicated by the label " ffn". An example
   !! of the file that is read here is given by:
   !! \file{
   !! neg. daughter Sc45 z=21 n=24 a=45 Q=  0.7678
   !! pos. daughter Ca45 z=20 n=25 a=45 Q= -0.7678
   !!                      +++ Sc45 --> Ca45 +++      --- Ca45 --> Sc45 ---
   !!   t9   lrho    uf    lbeta+    lfte-     lnu    lbeta-    lfte+   lnubar
   !!  0.01  1.0  -0.003  -99.999  -99.999  -99.999   -7.302  -99.999   -0.746
   !!  0.10  1.0  -0.058  -91.517    5.796   -1.577   -7.302    7.565   -0.746
   !!  0.20  1.0  -0.134  -49.637    5.792   -1.266   -7.302    7.083   -0.746
   !!  ... }
   !!
   !! @remark Some rates in the file will replace rates from the reaclib, changing them there may not
   !!         have an effect as other rates are used then.
   !!
   !! @see parameter_class::temp_reload_exp_weak_rates, parameter_class::iwformat
   !!
   !! \b Edited:
   !!        - 12.01.14
   !!        - 11.03.22 (MR: extended routine for more general grids)
   !! .
   subroutine readweak_logft(sourcefile,cnt)
     use mergesort_module,only: bubblesort,reorder
     use global_class,    only: isotope
     use benam_class
     use error_msg_class

     implicit none

     integer, intent(in)                   :: sourcefile !< file to read weak rates from
     integer, intent(out)                  :: cnt        !< total count of weak rates
     !
     real(r_kind)                          :: qfwd, qbwd !< q-values of forward and backward rates
     real(r_kind)                          :: min_rho,max_rho,min_temp,max_temp
     character*4,dimension(2)              :: wnam
     character*5,dimension(2)              :: wparts
     integer,dimension(2)                  :: parts_index
     real(r_kind),dimension(:),allocatable :: beta_bwd,ec_ft,nu_loss,beta_fwd,mue_tmp
     real(r_kind),dimension(:),allocatable :: pc_ft,anu_loss,t_grid,r_grid
     integer,dimension(:),allocatable      :: sort_keys
     character(len=300)                    :: null
     integer                               :: i,j,k,read_stat,alloc_stat,len_rate
     integer                               :: count
     logical                               :: cycle_loop

     INFO_ENTRY("readweak_logft")

     k=0
     walloc_loop: do
   ! read weak rate from file
       read(sourcefile,101, iostat = read_stat) wnam(1),qfwd
       if (read_stat /= 0) exit
       read(sourcefile,101) wnam(2),qbwd
       read(sourcefile,*)
       read(sourcefile,*)null
       call convert(wnam,wparts)

        len_rate = 0
        len_rate_loop: do
           read(sourcefile,"(a)", iostat = read_stat)null
           if (trim(adjustl(null)) .eq. "end") exit len_rate_loop
           if (read_stat /= 0) exit len_rate_loop
           len_rate = len_rate+1
        end do len_rate_loop

        parts_index = 0
        do j=1,2
          parts_index(j) = benam(wparts(j))
          if (parts_index(j) .eq. 0) then
            cycle walloc_loop
          end if
        end do

        k=k+2
     end do walloc_loop
     nweak = k
     allocate(weak_rate(nweak),stat=alloc_stat)
     if ( alloc_stat /= 0) call raise_exception('Allocation of "weak_rate" failed',&
                                                "readweak_logft",440001)
     rewind(sourcefile)

    ! Count the grid length
    k=0
    wlength_loop: do
      ! read weak rate from file
      read(sourcefile,101, iostat = read_stat) wnam(1),qfwd
      if (read_stat /= 0) exit
      read(sourcefile,101) wnam(2),qbwd
      read(sourcefile,*)
      read(sourcefile,*)null
      call convert(wnam,wparts)

       len_rate = 0
       len_rate_loop2: do
          read(sourcefile,"(a)", iostat = read_stat)null
          if (trim(adjustl(null)) .eq. "end") exit len_rate_loop2
          if (read_stat /= 0) exit len_rate_loop2
          len_rate = len_rate+1
       end do len_rate_loop2

       parts_index = 0
       do j=1,2
         parts_index(j) = benam(wparts(j))
         if (parts_index(j) .eq. 0) then
           cycle wlength_loop
         end if
       end do

       weak_rate(k+1)%n_points = len_rate
       weak_rate(k+2)%n_points = len_rate

       k=k+2
    end do wlength_loop
    rewind(sourcefile)


     k=1
     wouter_loop: do
   ! read weak rate from file
        read(sourcefile,101, iostat = read_stat) wnam(1),qfwd
        if (read_stat /= 0) exit
        read(sourcefile,101) wnam(2),qbwd
        read(sourcefile,*)
        read(sourcefile,*)null
        call convert(wnam,wparts)

        cycle_loop = .False.
        parts_index = 0
        do j=1,2
          parts_index(j) = benam(wparts(j))
          if (parts_index(j) .eq. 0) then
            cycle_loop = .True.
          end if
        end do
        if (cycle_loop) then
          len_rate_loop3: do
             read(sourcefile,"(a)", iostat = read_stat)null
             if (trim(adjustl(null)) .eq. "end") exit len_rate_loop3
             if (read_stat /= 0) exit len_rate_loop3
          end do len_rate_loop3
          cycle wouter_loop
        end if

        ! Allocate everything.... a bit lengthy
        if (allocated(beta_bwd)) then
          if (size(beta_bwd) .ne. weak_rate(k)%n_points) then
            deallocate(beta_bwd,ec_ft,nu_loss,beta_fwd,pc_ft,anu_loss,t_grid,r_grid,sort_keys,mue_tmp,stat=read_stat)
            if (read_stat /= 0) call raise_exception("Deallocation of rate variables failed!",&
                                                     "readweak_logft",440002)

            allocate(beta_bwd(weak_rate(k)%n_points),ec_ft(weak_rate(k)%n_points),&
                     nu_loss(weak_rate(k)%n_points),beta_fwd(weak_rate(k)%n_points),&
                     pc_ft(weak_rate(k)%n_points),anu_loss(weak_rate(k)%n_points),&
                     t_grid(weak_rate(k)%n_points),r_grid(weak_rate(k)%n_points),&
                     sort_keys(weak_rate(k)%n_points),mue_tmp(weak_rate(k)%n_points),stat=read_stat)
            if (read_stat /= 0) call raise_exception("Allocation of rate variables failed!",&
                                                     "readweak_logft",440001)
          end if
        else
          allocate(beta_bwd(weak_rate(k)%n_points),ec_ft(weak_rate(k)%n_points),&
                   nu_loss(weak_rate(k)%n_points),beta_fwd(weak_rate(k)%n_points),&
                   pc_ft(weak_rate(k)%n_points),anu_loss(weak_rate(k)%n_points),&
                   t_grid(weak_rate(k)%n_points),r_grid(weak_rate(k)%n_points),&
                   sort_keys(weak_rate(k)%n_points),mue_tmp(weak_rate(k)%n_points),stat=read_stat)
           if (read_stat /= 0) call raise_exception("Allocation of rate variables failed!",&
                                                    "readweak_logft",440001)
        end if

        ! Get minimum and maximum grid values
        min_temp=huge(min_temp)
        max_temp=-1*huge(min_temp)
        min_rho =huge(min_temp)
        max_rho =-1*huge(min_temp)
        ! Read everything
        do i=1,weak_rate(k)%n_points
            read(sourcefile,*)t_grid(i),r_grid(i),mue_tmp(i),beta_bwd(i),ec_ft(i),nu_loss(i),beta_fwd(i),pc_ft(i),anu_loss(i)
            if (t_grid(i)<min_temp)min_temp=t_grid(i)
            if (t_grid(i)>max_temp)max_temp=t_grid(i)
            if (r_grid(i)<min_rho) min_rho =r_grid(i)
            if (r_grid(i)>max_rho) max_rho =r_grid(i)
        end do
        read(sourcefile,'(a1)')null


        ! Sort according to temperature
        call bubblesort(1,weak_rate(k)%n_points,t_grid,sort_keys)
        ! Count amount of temperature
        count=1
        do i=2,weak_rate(k)%n_points
          if (t_grid(i) .ne. t_grid(i-1)) count = count + 1
        end do
        weak_rate(k)%n_temp_grid = count
        allocate(weak_rate(k)%temp_grid(count))
        ! Get the unique temperature points
        weak_rate(k)%temp_grid(1) = t_grid(1)
        count=1
        do i=2,weak_rate(k)%n_points
          if (t_grid(i) .ne. t_grid(i-1)) then
            count = count + 1
            weak_rate(k)%temp_grid(count) = t_grid(i)
          end if
        end do
        ! Sort all values accordingly
        call reorder(r_grid,sort_keys,weak_rate(k)%n_points)
        call reorder(beta_bwd,sort_keys,weak_rate(k)%n_points)
        call reorder(ec_ft,sort_keys,weak_rate(k)%n_points)
        call reorder(nu_loss,sort_keys,weak_rate(k)%n_points)
        call reorder(beta_fwd,sort_keys,weak_rate(k)%n_points)
        call reorder(mue_tmp,sort_keys,weak_rate(k)%n_points)
        call reorder(pc_ft,sort_keys,weak_rate(k)%n_points)
        call reorder(anu_loss,sort_keys,weak_rate(k)%n_points)

        ! Now sort according to density
        call bubblesort(1,weak_rate(k)%n_points,r_grid,sort_keys)
        ! Count amount of density
        count=1
        do i=2,weak_rate(k)%n_points
          if (r_grid(i) .ne. r_grid(i-1)) count = count + 1
        end do
        weak_rate(k)%n_rho_grid = count
        allocate(weak_rate(k)%rho_grid(count))
        ! Get unique densty values
        weak_rate(k)%rho_grid(1) = r_grid(1)
        count=1
        do i=2,weak_rate(k)%n_points
          if (r_grid(i) .ne. r_grid(i-1)) then
            count = count + 1
            weak_rate(k)%rho_grid(count) = r_grid(i)
          end if
        end do

        call reorder(t_grid,sort_keys,weak_rate(k)%n_points)
        call reorder(beta_bwd,sort_keys,weak_rate(k)%n_points)
        call reorder(ec_ft,sort_keys,weak_rate(k)%n_points)
        call reorder(nu_loss,sort_keys,weak_rate(k)%n_points)
        call reorder(beta_fwd,sort_keys,weak_rate(k)%n_points)
        call reorder(mue_tmp,sort_keys,weak_rate(k)%n_points)
        call reorder(pc_ft,sort_keys,weak_rate(k)%n_points)
        call reorder(anu_loss,sort_keys,weak_rate(k)%n_points)

        ! Allocate the things in the weak rate type
        allocate(weak_rate(k)%beta_rate(weak_rate(k)%n_temp_grid,weak_rate(k)%n_rho_grid),&
                 weak_rate(k)%mue_kin(weak_rate(k)%n_temp_grid,weak_rate(k)%n_rho_grid),&
                 weak_rate(k)%ft_rate(weak_rate(k)%n_temp_grid,weak_rate(k)%n_rho_grid),&
                 weak_rate(k)%nu_loss(weak_rate(k)%n_temp_grid,weak_rate(k)%n_rho_grid),stat=read_stat )
        if (read_stat /= 0) call raise_exception("Allocation of rate variables failed!",&
                                                 "readweak_logft",440001)

        count=1

        if (weak_rate(k)%n_rho_grid*weak_rate(k)%n_temp_grid .ne. weak_rate(k)%n_points) then
          call raise_exception("The weak rate grid was not a regular grid! Check weak rate:"//&
                                isotope(parts_index(1))%name//" ---> "//isotope(parts_index(2))%name//&
                                NEW_LINE('A')//"Number of rho gridpoints  : "//int_to_str(weak_rate(k)%n_rho_grid)//&
                                NEW_LINE('A')//"Number of T gridpoints    : "//int_to_str(weak_rate(k)%n_temp_grid)//&
                                NEW_LINE('A')//"Total number of gridpoints: "//int_to_str(weak_rate(k)%n_points),&
                                "readweak_logft",440003)
        end if

        do j=1, weak_rate(k)%n_rho_grid
          do i=1, weak_rate(k)%n_temp_grid
            weak_rate(k)%beta_rate(i,j) = beta_fwd(count)
            weak_rate(k)%mue_kin(i,j)   = mue_tmp(count)
            weak_rate(k)%ft_rate(i,j)   = pc_ft(count)
            weak_rate(k)%nu_loss(i,j)   = anu_loss(count)
            count = count + 1
          end do
        end do
        weak_rate(k)%parts(1)  = parts_index(2)
        weak_rate(k)%parts(2)  = parts_index(1)
        weak_rate(k)%q_value   = qfwd
        weak_rate(k)%source    = ' ffn'
        weak_rate(k)%min_temp  = min_temp
        weak_rate(k)%max_temp  = max_temp
        weak_rate(k)%min_rho   = min_rho
        weak_rate(k)%max_rho   = max_rho
        weak_rate(k)%is_ec     = .false.

        k=k+1
        ! Assign same grid length
        weak_rate(k)%n_temp_grid = weak_rate(k-1)%n_temp_grid
        weak_rate(k)%n_rho_grid  = weak_rate(k-1)%n_rho_grid
        ! And the same grid
        allocate(weak_rate(k)%temp_grid(weak_rate(k)%n_temp_grid ))
        allocate(weak_rate(k)%rho_grid(weak_rate(k)%n_rho_grid ))
        weak_rate(k)%temp_grid = weak_rate(k-1)%temp_grid
        weak_rate(k)%rho_grid = weak_rate(k-1)%rho_grid

        ! Allocate the things in the weak rate type
        allocate(weak_rate(k)%beta_rate(weak_rate(k)%n_temp_grid,weak_rate(k)%n_rho_grid),&
                 weak_rate(k)%mue_kin(weak_rate(k)%n_temp_grid,weak_rate(k)%n_rho_grid),&
                 weak_rate(k)%ft_rate(weak_rate(k)%n_temp_grid,weak_rate(k)%n_rho_grid),&
                 weak_rate(k)%nu_loss(weak_rate(k)%n_temp_grid,weak_rate(k)%n_rho_grid),stat=read_stat )
        if (read_stat /= 0) call raise_exception("Allocation of rate variables failed!",&
                                                 "readweak_logft",440001)

        count=1
        do j=1, weak_rate(k)%n_rho_grid
          do i=1, weak_rate(k)%n_temp_grid
            weak_rate(k)%beta_rate(i,j) = beta_bwd(count)
            weak_rate(k)%mue_kin(i,j)   = mue_tmp(count)
            weak_rate(k)%ft_rate(i,j)   = ec_ft(count)
            weak_rate(k)%nu_loss(i,j)   = nu_loss(count)
            count = count + 1
          end do
        end do

        weak_rate(k)%parts(1)  = parts_index(1)
        weak_rate(k)%parts(2)  = parts_index(2)
        weak_rate(k)%q_value   = qbwd
        weak_rate(k)%source    = ' ffn'
        weak_rate(k)%min_temp  = min_temp
        weak_rate(k)%max_temp  = max_temp
        weak_rate(k)%min_rho   = min_rho
        weak_rate(k)%max_rho   = max_rho
        weak_rate(k)%is_ec     = .true.
        cnt=k
        k=k+1

     end do wouter_loop
     call sort(cnt)

     INFO_EXIT("readweak_logft")

   101 format(t16,a4,t38,f8.4)

   end subroutine readweak_logft





   !>
   !! Set wk_index for given T, @f$\rho@f$ and @f$Y_e@f$
   !! The index is written into module variable \ref wk_index
   !!
   !! After the subroutine has been called the correct indices for a two
   !! dimensional interpolation are stored in \ref wk_index.
   !!
   !! @see weak_inter, readweak_logft, parameter_class::iwinterp, inter_module::get_indice_2D
   !!
   !! \b Edited:
   !!     - MR : 19.01.2021 - Removed hard coded temperatures
   !!     - MR : 11.03.2022 - Rates are now on two dimensional grid, rewrote this routine
   !!     - MR : 06.07.2022 - Moved things to inter_module, now it is only an interface
   !! .
   subroutine weak_index (ltemp, rho, Ye, wrt)
     use inter_module, only:get_indice_2D
      implicit none
      real(r_kind),intent(in)        :: ltemp         !< Temperature [GK]
      real(r_kind),intent(in)        :: rho           !< Density [g/ccm]
      real(r_kind),intent(in)        :: Ye            !< Electron fraction
      type(weakrate_type),intent(in) :: wrt           !< Input weak rate type
      real(r_kind)                   :: lrho          !< Logarithm of the \f$ \rho \cdot Y_e \f$

      lrho = dlog10(rho*Ye)
      call get_indice_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,&
                         wrt%n_temp_grid,wrt%n_rho_grid,wk_index)

   end subroutine weak_index


   !> Interpolation interface.
   !!
   !! This function interpolates the weak rates either bilinear or bicubic
   !!
   !! @author: M. Reichert
   !! @date: 14.03.22
   !!
   !! @warning The bicubic interpolation fails at values where -99 is neighboring.
   !!          There, it produces artifacts. For mue, it is, however, okay.
   !!
   !! \b Edited:
   !!     - MR : 06.07.2022 - Moved things to interpolation module
   !! .
   function weak_inter(rate,ltemp,lrho,wrt,force_cubic) result (wr)
     use parameter_class, only: iwinterp
     use error_msg_class, only: raise_exception, int_to_str
     use inter_module,    only: cubinter_2D, lininter_2D
     implicit none
     real(r_kind),intent(in)               :: lrho        !< Log(\f$ \rho \cdot Y_e \f$)
     real(r_kind),intent(in)               :: ltemp       !< Temperature [GK]
     logical,intent(in),optional           :: force_cubic !< Whether or not a cubic interpolation should be applied
     integer                               :: rate        !< Index for rate. 1: ft_rate, 2: beta_rate, 3: mue, 4: nu_loss
     type(weakrate_type),intent(in)        :: wrt         !< Input weak rate type
     real(r_kind)                          :: wr          !< Weak reaction rate, interpolated at the desired temp and dens

     ! 0: bilinear interpolation, 1: bicubic interpolation
     if (.not. present(force_cubic)) then
       select case(iwinterp)
        case(0)
          if (rate .eq. ft_ident) then
            wr = lininter_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,wrt%ft_rate,&
                             wrt%n_temp_grid,wrt%n_rho_grid,wk_index,2)
          elseif (rate .eq. beta_ident) then
            wr = lininter_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,wrt%beta_rate,&
                             wrt%n_temp_grid,wrt%n_rho_grid,wk_index,2)
          elseif (rate .eq. mue_ident) then
            wr = lininter_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,wrt%mue_kin,&
                             wrt%n_temp_grid,wrt%n_rho_grid,wk_index,2)
          elseif (rate .eq. nu_loss_ident) then
            wr = lininter_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,wrt%nu_loss,&
                             wrt%n_temp_grid,wrt%n_rho_grid,wk_index,2)
          else
            call raise_exception("Unknown array to interpolate, got "//int_to_str(rate),&
                                 "weak_inter",440009)
          end if
        case(1)
          if (rate .eq. ft_ident) then
            wr = cubinter_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,wrt%ft_rate,&
                             wrt%dft_dT,wrt%dft_dR,wrt%dft_dT_dR,&
                             wrt%n_temp_grid,wrt%n_rho_grid,wk_index,2)
          elseif (rate .eq. beta_ident) then
            wr = cubinter_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,wrt%beta_rate,&
                             wrt%dbeta_dT,wrt%dbeta_dR,wrt%dbeta_dT_dR,&
                             wrt%n_temp_grid,wrt%n_rho_grid,wk_index,2)
          elseif (rate .eq. mue_ident) then
            wr = cubinter_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,wrt%mue_kin,&
                             wrt%dmue_kin_dT,wrt%dmue_kin_dR,wrt%dmue_kin_dT_dR,&
                             wrt%n_temp_grid,wrt%n_rho_grid,wk_index,2)
          elseif (rate .eq. nu_loss_ident) then
            call raise_exception("Cubic interpolation not implemented for nu_loss",&
                                 "weak_inter") ! TODO implement error code
          end if
        case default
          call raise_exception("Unknown iwinterp type, got "//int_to_str(iwinterp),&
                               "weak_inter",440004)
        end select
      else
        if (rate .eq. ft_ident) then
          wr = cubinter_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,wrt%ft_rate,&
                            wrt%dft_dT,wrt%dft_dR,wrt%dft_dT_dR,&
                            wrt%n_temp_grid,wrt%n_rho_grid,wk_index,2)
        elseif (rate .eq. beta_ident) then
          wr = cubinter_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,wrt%beta_rate,&
                            wrt%dbeta_dT,wrt%dbeta_dR,wrt%dbeta_dT_dR,&
                            wrt%n_temp_grid,wrt%n_rho_grid,wk_index,2)
        elseif (rate .eq. mue_ident) then
          wr = cubinter_2D(ltemp,lrho,wrt%temp_grid,wrt%rho_grid,wrt%mue_kin,&
                            wrt%dmue_kin_dT,wrt%dmue_kin_dR,wrt%dmue_kin_dT_dR,&
                            wrt%n_temp_grid,wrt%n_rho_grid,wk_index,2)
        elseif (rate .eq. nu_loss_ident) then
            call raise_exception("Cubic interpolation not implemented for nu_loss",&
                                 "weak_inter") ! TODO implement error code
            end if
      end if

      ! Values are in logarithmic scale, convert back
      if (wr .lt. -30) then
         wr = 0.d0
      else
         wr = 1.d1**(wr)
      end if

      ! Check that it is not NaN
      if (wr .ne. wr) call raise_exception("Value got NaN in interpolation.",&
                                           "weak_inter")

      return
   end function weak_inter


   !>
   !! This routine is used to replace theoretical weak rates with experimental ones
   !! it's called when T9 < parameter_class::temp_reload_exp_weak_rates (default 1d-2)
   !!
   !! @remark Some rates will replace rates from the reaclib, changing them there may not
   !!         have an effect for low temperatures as other rates are used then.
   !!
   !! @see parameter_class::temp_reload_exp_weak_rates, parameter_class::iwformat
   !!      weak_index, weak_inter, temp_grid_weak
   !!
   !! \b Edited:
   !!     - 27.01.21, M.R.
   !!     - 23.03.23, M.R.
   !! .
   subroutine reload_exp_weak_rates()

     use parameter_class
     use global_class
     use benam_class
     use format_class
     use mergesort_module

     implicit none
     integer :: i, j, k, elements_still_to_replace,elements_still_to_remove
     logical :: replaced

     INFO_ENTRY("reload_exp_weak_rates")

     elements_still_to_replace = common_weak_rates
     elements_still_to_remove  = only_theo_weak_rates

     i = 1
     do while ((elements_still_to_replace .gt. 0).or.(elements_still_to_remove .gt. 0))
        if (rrate(i)%reac_src .eq. rrs_twr) then
           replaced = .false.
           do j = 1, common_weak_rates
              if ((rrate(i)%parts(1) .eq. rrate_weak_exp(j)%parts(1)) .and. &
                   (rrate(i)%parts(2) .eq. rrate_weak_exp(j)%parts(2))) then

                 rrate(i)%source      = rrate_weak_exp(j)%source
                 rrate(i)%q_value     = rrate_weak_exp(j)%q_value
                 rrate(i)%param       = rrate_weak_exp(j)%param
                 rrate(i)%reac_src    = rrate_weak_exp(j)%reac_src
                 rrate(i)%is_weak     = .true.
                 rrate(i)%reac_type   = rrate_weak_exp(j)%reac_type
                 rrate(i)%is_resonant = rrate_weak_exp(j)%is_resonant
                 rrate(i)%is_reverse  = rrate_weak_exp(j)%is_reverse
                 rrate(i)%is_const    = rrate_weak_exp(j)%is_const
                 rrate(i)%nu_frac     = rrate_weak_exp(j)%nu_frac
                 rrate(i)%ch_amount(1) = -1.d0
                 do k=2,6
                    if (rrate(i)%parts(k).ne.0) rrate(i)%ch_amount(k) = 1.d0
                 end do
                 elements_still_to_replace = elements_still_to_replace - 1
                 replaced = .true.
                 exit
              end if
           end do
           if (.not.replaced) then
             ! Put empty rate in order to not shift indices
             rrate(i)%is_reverse  = .false.
             rrate(i)%is_resonant = .false.
             rrate(i)%is_weak     = .true.
             rrate(i)%reac_type   = rrt_o
             rrate(i)%reac_src    = rrs_reacl
             rrate(i)%ch_amount(:)= 0.d0
             rrate(i)%param(:)    = 0
             rrate(i)%param(1)    = -99
             rrate(i)%nu_frac     = 0
             !MR: changed this to nreac -1 to stay in array bounds
             ! The following should not happen as it shifts indices that
             ! are needed by other rates
             ! (test with -check bounds compiler flag)
             !  do j = i, nreac-1
             !      rrate(j) = rrate(j+1)
             !   end do
             !   i = i-1
             elements_still_to_remove = elements_still_to_remove - 1
           end if
        end if
        i = i+1
        if (i .gt. nreac) exit
     end do

    ! As empty rates are introduced now, this is not necessary.
    !  nreac = nreac - only_theo_weak_rates
     deallocate(rrate_weak_exp)

     if (elements_still_to_replace .gt. 0) call raise_exception("Couldn't replace all weak rates!",&
                                                                "reload_exp_weak_rates",&
                                                                440005)
     if (elements_still_to_replace .lt. 0) call raise_exception("Replaced too many rates!",&
                                                                "reload_exp_weak_rates",&
                                                                440006)
     if (elements_still_to_remove .gt. 0)  call raise_exception("Couldn't remove all weak rates!",&
                                                                "reload_exp_weak_rates",&
                                                                440007)
     if (elements_still_to_remove .lt. 0)  call raise_exception("Removed too many rates!",&
                                                                "reload_exp_weak_rates",&
                                                                440008)
     INFO_EXIT("reload_exp_weak_rates")

     return

   end subroutine reload_exp_weak_rates


   !< Debug routine to output n <-> p reactions
   !!
   !! This routine creates the file "np_twr.dat", containing a variety of
   !! interpolated <ft> values. This can be useful to test the implemented
   !! interpolations.
   !!
   !! @author: M. Reichert
   !! @date: 05.07.22
   subroutine output_n_p
     use global_class, only: ineu,ipro
     use file_handling_class, only: open_outfile,close_io_file
     implicit none
     integer                   :: i,j                                !< Loop variable
     integer                   :: ind_n_p                            !< Index of the rate n -> p
     integer                   :: ind_p_n                            !< Index of the rate p -> n
     integer                   :: funit                              !< Unit for weak rates
     real(r_kind)              :: temp_min,temp_max,temp_incr        !< Minimum, Maximum, and increment of temperature for output
     real(r_kind)              :: rho_min ,rho_max ,rho_incr,rholog  !< Minimum, Maximum, and increment of electron density for output
     real(r_kind)              :: temp_tmp,rho_tmp,ye_tmp
     real(r_kind)              :: rat_pc,rat_ec,rat_bm,rat_bp

     if (VERBOSE_LEVEL .ge. 2) then
       temp_min = 0.0; temp_max = 101; temp_incr = 0.01
       rho_min  = 0.5; rho_max  = 10.5; rho_incr  = 0.25
       ind_n_p = -1
       ind_p_n = -1
       do i=1, nweak
         if ((weak_rate(i)%parts(1) .eq. ineu) .and. (weak_rate(i)%parts(2) .eq. ipro)) then
           ind_n_p= i
         end if
         if ((weak_rate(i)%parts(2) .eq. ineu) .and. (weak_rate(i)%parts(1) .eq. ipro)) then
           ind_p_n= i
         end if
       end do

       ! Only continue if rates are included
       if ((ind_n_p .eq. -1) .or. (ind_p_n .eq. -1)) return
       ! Open the file for debugging
       funit= open_outfile("np_twr.dat")
       write(funit,*) "#Temp [GK]   rho*ye  bp  ec  bm  pc wk11 wk12 wk21 wk22"
       ye_tmp   = 1
       temp_tmp = temp_min
       do i=1, int((temp_max-temp_min)/temp_incr)
         rho_tmp  = rho_min
         do j=1, int((rho_max-rho_min)/rho_incr)
           rholog = 10**rho_tmp

           call weak_index (temp_tmp, rholog, ye_tmp, weak_rate(ind_n_p))
           rat_pc= weak_inter(ft_ident,temp_tmp,rho_tmp,weak_rate(ind_n_p))
           rat_ec= weak_inter(ft_ident,temp_tmp,rho_tmp,weak_rate(ind_p_n))
           rat_bm= weak_inter(beta_ident,temp_tmp,rho_tmp,weak_rate(ind_n_p))
           rat_bp= weak_inter(beta_ident,temp_tmp,rho_tmp,weak_rate(ind_p_n))

           write(funit,"((6es12.5,4I4))") temp_tmp,rho_tmp,rat_bp,rat_ec,rat_bm,rat_pc,wk_index(1,1),wk_index(1,2),wk_index(2,1),wk_index(2,2)
           rho_tmp = rho_tmp+rho_incr
         end do
         temp_tmp = temp_tmp+temp_incr
       end do
       call close_io_file(funit, "np_twr.dat")
     end if
   end subroutine output_n_p


   !----------------------------Auxiliary helper functions---------------------!
   !---------------------------------------------------------------------------!




   !> Sorts the entries in weak_rate in increasing order of the decaying nuclide
   !!
   subroutine sort(cnt)

     use error_msg_class
     implicit none
     integer, intent(in)   :: cnt
     integer :: min_index, alloc_stat
     integer :: i,j,k
     type(weakrate_type),dimension(:),allocatable :: weak_temp

     allocate(weak_temp(nweak),stat=alloc_stat)
     if(alloc_stat /= 0) call raise_exception('Allocation of "weak_temp" failed',&
                                              "sort",440001)

     do i = 1,cnt
        do k=1,cnt
           if(weak_rate(k)%parts(1).gt.0) min_index = k
        end do
        do j=1,cnt
           if ((weak_rate(j)%parts(1).gt.0).and.                                   &
                (weak_rate(j)%parts(1).eq.weak_rate(min_index)%parts(1)).and.         &
                (weak_rate(j)%parts(2).lt.weak_rate(min_index)%parts(2))) then
              min_index = j
           else if((weak_rate(j)%parts(1).gt.0).and.                                    &
                (weak_rate(j)%parts(1).lt.weak_rate(min_index)%parts(1))) then
              min_index = j
           end if
        end do
        weak_temp(i) = weak_rate(min_index)
        weak_rate(min_index)%parts(1) = -1
     end do

     weak_rate = weak_temp
     deallocate(weak_temp)

   end subroutine sort


end module tw_rate_module
