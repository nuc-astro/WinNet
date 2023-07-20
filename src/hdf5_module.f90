!> @file hdf5_module.f90
!! The error file code for this file is ***W24***.
!! Contains the module \ref hdf5_module.

!> @brief Module that contains necessary subroutines to write
!! abundances to an hdf5 file.
!!
!! It is only included when the compilerflag "USE_HDF5" is set.
!! For more information about the hdf5 format
!! see [What is HDF5?](https://support.hdfgroup.org/HDF5/doc/H5.intro.html).
!! There is a tradeoff between having the data in the memory
!! and making an IO operation. Therefore, each output is first stored in a
!! temporary array of a specific size (or chunksize)
!! to prevent too much IO going on.
!!
!! @remark It is recommended, that the chunk sizes should
!!         be around the size of 1MB, however, here are written many 1d real(double)
!!         arrays and 1MB would be way more than storing the whole simulation temporary.
!!         Therefore the chunk sizes are smaller (in the order of a couple of hundreds).
!!
!! @author Moritz Reichert
!! @date 28.01.21
#include "macros.h"
module hdf5_module
! Only include the module with correct compiler
! Note that this is bad practice as it makes the tests
! more difficult, but it has the advantage that the code can still
! work even when no hdf5 is installed and used.
#ifdef USE_HDF5

   use hdf5
   use error_msg_class
   use global_class, only: net_size
   implicit none

   ! File related variables
   integer(HID_T),private    :: file_id                        !< File identifier
   character(len=20),private :: hdf5_filename="WinNet_data.h5" !< Name of the hdf5 file

   ! Flag that indicates if there will be an hdf5 output
   logical,private           :: hdf5_output=.False.

   character(LEN=*), parameter,private :: dsetname_u = "Units" !< Name of the units dataset

   ! General helper variables
   integer(HSIZE_T), dimension(1),private   :: dims_1d=(/1/) !< Helper variable to specify 1 dimension
   integer(HID_T),private                   :: dataspace     !< Dataspace identifier
   integer(HID_T),private                   :: crp_list      !< Dataset creation property identifier
   integer(HSIZE_T), dimension(1:2),private :: maxdims       !< Maximum dimensions (later set to unlimited)


   !-------- Snapshots variables -------------!
   INTEGER(HID_T),private  :: snaps_group_id !< The ID of the group
   integer,private         :: iter_snaps = 0 !< Iteration count,
                                             !< how often was already
                                             !< extended to the datasets?
   integer,private           :: chunk_counter_snaps = 0 !< How full is the chunk already?
   integer,parameter,private :: chunk_size_snaps=200    !< Chunk size
   real(r_kind),private,dimension(chunk_size_snaps)  :: chunk_store_snaps_time !< Storage for the
                                                                               !< chunk which is later
                                                                               !< written to the file.
   real(r_kind),private,dimension(:,:),allocatable :: chunk_store_snaps_Y !< Storage for the
                                                                          !< chunk which is later
                                                                          !< written to the file.

   ! Abundances
   integer(HID_T),private                   :: dset_id_Y                 !< Dataset identifier
   character(len=*), parameter,private      :: dsetname_Y = "Y"          !< Name of the dataset
   integer(HSIZE_T), dimension(1:2),private :: dims_Y                    !< Array dimensions (net_size)
   ! Time
   character(len=*), parameter,private  :: snaps_dsetname_t = "time" !< Name of the dataset
   integer(HID_T),private               :: snaps_dset_id_t           !< Dataset identifier


   !---------- Mainout variables -------------!
   INTEGER(HID_T),private  :: mainout_group_id !< The ID of the group
   integer,private         :: iter_mainout = 0 !< Iteration count,
                                               !< how often was already
                                               !< extended to the datasets?
   integer,private           :: chunk_counter_mainout = 0 !< How full is the chunk already?
   integer,parameter,private :: chunk_size_mainout=500    !< Chunk size
   real(r_kind),private,dimension(chunk_size_mainout,13)  :: chunk_store_mainout !< Storage for the
                                                                                 !< chunk which is later
                                                                                 !< written to the file.
   integer,private,dimension(chunk_size_mainout)  :: chunk_store_int_mainout !< Storage for integers in chunk
                                                                             !< which is later
                                                                             !< written to the file.
   ! Iteration
   character(len=*), parameter,private  :: mainout_dsetname_iter = "iteration" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_iter                !< Dataset identifier
   ! Time
   character(len=*), parameter,private  :: mainout_dsetname_t = "time" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_t           !< Dataset identifier
   ! Temperature
   character(len=*), parameter,private  :: mainout_dsetname_temp = "temp" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_temp           !< Dataset identifier
   ! Density
   character(len=*), parameter,private  :: mainout_dsetname_dens = "dens" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_dens           !< Dataset identifier
   ! Entropy
   character(len=*), parameter,private  :: mainout_dsetname_entr = "entr" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_entr           !< Dataset identifier
   ! Radius
   character(len=*), parameter,private  :: mainout_dsetname_rad = "rad"   !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_rad            !< Dataset identifier
   ! Electron fraction
   character(len=*), parameter,private  :: mainout_dsetname_ye = "ye" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_ye         !< Dataset identifier
   ! Neutron abundance
   character(len=*), parameter,private  :: mainout_dsetname_yn = "yn" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_yn         !< Dataset identifier
   ! Proton abundance
   character(len=*), parameter,private  :: mainout_dsetname_yp = "yp" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_yp         !< Dataset identifier
   ! Alpha abundance
   character(len=*), parameter,private  :: mainout_dsetname_ya = "ya" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_ya         !< Dataset identifier
   ! Ylight
   character(len=*), parameter,private  :: mainout_dsetname_ylight = "ylight" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_ylight             !< Dataset identifier
   ! Yheavy
   character(len=*), parameter,private  :: mainout_dsetname_yheavy = "yheavy" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_yheavy             !< Dataset identifier
   ! abar
   character(len=*), parameter,private  :: mainout_dsetname_abar = "abar" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_abar           !< Dataset identifier
   ! zbar
   character(len=*), parameter,private  :: mainout_dsetname_zbar = "zbar" !< Name of the dataset
   integer(HID_T),private               :: mainout_dset_id_zbar           !< Dataset identifier


   !---------- Track nuclei variables -------------!
   INTEGER(HID_T),private  :: track_group_id   !< The ID of the group
   integer,private         :: iter_tracked = 0 !< Iteration count,
                                               !< how often was already
                                               !< extended to the datasets?
   integer,private           :: chunk_counter_track = 0 !< How full is the chunk already?
   integer,parameter,private :: chunk_size_track=500    !< Chunk size
   real(r_kind),private,dimension(:,:),allocatable  :: chunk_store_track !< Temporary storage for
                                                                         !< tracked nuclei

   ! Time
   character(len=*), parameter,private                 :: tracked_dsetname_t = "time" !< Name of the dataset
   integer(HID_T),private                              :: tracked_dset_id_t           !< Dataset identifier
   ! Other
   integer(HSIZE_T), dimension(1:2),private 		   :: dims_track      !< Array dimensions
   character(len=5),dimension(:),allocatable, private  :: tracked_names   !< Names of the tracked nuclei
   integer(HID_T), private                             :: tracked_abu_id  !< Dataset identifier
   character(len=*), parameter,private                 :: tracked_dsetname_Y = "Y" !< Name of the dataset

   !---------- Timescale variables -------------!
   integer(HID_T),private  :: ts_group_id !< The ID of the group
   integer,private         :: iter_ts = 0 !< Iteration count,
                                          !< how often was already
                                          !< extended to the datasets?
   integer,private           :: chunk_counter_ts = 0 !< How full is the chunk already?
   integer,parameter,private :: chunk_size_ts=500    !< Chunk size
   real(r_kind),private,dimension(chunk_size_ts,18)  :: chunk_store_ts !< Temporary storage for
                                                                       !< average timescales

   ! Time
   character(len=*), parameter,private      :: ts_dsetname_t = "time" !< Name of the dataset
   integer(HID_T),private                   :: ts_dset_id_t           !< Dataset identifier
   ! Other
   character(len=9),dimension(17), private  :: ts_names = &
                 (/"tau_ga   ","tau_ag   ","tau_ng   ","tau_gn   ","tau_pg   ","tau_gp   ","tau_np   ",&
                   "tau_pn   ","tau_an   ","tau_na   ","tau_ap   ","tau_pa   ","tau_beta ","tau_alpha", &
                   "tau_nfiss","tau_sfiss","tau_bfiss"/)
                   !< Names of the timescales
   integer(HID_T),dimension(17), private    :: ts_reac_ids !< Dataset identifier


   !---------- Energy variables -------------!
   integer(HID_T),private  :: en_group_id !< The ID of the group
   integer,private         :: iter_en = 0 !< Iteration count,
                                          !< how often was already
                                          !< extended to the datasets?
   integer,private           :: chunk_counter_en = 0 !< How full is the chunk already?
   integer,parameter,private :: chunk_size_en=500    !< Chunk size
   real(r_kind),private,dimension(chunk_size_en,11)  :: chunk_store_en !< Temporary storage for
                                                                       !< energy generation
   real(r_kind),private,dimension(:,:),allocatable :: & !< Storage for the chunk which is later written to the file.
                                chunk_store_en_bet, chunk_store_en_ag,    chunk_store_en_pg, &
                                chunk_store_en_ng , chunk_store_en_pn,    chunk_store_en_ap, &
                                chunk_store_en_an , chunk_store_en_other, chunk_store_en_fiss, &
                                chunk_store_en_tot

   ! Time
   character(len=*), parameter,private      :: en_dsetname_t = "time" !< Name of the dataset
   integer(HID_T),private                   :: en_dset_id_t           !< Dataset identifier
   ! Other
   character(len=11),dimension(10), private  :: en_names = &
                 (/"engen_tot  ","S_src      ","engen_ng_gn","engen_pg_gp","engen_ag_ga",&
                   "engen_np_pn","engen_an_na","engen_ap_pa","engen_beta ", &
                   "engen_fiss "/)
                   !< Name of the energy generations
   integer(HID_T),dimension(10), private    :: en_reac_ids !< Dataset identifier

   character(len=14), private  :: & !< Name of the decay energy splitted by parent
                                 en_det_bet_name="detailed decay", en_det_ag_name   ="detailed (a,g)",   en_det_pg_name="detailed (p,g)", &
                                 en_det_ng_name ="detailed (n,g)", en_det_pn_name   ="detailed (p,n)",   en_det_ap_name="detailed (a,p)", &
                                 en_det_an_name ="detailed (a,n)", en_det_other_name="detailed other", en_det_fiss_name="detailed fiss ", &
                                 en_det_tot_name="detailed total"

   integer(HID_T), private     :: &                              !< Dataset identifier
                  en_det_bet_id,en_det_ag_id,    en_det_pg_id, &
                  en_det_ng_id, en_det_pn_id,    en_det_ap_id, &
                  en_det_an_id, en_det_other_id, en_det_fiss_id,&
                  en_det_tot_id

   !----------- Flow variables --------------!
   integer(HID_T),private  :: flow_group_id !< The ID of the group
   integer,private         :: iter_flow = 0 !< Iteration count,
                                            !< how often was already
                                            !< extended to the datasets?
   ! Time
   character(len=*), parameter,private      :: flow_dsetname_t = "time" !< Name of the dataset
   integer(HID_T),private                   :: flow_dset_id_t           !< Dataset identifier


   !----------- Nu gain/loss variables --------------!
   integer(HID_T),private  :: nuloss_group_id !< The ID of the group
   integer,private         :: iter_nuloss = 0 !< Iteration count,
                                              !< how often was already
                                              !< extended to the datasets?
   integer,private           :: chunk_counter_nuloss = 0 !< How full is the chunk already?
   integer,parameter,private :: chunk_size_nuloss=500    !< Chunk size
   real(r_kind),private,dimension(chunk_size_nuloss)  :: chunk_store_nuloss_t,    &
                                                         chunk_store_nuloss_temp, &
                                                         chunk_store_nuloss_dens, &
                                                         chunk_store_nuloss_rad , &
                                                         chunk_store_nuloss_nut , &
                                                         chunk_store_nuloss_bet , &
                                                         chunk_store_nuloss_the , &
                                                         chunk_store_nuloss_heat

   ! Time
   character(len=*), parameter,private      :: nuloss_dsetname_t = "time" !< Name of the dataset
   integer(HID_T),private                   :: nuloss_dset_id_t           !< Dataset identifier
   ! Temperature
   character(len=*), parameter,private      :: nuloss_dsetname_temp = "temp" !< Name of the dataset
   integer(HID_T),private                   :: nuloss_dset_id_temp           !< Dataset identifier
   ! Density
   character(len=*), parameter,private      :: nuloss_dsetname_dens = "dens" !< Name of the dataset
   integer(HID_T),private                   :: nuloss_dset_id_dens           !< Dataset identifier
   ! Radius
   character(len=*), parameter,private      :: nuloss_dsetname_rad = "rad"  !< Name of the dataset
   integer(HID_T),private                   :: nuloss_dset_id_rad           !< Dataset identifier
   ! nu loss total
   character(len=*), parameter,private      :: nuloss_dsetname_nut = "nu_total" !< Name of the dataset
   integer(HID_T),private                   :: nuloss_dset_id_nut               !< Dataset identifier
   ! nu loss beta
   character(len=*), parameter,private      :: nuloss_dsetname_bet = "nu_beta" !< Name of the dataset
   integer(HID_T),private                   :: nuloss_dset_id_bet              !< Dataset identifier
   ! nu loss thermal
   character(len=*), parameter,private      :: nuloss_dsetname_the = "nu_thermal" !< Name of the dataset
   integer(HID_T),private                   :: nuloss_dset_id_the                 !< Dataset identifier
   ! nu gain neutrino heating
   character(len=*), parameter,private      :: nuloss_dsetname_heat = "nu_heat" !< Name of the dataset
   integer(HID_T),private                   :: nuloss_dset_id_heat              !< Dataset identifier


   !----------- Finab variables --------------!
   integer(HID_T),private  :: finab_group_id !< The ID of the group



   !
   ! Public and private fields and methods of the module
   !
   public:: &
     init_hdf5_module, hdf5_module_finalize, extend_snaps, extend_mainout, &
     extend_track_nuclei, extend_av_timescales, extend_flow, write_finab, &
     extend_engen, extend_nu_loss_gain
   private:: &
     write_units, create_1d_dataset, create_constant_1d_int_arrays,&
     extend_1d_dataset, extend_hdf5_Y, write_snaps, write_track_data,&
     write_mainout_data, write_av_timescales, extend_1d_int_dataset,&
     create_1d_int_dataset, init_snaps, init_av_timescales, init_mainout,&
     init_track_nuclei, init_flow, create_constant_1d_arrays, init_engen,&
     write_engen, write_nu_loss_gain, init_nu_loss_gain


contains

   !> @brief Initialize the hdf5 module
   !!
   !! Create the file and the datasets. Furthermore, initialize arrays.
   !!
   !! @author Moritz Reichert
   !! @date 28.01.21
   subroutine init_hdf5_module()
      use parameter_class, only: h_snapshot_every,h_mainout_every,h_custom_snapshots,&
                                 h_track_nuclei_every,h_timescales_every,h_flow_every,&
                                 h_finab,h_engen_every,h_engen_detailed,h_nu_loss_every

      implicit none
      integer                     :: i_stat  !< status variable


      ! Check if a file has to get created
      if ((h_snapshot_every .gt. 0) .or. (h_custom_snapshots) .or. &
          (h_mainout_every .gt. 0) .or. (h_track_nuclei_every .gt. 0) .or. &
          (h_timescales_every .gt. 0) .or. (h_flow_every .gt. 0) .or.&
          (h_finab) .or. (h_engen_every .gt. 0) .or. &
          ((h_nu_loss_every .gt. 0))) then
          hdf5_output = .True.
      end if
      if (.not. hdf5_output) return

      !--- Create and initialize the file
      !
      ! Initialize FORTRAN predefined datatypes
      call h5open_f(i_stat)

      if (i_stat .ne. 0) call raise_exception("Could not initialize fortran predefined datatypes",&
                                              "init_hdf5_module",240003)

      ! Create a new file
      call h5fcreate_f(hdf5_filename, H5F_ACC_TRUNC_F, file_id, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Could not open HDF5 file."//NEW_LINE("A")//&
                              "Close the file if you have opened it.",&
                              "init_hdf5_module",240004)
      end if

      ! Snapshot things:
      if ((h_snapshot_every .gt. 0) .or. (h_custom_snapshots)) then
         call init_snaps
      end if

      ! Mainout things:
      if (h_mainout_every .gt. 0) then
         call init_mainout
      end if

      ! Open the track_nuclei file and write a header
      if (h_track_nuclei_every .gt. 0) then
         call init_track_nuclei
      end if

      ! Store average timescales
      if (h_timescales_every .gt. 0) then
         call init_av_timescales
      end if

      ! Initialize flow output
      if (h_flow_every .gt. 0) then
         call init_flow
      end if

      if ((h_nu_loss_every .gt. 0)) then
         call init_nu_loss_gain
      end if

      ! Initialize flow output
      if (h_engen_every .gt. 0) then
         call init_engen(h_engen_detailed)
      end if

   end subroutine init_hdf5_module



   !> Initialize the flow output
   subroutine init_flow
      implicit none
      integer                          :: i_stat  !< status variable

      ! Abundance dataset
      ! Create a group for the snapshots
      call h5gcreate_f(file_id, "/flows", flow_group_id, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create flows group.",&
                              "init_flow",240005)
      end if

   end subroutine init_flow


   !> Initialize neutrino loss/gain output
   subroutine init_nu_loss_gain
     implicit none
     ! Helper variables
     integer    :: i_stat  !< status variable

     ! Create group
     call h5gcreate_f(file_id, "/nu_loss_gain", nuloss_group_id, i_stat)
     if (i_stat .ne. 0) then
        call raise_exception("Unable to create nu_loss_gain group.",&
                             "init_nu_loss_gain",240005)
     end if

     ! Also make space for the time and all other outputs
     nuloss_dset_id_t    = create_1d_dataset(nuloss_dsetname_t,nuloss_group_id)
     nuloss_dset_id_temp = create_1d_dataset(nuloss_dsetname_temp,nuloss_group_id)
     nuloss_dset_id_dens = create_1d_dataset(nuloss_dsetname_dens,nuloss_group_id)
     nuloss_dset_id_rad  = create_1d_dataset(nuloss_dsetname_rad,nuloss_group_id)
     nuloss_dset_id_nut  = create_1d_dataset(nuloss_dsetname_nut,nuloss_group_id)
     nuloss_dset_id_the  = create_1d_dataset(nuloss_dsetname_the,nuloss_group_id)
     nuloss_dset_id_bet  = create_1d_dataset(nuloss_dsetname_bet,nuloss_group_id)
     nuloss_dset_id_heat = create_1d_dataset(nuloss_dsetname_heat,nuloss_group_id)

  end subroutine init_nu_loss_gain

   !>
   !! Write neutrino gain and loss
   !!
   !! These values are stored in a separate group (called snapshots) that is initialized
   !! in \ref init_hdf5_module.
   !!
   !! @author Moritz Reichert
   !! @date 12.04.23
  subroutine extend_nu_loss_gain(time,temp,dens,rad,the,heat,bet)
    use global_class, only: net_size
    implicit none
    real(r_kind),intent(in)                     :: time    !< Time [s]
    real(r_kind),intent(in)                     :: temp    !< Temperature [GK]
    real(r_kind),intent(in)                     :: dens    !< Density [g/cm^3]
    real(r_kind),intent(in)                     :: rad     !< Radius [km]
    real(r_kind),intent(in)                     :: the     !< Neutrino gain/loss thermal [MeV/baryon/s]
    real(r_kind),intent(in)                     :: heat    !< Neutrino gain/loss heating [MeV/baryon/s]
    real(r_kind),intent(in)                     :: bet     !< Neutrino gain/loss beta decay [MeV/baryon/s]
    real(r_kind)                                :: nut     !< Neutrino gain/loss total [MeV/baryon/s]
    ! Total energy
    nut = the + heat + bet
    ! Store the time
    chunk_counter_nuloss = chunk_counter_nuloss + 1
    chunk_store_nuloss_t(chunk_counter_nuloss)    = time
    chunk_store_nuloss_temp(chunk_counter_nuloss) = temp
    chunk_store_nuloss_dens(chunk_counter_nuloss) = dens
    chunk_store_nuloss_rad(chunk_counter_nuloss)  = rad
    chunk_store_nuloss_nut(chunk_counter_nuloss)  = nut
    chunk_store_nuloss_the(chunk_counter_nuloss)  = the
    chunk_store_nuloss_heat(chunk_counter_nuloss) = heat
    chunk_store_nuloss_bet(chunk_counter_nuloss)  = bet

    ! Write to the file
    if (chunk_counter_nuloss .eq. chunk_size_nuloss) then
       call write_nu_loss_gain()
    end if
 end subroutine extend_nu_loss_gain


   !>
   !! @brief Write the neutrino loss/gain into the hdf5 file
   !!
   !! These values are stored in a separate group (called nu_loss_gain) that is initialized
   !! in \ref init_hdf5_module.
   !!
   !! @author Moritz Reichert
   !! @date 12.04.23
 subroutine write_nu_loss_gain
    implicit none
    ! Helper variables

    if (chunk_counter_nuloss .lt. 1) return
    ! Save the data
    call extend_1d_dataset(nuloss_dset_id_t,chunk_store_nuloss_t(1:chunk_counter_nuloss),iter_nuloss,chunk_counter_nuloss)
    call extend_1d_dataset(nuloss_dset_id_temp,chunk_store_nuloss_temp(1:chunk_counter_nuloss),iter_nuloss,chunk_counter_nuloss)
    call extend_1d_dataset(nuloss_dset_id_dens,chunk_store_nuloss_dens(1:chunk_counter_nuloss),iter_nuloss,chunk_counter_nuloss)
    call extend_1d_dataset(nuloss_dset_id_rad,chunk_store_nuloss_rad(1:chunk_counter_nuloss),iter_nuloss,chunk_counter_nuloss)
    call extend_1d_dataset(nuloss_dset_id_nut,chunk_store_nuloss_nut(1:chunk_counter_nuloss),iter_nuloss,chunk_counter_nuloss)
    call extend_1d_dataset(nuloss_dset_id_the,chunk_store_nuloss_the(1:chunk_counter_nuloss),iter_nuloss,chunk_counter_nuloss)
    call extend_1d_dataset(nuloss_dset_id_heat,chunk_store_nuloss_heat(1:chunk_counter_nuloss),iter_nuloss,chunk_counter_nuloss)
    call extend_1d_dataset(nuloss_dset_id_bet,chunk_store_nuloss_bet(1:chunk_counter_nuloss),iter_nuloss,chunk_counter_nuloss)

    iter_nuloss = iter_nuloss+chunk_counter_nuloss
    chunk_counter_nuloss = 0
 end subroutine write_nu_loss_gain





   !> Initialize energy generation output
   subroutine init_engen(include_decay)
     use global_class, only: net_size,isotope
      implicit none
      logical,intent(in) :: include_decay !< Flag that determines if decay energy per parent is written
      ! Helper variables
      integer    :: i_stat  !< status variable
      integer    :: i       !< Loop variable
      logical    :: avail
      integer    :: filter_info,filter_info_both
      integer,dimension(net_size) :: helper  !< Helper array to store A,Z,N


      ! Create group
      call h5gcreate_f(file_id, "/energy", en_group_id, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create energies group.",&
                              "init_engen",240005)
      end if
      ! Store dataset Ids
      do i=1, 10
         en_reac_ids(i) = create_1d_dataset(trim(adjustl(en_names(i))),en_group_id)
      end do

      ! Also make space for the time
      en_dset_id_t = create_1d_dataset(en_dsetname_t,en_group_id)

      if (include_decay .eqv. .true.) then
         ! Set the initial dimensions
         dims_Y(1:2) = (/net_size,1/)
         ! Create the data space with unlimited dimensions.
         maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
         call h5screate_simple_f(2, dims_Y, dataspace, i_stat, maxdims)
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create necessary dataspace.",&
                                 "init_engen",240006)
         end if


         ! Check if gzip compression is available and can be used for both
         ! compression and decompression.
         ! The compression can greatly reduce the used space as
         ! The energy generation has many zeros included. In a small
         ! test, the filesize got reduced by a factor of 10.
         ! However, it consumes some time to compress whenever data is written...
         CALL h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, avail, i_stat)

         if (.not. avail) then
           call raise_exception("gzip filter not available.","init_engen",240008)
         end if
         call h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, filter_info, i_stat)

         filter_info_both=ior(H5Z_FILTER_ENCODE_ENABLED_F,H5Z_FILTER_DECODE_ENABLED_F)
         if (filter_info .ne. filter_info_both) then
           call raise_exception("gzip filter not available for encoding and decoding.",&
                                "init_engen",240009)
         endif

         ! Modify dataset properties
         call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, i_stat)
         if (i_stat .ne. 0) then
            call raise_exception("Unable to modify dataset properties.",&
                                 "init_engen",240010)
         end if
         ! Apply the gzip filter
         call h5pset_deflate_f(crp_list, 9, i_stat)
         if (i_stat .ne. 0) then
           call raise_exception("Unable to deflate and apply gzip filter.",&
                                "init_engen",240011)
         end if
         ! Enable chunking
         call h5pset_chunk_f(crp_list, 2, dims_Y, i_stat)
         if (i_stat .ne. 0) then
            call raise_exception("Unable to enable chunking.",&
                                 "init_engen",240012)
         end if

         ! Create a dataset using cparms creation property
         call h5dcreate_f(en_group_id, en_det_bet_name, H5T_NATIVE_DOUBLE, dataspace, &
             en_det_bet_id, i_stat, crp_list )
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create dataspace (en_det_bet_id).",&
                                 "init_engen",240006)
         end if
         call h5dcreate_f(en_group_id, en_det_ag_name, H5T_NATIVE_DOUBLE, dataspace, &
             en_det_ag_id, i_stat, crp_list )
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create dataspace (en_det_ag_id).",&
                                 "init_engen",240006)
         end if
         call h5dcreate_f(en_group_id, en_det_pg_name, H5T_NATIVE_DOUBLE, dataspace, &
             en_det_pg_id, i_stat, crp_list )
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create dataspace (en_det_pg_id).",&
                                 "init_engen",240006)
         end if
         call h5dcreate_f(en_group_id, en_det_ng_name, H5T_NATIVE_DOUBLE, dataspace, &
             en_det_ng_id, i_stat, crp_list )
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create dataspace (en_det_ng_id).",&
                                 "init_engen",240006)
         end if
         call h5dcreate_f(en_group_id, en_det_pn_name, H5T_NATIVE_DOUBLE, dataspace, &
             en_det_pn_id, i_stat, crp_list )
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create dataspace (en_det_pn_id).",&
                                 "init_engen",240006)
         end if
         call h5dcreate_f(en_group_id, en_det_ap_name, H5T_NATIVE_DOUBLE, dataspace, &
             en_det_ap_id, i_stat, crp_list )
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create dataspace (en_det_ap_id).",&
                                 "init_engen",240006)
         end if
         call h5dcreate_f(en_group_id, en_det_an_name, H5T_NATIVE_DOUBLE, dataspace, &
             en_det_an_id, i_stat, crp_list )
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create dataspace (en_det_an_id).",&
                                 "init_engen",240006)
         end if
         call h5dcreate_f(en_group_id, en_det_other_name, H5T_NATIVE_DOUBLE, dataspace, &
             en_det_other_id, i_stat, crp_list )
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create dataspace (en_det_other_id).",&
                                 "init_engen",240006)
         end if
         call h5dcreate_f(en_group_id, en_det_fiss_name, H5T_NATIVE_DOUBLE, dataspace, &
             en_det_fiss_id, i_stat, crp_list )
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create dataspace (en_det_fiss_id).",&
                                 "init_engen",240006)
         end if
         call h5dcreate_f(en_group_id, en_det_tot_name, H5T_NATIVE_DOUBLE, dataspace, &
             en_det_tot_id, i_stat, crp_list )
         if (i_stat .ne. 0) then
            call raise_exception("Unable to create dataspace (en_det_tot_id).",&
                                 "init_engen",240006)
         end if
         call h5sclose_f(dataspace, i_stat) ! close the dataspace again
         if (i_stat .ne. 0) then
            call raise_exception("Unable to close the dataspace.",&
                                 "init_engen",240007)
         end if

         !--- Write constant information into the hdf5
         ! Write A, Z, N into the hdf5
         ! Necessary to identify nuclei in detailed output
         do i=1,net_size
            helper(i) = isotope(i)%mass
         end do
         call create_constant_1d_int_arrays(net_size,helper,"A",en_group_id)
         ! Write Z into hdf5
         do i=1,net_size
            helper(i) = isotope(i)%p_nr
         end do
         call create_constant_1d_int_arrays(net_size,helper,"Z",en_group_id)
         ! Write N into hdf5
         do i=1,net_size
            helper(i) = isotope(i)%n_nr
         end do
         call create_constant_1d_int_arrays(net_size,helper,"N",en_group_id)

         ! Allocate the array to store the chunks in
         allocate(chunk_store_en_bet(net_size,chunk_size_en),chunk_store_en_ag(net_size,chunk_size_en), &
                   chunk_store_en_pg(net_size,chunk_size_en),chunk_store_en_ng(net_size,chunk_size_en), &
                   chunk_store_en_pn(net_size,chunk_size_en),chunk_store_en_ap(net_size,chunk_size_en), &
                   chunk_store_en_an(net_size,chunk_size_en),chunk_store_en_other(net_size,chunk_size_en), &
                   chunk_store_en_fiss(net_size,chunk_size_en),chunk_store_en_tot(net_size,chunk_size_en),stat=i_stat)

         if (i_stat /= 0) call raise_exception('Allocation failed.',"init_engen",&
                                               240001)
      end if

   end subroutine init_engen


   !>
   !! @brief Subroutine to write the energy generation per reaction type
   !!
   !! @author Moritz Reichert
   !! @date 22.03.21
   subroutine extend_engen(time,engen_tot, heat,engen_ng_gn, engen_pg_gp, &
                           engen_ag_ga, engen_np_pn, engen_an_na, &
                           engen_ap_pa, engen_beta, engen_fiss, &
                           det_bet_engen, det_ag_engen, det_pg_engen, &
                           det_ng_engen, det_pn_engen, det_ap_engen, &
                           det_an_engen, det_other_engen, det_fiss_engen, &
                           det_tot_engen, &
                           include_decay)
      use global_class, only: net_size
      implicit none
      real(r_kind),intent(in)  :: time                                       !< Time [s]
      real(r_kind),intent(in)  :: engen_tot, heat, engen_ng_gn, engen_pg_gp  !< Energy generation [erg/g/s]
      real(r_kind),intent(in)  :: engen_ag_ga, engen_np_pn, engen_an_na      !< Energy generation [erg/g/s]
      real(r_kind),intent(in)  :: engen_ap_pa, engen_beta, engen_fiss        !< Energy generation [erg/g/s]
      real(r_kind),dimension(net_size),intent(in) :: det_bet_engen
      real(r_kind),dimension(net_size),intent(in) :: det_ag_engen
      real(r_kind),dimension(net_size),intent(in) :: det_pg_engen
      real(r_kind),dimension(net_size),intent(in) :: det_ng_engen
      real(r_kind),dimension(net_size),intent(in) :: det_pn_engen
      real(r_kind),dimension(net_size),intent(in) :: det_ap_engen
      real(r_kind),dimension(net_size),intent(in) :: det_an_engen
      real(r_kind),dimension(net_size),intent(in) :: det_other_engen
      real(r_kind),dimension(net_size),intent(in) :: det_fiss_engen      !< Energy generation [erg/g/s]
      real(r_kind),dimension(net_size),intent(in) :: det_tot_engen       !< Energy generation [erg/g/s]
      logical,intent(in)                          :: include_decay           !< Write engen_decay_parents?
      !
      integer                   :: i  !< Loop variable

      ! Count how large the space has to be
      chunk_counter_en = chunk_counter_en + 1
      !< Save the data temporary in an array
      chunk_store_en(chunk_counter_en,:) = (/time,engen_tot, heat,engen_ng_gn, engen_pg_gp, &
                                             engen_ag_ga, engen_np_pn, engen_an_na, &
                                             engen_ap_pa, engen_beta, engen_fiss/)

      if (include_decay .eqv. .true.) then
         ! Extend the decay energy generation splitted in parent nuclei
         chunk_store_en_bet(:,chunk_counter_en)  = det_bet_engen
         chunk_store_en_ag(:,chunk_counter_en)   = det_ag_engen
         chunk_store_en_pg(:,chunk_counter_en)   = det_pg_engen
         chunk_store_en_ng(:,chunk_counter_en)   = det_ng_engen
         chunk_store_en_pn(:,chunk_counter_en)   = det_pn_engen
         chunk_store_en_ap(:,chunk_counter_en)   = det_ap_engen
         chunk_store_en_an(:,chunk_counter_en)   = det_an_engen
         chunk_store_en_other(:,chunk_counter_en)= det_other_engen
         chunk_store_en_fiss(:,chunk_counter_en) = det_fiss_engen
         chunk_store_en_tot(:,chunk_counter_en)  = det_tot_engen
      end if

      ! Write the chunk to the file
      if (chunk_counter_en .eq. chunk_size_en) then
         call write_engen(include_decay)
      end if

   end subroutine extend_engen

   !>
   !! @brief Write energy generations and time into the hdf5 file
   !!
   !! These values are stored in a separate group (called energy) that is initialized
   !! in \ref init_hdf5_module.
   !!
   !! @author Moritz Reichert
   !! @date 22.03.21
   subroutine write_engen(include_decay)
      use global_class, only: net_size
      implicit none
      logical, intent(in) :: include_decay !< Write engen_decay?
      ! Helper variables
      integer :: i !< Loop variable

      if (chunk_counter_en .lt. 1) return
      ! Save the time
      call extend_1d_dataset(en_dset_id_t,chunk_store_en(1:chunk_counter_en,1),iter_en,chunk_counter_en)
      ! Save timescales
      do i=1, 10
         call extend_1d_dataset(en_reac_ids(i),chunk_store_en(1:chunk_counter_en,i+1),iter_en,chunk_counter_en)
      end do

      if (include_decay .eqv. .true.) then
         ! Decay energy splitted into parent nuclei
         call extend_hdf5_Y(en_det_bet_id  ,chunk_store_en_bet  ,net_size,iter_en,chunk_counter_en)
         call extend_hdf5_Y(en_det_ag_id   ,chunk_store_en_ag   ,net_size,iter_en,chunk_counter_en)
         call extend_hdf5_Y(en_det_pg_id   ,chunk_store_en_pg   ,net_size,iter_en,chunk_counter_en)
         call extend_hdf5_Y(en_det_ng_id   ,chunk_store_en_ng   ,net_size,iter_en,chunk_counter_en)
         call extend_hdf5_Y(en_det_pn_id   ,chunk_store_en_pn   ,net_size,iter_en,chunk_counter_en)
         call extend_hdf5_Y(en_det_ap_id   ,chunk_store_en_ap   ,net_size,iter_en,chunk_counter_en)
         call extend_hdf5_Y(en_det_an_id   ,chunk_store_en_an   ,net_size,iter_en,chunk_counter_en)
         call extend_hdf5_Y(en_det_other_id,chunk_store_en_other,net_size,iter_en,chunk_counter_en)
         call extend_hdf5_Y(en_det_fiss_id ,chunk_store_en_fiss ,net_size,iter_en,chunk_counter_en)
         call extend_hdf5_Y(en_det_tot_id  ,chunk_store_en_tot  ,net_size,iter_en,chunk_counter_en)
      end if

      iter_en = iter_en+chunk_counter_en
      chunk_counter_en = 0
   end subroutine write_engen




   !> Extend the flow in the file
   !! With every call a new group is created and the flows are written in this
   !! group. The group is named after the calls (the third call is then located
   !! in flows/3)
   subroutine extend_flow(time,temp,dens,flow,n_flows,Y)
      use global_class, only: flow_vector,isotope,net_size
      implicit none
      real(r_kind),intent(in)   :: time                          !< Value of the time
      real(r_kind),intent(in)   :: temp                          !< Value of the time
      real(r_kind),intent(in)   :: dens                          !< Value of the time
      integer, intent(in)       :: n_flows                       !< Maximum number of flows
      type(flow_vector),dimension(n_flows), intent(in)  :: flow  !< Array with already calculated flows
      real(r_kind),dimension(net_size), intent(in)      :: Y     !< Array with abundances
      integer                                :: i                !< Loop variable
      integer                                :: count            !< Amount of non zero flows
      integer                                :: i_stat           !< IO status variable
      integer                                :: alloc_stat       !< Allocation status variable
      real(r_kind)                           :: flow_diff        !< Value of the flow
      integer,dimension(:),allocatable       :: n_in_nr,p_in_nr  !< Amount of ingoing neutrons and protons
      integer,dimension(:),allocatable       :: n_out_nr,p_out_nr!< Amount of outgoing neutrons and protons
      real(r_kind),dimension(:),allocatable  :: flow_in,flow_out !< In/Out-going flow
      real(r_kind),dimension(:),allocatable  :: Y_in,Y_out       !< In/Out-going abundances
      integer(HID_T)                         :: tmp_group_id     !< Temporary group id of the flows
      real(r_kind),dimension(1)              :: helper_time      !< Helper array to store the time
      real(r_kind),dimension(1)              :: helper_temp      !< Helper array to store the temp
      real(r_kind),dimension(1)              :: helper_dens      !< Helper array to store the dens

      iter_flow = iter_flow + 1 ! Count the amount of calls

      ! Initialize count variable
      count = 0

      ! Loop over the flows and count non zero ones
      do i=1,n_flows
         flow_diff = flow(i)%fwd - flow(i)%bwd
         ! Count non zero flows
         if (flow_diff .ne. 0) count = count + 1
      end do

      ! Allocate the arrays
      allocate(n_out_nr(count),stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Allocation of 'n_out_nr' failed!","extend_flow",240001)
      allocate(p_out_nr(count),stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Allocation of 'p_out_nr' failed!","extend_flow",240001)
      allocate(flow_out(count),stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Allocation of 'flow_out' failed!","extend_flow",240001)
      allocate( n_in_nr(count),stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Allocation of 'n_in_nr' failed!","extend_flow",240001)
      allocate( p_in_nr(count),stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Allocation of 'p_in_nr' failed!","extend_flow",240001)
      allocate( flow_in(count),stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Allocation of 'flow_in' failed!","extend_flow",240001)
      allocate( Y_in(count)   ,stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Allocation of 'Y_in' failed!","extend_flow",240001)
      allocate( Y_out(count)  ,stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Allocation of 'Y_out' failed!","extend_flow",240001)

      ! Loop over the flows and fill the arrays with non zero flows
      count = 0
      do i=1,n_flows
         flow_diff = flow(i)%fwd - flow(i)%bwd
         ! Count non zero flows
         if (flow_diff .eq. 0) cycle
         count = count + 1
         ! Store ingoing and outgoing quantities
         n_in_nr(count)  = isotope(flow(i)%iin)%n_nr
         p_in_nr(count)  = isotope(flow(i)%iin)%p_nr
         n_out_nr(count) = isotope(flow(i)%iout)%n_nr
         p_out_nr(count) = isotope(flow(i)%iout)%p_nr
         flow_out(count) = flow(i)%bwd
         flow_in(count)  = flow(i)%fwd
         Y_in(count)     = Y(flow(i)%iin)
         Y_out(count)    = Y(flow(i)%iout)
      end do
      ! Store the time in form of an array with length 1
      helper_time(1) = time
      helper_temp(1) = temp
      helper_dens(1) = dens

      ! Abundance dataset
      ! Create a group for the snapshots
      call h5gcreate_f(file_id, "/flows/"//int_to_str(iter_flow), tmp_group_id, i_stat)
      call create_constant_1d_int_arrays(count,n_in_nr ,"n_in",tmp_group_id)
      call create_constant_1d_int_arrays(count,p_in_nr ,"p_in",tmp_group_id)
      call create_constant_1d_int_arrays(count,n_out_nr,"n_out",tmp_group_id)
      call create_constant_1d_int_arrays(count,p_out_nr,"p_out",tmp_group_id)
      call create_constant_1d_arrays(count,flow_in-flow_out,"flow",tmp_group_id)
      call create_constant_1d_arrays(1,helper_time,"time",tmp_group_id)
      call create_constant_1d_arrays(1,helper_temp,"temp",tmp_group_id)
      call create_constant_1d_arrays(1,helper_dens,"dens",tmp_group_id)
      call create_constant_1d_arrays(count,Y_in,"y_in",tmp_group_id)
      call create_constant_1d_arrays(count,Y_out,"y_out",tmp_group_id)
      call h5gclose_f(tmp_group_id  , i_stat)


      ! Deallocate the arrays
      deallocate(n_out_nr,stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Deallocation of 'n_out_nr' failed!","extend_flow",240002)
      deallocate(p_out_nr,stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Deallocation of 'p_out_nr' failed!","extend_flow",240002)
      deallocate(flow_out,stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Deallocation of 'flow_out' failed!","extend_flow",240002)
      deallocate( n_in_nr,stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Deallocation of 'n_in_nr' failed!","extend_flow",240002)
      deallocate( p_in_nr,stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Deallocation of 'p_in_nr' failed!","extend_flow",240002)
      deallocate( flow_in,stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Deallocation of 'flow_in' failed!","extend_flow",240002)
      deallocate(Y_in ,stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Deallocation of 'Y_in' failed!","extend_flow",240002)
      deallocate(Y_out,stat=alloc_stat)
      if (alloc_stat .ne. 0) call raise_exception("Deallocation of 'Y_out' failed!","extend_flow",240002)

   end subroutine extend_flow


   !> Initialize the snapshot hdf5 group.
   !!
   !! This subroutine creates a snapshot group and writes mass number
   !! proton number and neutron number into it. Furthermore, it creates a dataspace
   !! for the abundances and the time. In addition a array that stores
   !! the abundances is initialized. This array is used to create larger chunks
   !! and reduce the IO when writing to the file on the cost of memory usage.
   subroutine init_snaps
      use global_class, only: net_size,isotope
      implicit none
      integer                     :: i_stat  !< status variable
      integer                     :: i       !< Loop variable
      integer,dimension(net_size) :: helper  !< Helper array to store A,Z,N
      logical                     :: avail   !< For gzip compression
      integer                     :: filter_info,filter_info_both
      ! Abundance dataset
      ! Create a group for the snapshots
      call h5gcreate_f(file_id, "/snapshots", snaps_group_id, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create snapshot group.",&
                              "init_snaps",240005)
      end if

      ! Set the initial dimensions
      dims_Y(1:2) = (/net_size,1/)
      ! Create the data space with unlimited dimensions.
      maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
      call h5screate_simple_f(2, dims_Y, dataspace, i_stat, maxdims)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create necessary dataspace.",&
                              "init_snaps",240006)
      end if

      ! Check if gzip compression is available and can be used for both
      ! compression and decompression.
      call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, avail, i_stat)

      if (.not. avail) then
        call raise_exception("gzip filter not available.","init_snaps",240008)
      end if
      call h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, filter_info, i_stat)

      filter_info_both=ior(H5Z_FILTER_ENCODE_ENABLED_F,H5Z_FILTER_DECODE_ENABLED_F)
      if (filter_info .ne. filter_info_both) then
        call raise_exception("gzip filter not available for encoding and decoding.","init_snaps",240009)
      endif

      ! Modify dataset properties
      call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to modify dataset properties.",&
                              "init_snaps",240010)
      end if

      ! Apply the gzip filter
      call h5pset_deflate_f(crp_list, 9, i_stat)
      if (i_stat .ne. 0) then
        call raise_exception("Unable to deflate and apply gzip filter.",&
                             "init_snaps",240011)
      end if

      call h5pset_chunk_f(crp_list, 2, dims_Y, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to enable chunking.",&
                              "init_snaps",240012)
      end if

      ! Create a dataset using cparms creation property
      call h5dcreate_f(snaps_group_id, dsetname_Y, H5T_NATIVE_DOUBLE, dataspace, &
          dset_id_Y, i_stat, crp_list )
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create dataspace.",&
                              "init_snaps",240006)
      end if
      call h5sclose_f(dataspace, i_stat) ! close the dataspace again
      if (i_stat .ne. 0) then
         call raise_exception("Unable to close the dataspace.",&
                              "init_snaps",240007)
      end if

      !--- Initialize 1D datasets
      ! Time dataset
      snaps_dset_id_t = create_1d_dataset(snaps_dsetname_t,snaps_group_id)

      !--- Write constant information into the hdf5
      ! Write A into the hdf5
      do i=1,net_size
         helper(i) = isotope(i)%mass
      end do
      call create_constant_1d_int_arrays(net_size,helper,"A",snaps_group_id)
      ! Write Z into hdf5
      do i=1,net_size
         helper(i) = isotope(i)%p_nr
      end do
      call create_constant_1d_int_arrays(net_size,helper,"Z",snaps_group_id)
      ! Write N into hdf5
      do i=1,net_size
         helper(i) = isotope(i)%n_nr
      end do
      call create_constant_1d_int_arrays(net_size,helper,"N",snaps_group_id)

      ! Allocate the array to store the chunks in
      allocate(chunk_store_snaps_Y(net_size,chunk_size_snaps))

   end subroutine init_snaps


   !> Initialize the mainout hdf5 group.
   !!
   !! This subroutine creates a mainout group and creates the dataspaces for
   !! all included quantities.
   subroutine init_mainout
      implicit none
      integer                     :: i_stat  !< status variable

      ! Create a group for the mainout
      call h5gcreate_f(file_id, "/mainout", mainout_group_id, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create mainout group.",&
                              "init_mainout",240005)
      end if
      !--- Initialize 1D datasets
      ! Iteration (integer)dataset
      mainout_dset_id_iter  = create_1d_int_dataset(mainout_dsetname_iter,mainout_group_id)
      ! Time dataset
      mainout_dset_id_t     = create_1d_dataset(mainout_dsetname_t,mainout_group_id)
      ! Temp. dataset
      mainout_dset_id_temp  = create_1d_dataset(mainout_dsetname_temp,mainout_group_id)
      ! Density dataset
      mainout_dset_id_dens  = create_1d_dataset(mainout_dsetname_dens,mainout_group_id)
      ! Entropy dataset
      mainout_dset_id_entr  = create_1d_dataset(mainout_dsetname_entr,mainout_group_id)
      ! Radius dataset
      mainout_dset_id_rad   = create_1d_dataset(mainout_dsetname_rad,mainout_group_id)
      ! Electron fraction dataset
      mainout_dset_id_ye    = create_1d_dataset(mainout_dsetname_ye,mainout_group_id)
      ! neutron abundance dataset
      mainout_dset_id_yn    = create_1d_dataset(mainout_dsetname_yn,mainout_group_id)
      ! proton abundance dataset
      mainout_dset_id_yp    = create_1d_dataset(mainout_dsetname_yp,mainout_group_id)
      ! alpha abundance dataset
      mainout_dset_id_ya    = create_1d_dataset(mainout_dsetname_ya,mainout_group_id)
      ! ylight abundance dataset
      mainout_dset_id_ylight= create_1d_dataset(mainout_dsetname_ylight,mainout_group_id)
      ! yheavy abundance dataset
      mainout_dset_id_yheavy= create_1d_dataset(mainout_dsetname_yheavy,mainout_group_id)
      ! abar abundance dataset
      mainout_dset_id_abar  = create_1d_dataset(mainout_dsetname_abar,mainout_group_id)
      ! zbar abundance dataset
      mainout_dset_id_zbar  = create_1d_dataset(mainout_dsetname_zbar,mainout_group_id)
   end subroutine init_mainout


   !> Initialize the tracked_nuclei hdf5 group.
   !!
   !! This subroutine creates a tracked_nuclei group and creates the dataspaces for
   !! all included abundances and the time. Furthermore, it allocates a
   !! array to store the abundances over several iterations.
   subroutine init_track_nuclei
      use global_class, only: track_nuclei_nr,track_nuclei_indices,&
                              net_names,isotope
      implicit none
      integer                            :: i_stat  !< status variable
      integer                            :: i       !< Loop variable
      integer,dimension(track_nuclei_nr) :: helper  !< Helper array to store A,Z,N


      ! Create a group for the tracked nuclei
      call h5gcreate_f(file_id, "/tracked_nuclei", track_group_id, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create tracked_nuclei group.",&
                              "init_track_nuclei",240005)
      end if


      ! Set the initial dimensions
      dims_track(1:2) = (/track_nuclei_nr,1/)
      ! Create the data space with unlimited dimensions.
      maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
      call h5screate_simple_f(2, dims_track, dataspace, i_stat, maxdims)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create necessary dataspace.",&
                              "init_snaps",240006)
      end if

      ! Modify dataset properties
      call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to modify dataset properties.",&
                              "init_snaps",240010)
      end if

      call h5pset_chunk_f(crp_list, 2, dims_track, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to enable chunking.",&
                              "init_snaps",240012)
      end if

      ! Create a dataset using cparms creation property
      call h5dcreate_f(track_group_id, tracked_dsetname_Y, H5T_NATIVE_DOUBLE, dataspace, &
          tracked_abu_id, i_stat, crp_list )
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create dataspace.",&
                              "init_snaps",240006)
      end if
      call h5sclose_f(dataspace, i_stat) ! close the dataspace again
      if (i_stat .ne. 0) then
         call raise_exception("Unable to close the dataspace.",&
                              "init_snaps",240007)
      end if

      ! Store the names
      allocate(tracked_names(track_nuclei_nr),stat=i_stat)   ! Allocate the name array
      if (i_stat .ne. 0) then
         call raise_exception('Allocation of "tracked_names" failed.',&
                              "init_track_nuclei",240001)
      end if

      ! Convert indices to names
      do i=1, track_nuclei_nr
        tracked_names(i) = net_names(track_nuclei_indices(i))
      end do

      !--- Write constant information into the hdf5
      ! Write A into the hdf5
      do i=1,track_nuclei_nr
         helper(i) = isotope(track_nuclei_indices(i))%mass
      end do
      call create_constant_1d_int_arrays(track_nuclei_nr,helper,"A",track_group_id)
      ! Write Z into hdf5
      do i=1,track_nuclei_nr
         helper(i) = isotope(track_nuclei_indices(i))%p_nr
      end do
      call create_constant_1d_int_arrays(track_nuclei_nr,helper,"Z",track_group_id)
      ! Write N into hdf5
      do i=1,track_nuclei_nr
         helper(i) = isotope(track_nuclei_indices(i))%n_nr
      end do
      call create_constant_1d_int_arrays(track_nuclei_nr,helper,"N",track_group_id)

      ! Also make space for the time
      tracked_dset_id_t = create_1d_dataset(tracked_dsetname_t,track_group_id)
      allocate(chunk_store_track(track_nuclei_nr+1,chunk_size_track),stat=i_stat) ! +1 for the time
      if (i_stat .ne. 0) then
         call raise_exception('Allocation of "chunk_store_track" failed.',&
                              "init_track_nuclei",240001)
      end if
   end subroutine init_track_nuclei


   !> Initialize the timescales hdf5 group.
   !!
   !! This subroutine creates a timescales group and creates the dataspaces for
   !! all included timescales and the time.
   subroutine init_av_timescales
      implicit none
      integer                     :: i_stat  !< status variable
      integer                     :: i       !< Loop variable

      ! Create a group for the tracked nuclei
      call h5gcreate_f(file_id, "/timescales", ts_group_id, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create timescales group.",&
                              "init_av_timescales",240005)
      end if
      ! Store dataset Ids
      do i=1, 17
         ts_reac_ids(i) = create_1d_dataset(trim(adjustl(ts_names(i))),ts_group_id)
      end do
      ! Also make space for the time
      ts_dset_id_t = create_1d_dataset(ts_dsetname_t,ts_group_id)

   end subroutine init_av_timescales



   !> @brief Write the units of the individual entries into the hdf5 file
   !!
   !! For example: \n
   !! Temperature - GK\n
   !! Entropy     - kB/nuc
   !!
   !! @note This subroutine is never used
   !!
   !! @author Moritz Reichert
   !! @date 28.01.21
   subroutine write_units()
      implicit none
      integer,parameter                     :: length=10 !< length of entries
      character(len=20),dimension(length,2) :: quantity  !< Storage of the text to write
      integer(HID_T)                        :: dset_id   !< ID of the dataset
      integer(HSIZE_T), dimension(2)        :: dims_static = (/length,2/)  !< Dimensions (length,2)
      integer                               :: i_stat    !< status variable
      INTEGER(HID_T)  :: strtype
      integer(SIZE_T) :: typesize
      ! Write the units of every quantity into the hdf5 file
      quantity(1 ,1) = "Abundances"        ; quantity(1 ,2) = "-"
      quantity(2 ,1) = "Atomic number"     ; quantity(2 ,2) = "-"
      quantity(3 ,1) = "Density"           ; quantity(3 ,2) = "g/ccm"
      quantity(4 ,1) = "Electron fraction" ; quantity(4 ,2) = "-"
      quantity(5 ,1) = "Entropy"           ; quantity(5 ,2) = "kB/nuc"
      quantity(6 ,1) = "Mass number"       ; quantity(6 ,2) = "-"
      quantity(7 ,1) = "Neutron number"    ; quantity(7 ,2) = "-"
      quantity(8 ,1) = "Radius"            ; quantity(8 ,2) = "Km"
      quantity(9 ,1) = "Temperature"       ; quantity(9 ,2) = "GK"
      quantity(10,1) = "Time"              ; quantity(10,2) = "s"


      ! Generate a datatype that is a character with length of 20
      call h5tcopy_f(H5T_NATIVE_CHARACTER, strtype, i_stat)
      typesize = 20
      call h5tset_size_f(strtype, typesize, i_stat)

      ! Create an appropriate dataspace
      call h5screate_simple_f(2, dims_static, dataspace, i_stat)

      ! Create a dataset using cparms creation property
      call h5dcreate_f(file_id, dsetname_u, strtype, dataspace, &
                       dset_id, i_stat )
      call h5sclose_f(dataspace, i_stat)

      ! Write initial value to the dataset
      call h5dwrite_f(dset_id, strtype, quantity, dims_static, i_stat)

      ! Close the dataset
      call h5dclose_f(dset_id, i_stat)

   end subroutine write_units



   !> Write final abundances into hdf5 file
   !!
   !! This routine writes the final abundances of all nuclei (group:finab/finab),
   !! the final abundances summed over mass number (group:finab/finabsum),
   !! and the final abundances summed over equal proton numbers (group:finab/finabelem).
   !! These entries are similar to the ascii version the output (finabelem.dat,..)
   !!
   !! @author M. Reichert
   !! @date  05.03.21
   subroutine write_finab(Y)
     use global_class, only: isotope, net_size
     implicit none

     real(r_kind),dimension(net_size),intent(in)  :: Y                  !< Abundances
     real(r_kind),dimension(:),allocatable        :: Y_tmpA,X_tmpA      !< Temporary helper variables
     real(r_kind),dimension(:),allocatable        :: Y_tmp,X_tmp        !< Temporary helper variables
     integer,dimension(:),allocatable             :: A_tmp,Z_tmp        !< Temporary helper variables
     integer                                      :: i                  !< Loop variable
     integer                                      :: count              !< Helper variable
     integer                                      :: i_stat             !< Status variable
     integer                                      :: mval               !< helper variable

     ! Create a group for the finab
     call h5gcreate_f(file_id, "/finab/", finab_group_id, i_stat)
     if (i_stat .ne. 0) then
        call raise_exception("Unable to create finab group.",&
                             "write_finab",240005)
     end if
     call h5gclose_f(finab_group_id    , i_stat)
     call h5gcreate_f(file_id, "/finab/finab", finab_group_id, i_stat)
     if (i_stat .ne. 0) then
        call raise_exception("Unable to create finab group.",&
                             "write_finab",240005)
     end if


     ! Write analog entry to the ascii file finab.dat
     ! Count the amount of abundances above 1d-25 for the output
     count = 0
     do i=1,net_size
        if (Y(i) .gt. 1d-25) then
           count = count+1
        end if
     end do

     ! Allocate arrays
     allocate(Y_tmp(count))
     allocate(A_tmp(count))
     allocate(Z_tmp(count))

     ! Fill the temporary arrays
     count = 0
     do i=1,net_size
       if (Y(i) .gt. 1d-25) then
          count = count+1
          Y_tmp(count) = Y(i)
          A_tmp(count) = isotope(i)%mass
          Z_tmp(count) = isotope(i)%p_nr
       end if
     end do

     ! Write to the hdf5 file
     call create_constant_1d_int_arrays(count,A_tmp,"A",finab_group_id)
     call create_constant_1d_int_arrays(count,Z_tmp,"Z",finab_group_id)
     call create_constant_1d_arrays(count,Y_tmp,"Y",finab_group_id)
     call create_constant_1d_arrays(count,Y_tmp*A_tmp,"X",finab_group_id)

     ! Deallocate arrays
     deallocate(Y_tmp)
     deallocate(A_tmp)
     deallocate(Z_tmp)

     ! Close the group
     call h5gclose_f(finab_group_id    , i_stat)

     ! Create group for finabsum
     call h5gcreate_f(file_id, "/finab/finabsum", finab_group_id, i_stat)
     if (i_stat .ne. 0) then
      call raise_exception("Unable to create finabsum group.",&
                           "write_finab",240005)
     end if

     mval = maxval(isotope%mass)
     allocate(Y_tmpA(mval))

     ! Sum up abundances
     Y_tmpA(:) = 0
     do i=1,net_size
        Y_tmpA(isotope(i)%mass)=Y_tmpA(isotope(i)%mass)+Y(i)
     end do

     count = 0
     do i=1,mval
        if (Y_tmpA(i) .gt. 1d-20) then
           count = count+1
        end if
     end do

     ! Allocate arrays once more
     allocate(Y_tmp(count))
     allocate(A_tmp(count))

     ! Prepare helper arrays for the output
     count = 0
     do i=1,mval
       if (Y_tmpA(i) .gt. 1d-20) then
          count = count+1
          A_tmp(count) = i
          Y_tmp(count) = Y_tmpA(i)
       end if
     end do

     ! Write the arrays
     call create_constant_1d_int_arrays(count,A_tmp,"A",finab_group_id)
     call create_constant_1d_arrays(count,Y_tmp,"Y",finab_group_id)
     call create_constant_1d_arrays(count,Y_tmp*A_tmp,"X",finab_group_id)

     ! Deallocate again
     deallocate(Y_tmp)
     deallocate(A_tmp)
     deallocate(Y_tmpA)

     ! Close the group
     call h5gclose_f(finab_group_id    , i_stat)

     ! Create group for finabsum
     call h5gcreate_f(file_id, "/finab/finabelem", finab_group_id, i_stat)
     if (i_stat .ne. 0) then
      call raise_exception("Unable to create finabelem group.",&
                           "write_finab",240005)
     end if

     ! Prepare arrays to sum over Z
     mval = maxval(isotope%p_nr)+1 ! To also write neutrons add 1
     allocate(Y_tmpA(mval))
     allocate(X_tmpA(mval))

     ! Sum up abundances
     Y_tmpA(:) = 0
     X_tmpA(:) = 0
     do i=1,net_size
       Y_tmpA(isotope(i)%p_nr+1)=Y_tmpA(isotope(i)%p_nr+1)+Y(i)
       X_tmpA(isotope(i)%p_nr+1)=X_tmpA(isotope(i)%p_nr+1)+Y(i)*isotope(i)%mass
     end do

     count = 0
     do i=1,mval
       if (Y_tmpA(i) .gt. 1d-20) then
          count = count+1
       end if
     end do

     ! Allocate
     allocate(Y_tmp(count))
     allocate(X_tmp(count))
     allocate(Z_tmp(count))

     ! Write everything to helper variables
     Y_tmp(:) = 0
     X_tmp(:) = 0
     Z_tmp(:) = 0
     count = 0
     do i=1,mval
      if (Y_tmpA(i) .gt. 1d-20) then
          count = count+1
          Y_tmp(count) = Y_tmpA(i)
          X_tmp(count) = X_tmpA(i)
          Z_tmp(count) = i-1
      end if
     end do

     ! Write the arrays
     call create_constant_1d_int_arrays(count,Z_tmp,"Z",finab_group_id)
     call create_constant_1d_arrays(count,Y_tmp,"Y",finab_group_id)
     call create_constant_1d_arrays(count,X_tmp,"X",finab_group_id)

     ! Close the group
     call h5gclose_f(finab_group_id    , i_stat)

   end subroutine



   !> @brief Gets a dataset name and creates it
   !!
   !! Creates a 1d dataset in the hdf5 and returns its
   !! dataset ID.
   !!
   !! @author Moritz Reichert
   !! @date 28.01.21
   function create_1d_dataset(dsetname,group_id)
      implicit none
      character(len=*),intent(in)       :: dsetname           !< Name of the dataset
      integer(HID_T),intent(in)         :: group_id           !< Group identifier
      integer(HSIZE_T)                  :: create_1d_dataset  !< Dataset ID of the entry
      integer(HID_T)                    :: dset_id            !< ID of the dataset
      integer(HSIZE_T), dimension(1)    :: maxdims_1d         !< Maximum (1d) dimensions
      integer                           :: i_stat             !< Reading status

      ! Set the dimension to unlimited to be able to extend to it
      maxdims_1d = H5S_UNLIMITED_F

      ! Create an appropriate dataspace
      call h5screate_simple_f(1, dims_1d, dataspace, i_stat, maxdims_1d)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create the dataspace.",&
                              "create_1d_dataset",240006)
      end if

      ! Modify dataset properties
      call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to modify dataspace properties.",&
                              "create_1d_dataset",240010)
      end if

      call h5pset_chunk_f(crp_list, 1, dims_1d, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to enable chunking.",&
                              "create_1d_dataset",240012)
      end if

      ! Create a dataset using cparms creation property
      call h5dcreate_f(group_id, dsetname, H5T_NATIVE_DOUBLE, dataspace, &
         dset_id, i_stat, crp_list )
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create the dataspace.",&
                              "create_1d_dataset",240006)
      end if
      ! Close the dataspace
      call h5sclose_f(dataspace, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to close the dataspace again.",&
                              "create_1d_dataset",240023)
      end if

      ! ! Write initial value to the dataset
      ! call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, init_val, dims_1d, i_stat)

      ! Dset ID return value
      create_1d_dataset = dset_id

   end function create_1d_dataset


   !> @brief Gets a dataset name and creates an integer dataset
   !!
   !! Creates a 1d dataset in the hdf5 and returns its
   !! dataset ID.
   !!
   !! @author Moritz Reichert
   !! @date 28.01.21
   function create_1d_int_dataset(dsetname,group_id)
      implicit none
      character(len=*),intent(in)       :: dsetname              !< Name of the dataset
      integer(HID_T),intent(in)         :: group_id              !< Group identifier
      integer(HSIZE_T)                  :: create_1d_int_dataset !< Dataset ID of the entry
      integer(HID_T)                    :: dset_id               !< ID of the dataset
      integer(HSIZE_T), dimension(1)    :: maxdims_1d            !< Maximum (1d) dimensions
      integer                           :: i_stat                !< Reading status

      ! Set the dimension to unlimited to be able to extend to it
      maxdims_1d = H5S_UNLIMITED_F

      ! Create an appropriate dataspace
      call h5screate_simple_f(1, dims_1d, dataspace, i_stat, maxdims_1d)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create the dataspace.",&
                              "create_1d_int_dataset",240006)
      end if

      ! Modify dataset properties
      call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to modify dataspace properties.",&
                              "create_1d_int_dataset",240010)
      end if

      call h5pset_chunk_f(crp_list, 1, dims_1d, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to enable chunking.",&
                              "create_1d_int_dataset",240012)
      end if


      ! Create a dataset using cparms creation property
      call h5dcreate_f(group_id, dsetname, H5T_NATIVE_INTEGER, dataspace, &
         dset_id, i_stat, crp_list )
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create the dataspace.",&
                              "create_1d_int_dataset",240006)
      end if

      call h5sclose_f(dataspace, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to close the dataspace again.",&
                              "create_1d_int_dataset",240007)
      end if

      ! Dset ID return value
      create_1d_int_dataset = dset_id

   end function create_1d_int_dataset


   !> @brief Create a hdf5 entry with constant values
   !!
   !! This is used for mass number, proton number, and
   !! neutron number.
   !!
   !! @author Moritz Reichert
   !! @date 28.01.21
   subroutine create_constant_1d_int_arrays(length,data,dsetname,group_id)
      implicit none
      integer,intent(in)                        :: length      !< ID of the dataset
      integer,dimension(length),intent(in)      :: data        !< Data values
      character(len=*),intent(in)               :: dsetname    !< Name of the dataset
      integer(HID_T)                            :: dset_id     !< ID of the dataset
      integer(HID_T)                            :: group_id    !< ID of the group
      integer(HSIZE_T), dimension(1)            :: maxdims_1d  !< Maximum (1d) dimensions
      integer(HSIZE_T), dimension(1)            :: dims_static !< Maximum (1d) dimensions
      integer                                   :: i_stat      !< Reading status

      ! Set the dimension to unlimited to be able to extend to it
      dims_static = length

      ! Create an appropriate dataspace
      call h5screate_simple_f(1, dims_static, dataspace, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create the dataspace.",&
                              "create_constant_1d_int_arrays",240006)
      end if

      ! Create a dataset using cparms creation property
      call h5dcreate_f(group_id, dsetname, H5T_NATIVE_INTEGER, dataspace, &
         dset_id, i_stat )
      if (i_stat .ne. 0) then
         call raise_exception("Unable to modify dataspace properties.",&
                              "create_constant_1d_int_arrays",240010)
      end if

      call h5sclose_f(dataspace, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to close the dataspace again.",&
                              "create_constant_1d_int_arrays",240007)
      end if
      ! Write initial value to the dataset
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims_static, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to write to the dataset.",&
                              "create_constant_1d_int_arrays",240013)
      end if

      ! Close the dataset
      call h5dclose_f(dset_id, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to close the dataset.",&
                              "create_constant_1d_int_arrays",240007)
      end if

   end subroutine create_constant_1d_int_arrays


   !> @brief Create a hdf5 entry with constant values
   !!
   !! This is used for e.g., the flows
   !!
   !! @author Moritz Reichert
   !! @date 28.01.21
   subroutine create_constant_1d_arrays(length,data,dsetname,group_id)
      implicit none
      integer,intent(in)                        :: length      !< ID of the dataset
      real(r_kind),dimension(length),intent(in) :: data        !< Data values
      character(len=*),intent(in)               :: dsetname    !< Name of the dataset
      integer(HID_T)                            :: dset_id     !< ID of the dataset
      integer(HID_T)                            :: group_id    !< ID of the group
      integer(HSIZE_T), dimension(1)            :: maxdims_1d  !< Maximum (1d) dimensions
      integer(HSIZE_T), dimension(1)            :: dims_static !< Maximum (1d) dimensions
      integer                                   :: i_stat      !< Reading status

      ! Set the dimension to unlimited to be able to extend to it
      dims_static = length

      ! Create an appropriate dataspace
      call h5screate_simple_f(1, dims_static, dataspace, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create the dataspace.",&
                              "create_constant_1d_int_arrays",240006)
      end if

      ! Create a dataset using cparms creation property
      call h5dcreate_f(group_id, dsetname, H5T_NATIVE_DOUBLE, dataspace, &
         dset_id, i_stat )
      if (i_stat .ne. 0) then
         call raise_exception("Unable to modify dataspace properties.",&
                              "create_constant_1d_int_arrays",240010)
      end if

      call h5sclose_f(dataspace, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to close the dataspace again.",&
                              "create_constant_1d_int_arrays",240007)
      end if
      ! Write initial value to the dataset
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims_static, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to write to the dataset.",&
                              "create_constant_1d_int_arrays",240013)
      end if

      ! Close the dataset
      call h5dclose_f(dset_id, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to close the dataset.",&
                              "create_constant_1d_int_arrays",240007)
      end if

   end subroutine create_constant_1d_arrays


   !> @brief Extend a previously created 1D - dataset
   !!
   !! The necessary input is the value and the
   !! id of the dataset only.
   !!
   !! @author Moritz Reichert
   !! @date 28.01.21
   subroutine extend_1d_dataset(dset_id,value,in_size,chunksize)
      use global_class, only: net_size
      implicit none
      integer,intent(in)                           :: chunksize!< Size of the amount of variables written
      real(r_kind),dimension(chunksize),intent(in) :: value    !< New value that will get appended
      integer(HID_T),intent(in)                    :: dset_id  !< ID of the dataset
      integer,intent(in)                           :: in_size  !< New size of the array
      integer(HID_T)                               :: memspace !< Memory dataspace identifier
      integer                                      :: i_stat   !< Reading status
      integer(HSIZE_T), dimension(1)               :: offset   !< Offset, dont overwrite the old data
      integer(HSIZE_T), dimension(1)               :: count    !< New line
      integer(HSIZE_T), dimension(1)               :: new_size !< New total size of the dataset

      !Extend the dataset (new dimensions)
      new_size(1)   = in_size+chunksize
      dims_1d = (/chunksize/)
      call h5dset_extent_f(dset_id, new_size, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to extend the dataset.",&
                              "extend_1d_dataset",240014)
      end if

      call h5screate_simple_f (1, dims_1d, memspace, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create new memspace.",&
                              "extend_1d_dataset",240015)
      end if

      ! Offset by previous size
      offset(1) = in_size
      ! write a new row
      count(1)  = chunksize

      ! Create space, select hyperslab
      call h5dget_space_f(dset_id, dataspace, i_stat)
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          offset, count, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to select hyperslab.",&
                              "extend_1d_dataset",240016)
      end if
      ! Write to the new row
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, value, dims_1d, i_stat, &
         memspace, dataspace)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to write to the file.",&
                              "extend_1d_dataset",240013)
      end if
   end subroutine extend_1d_dataset


   !> @brief Extend a previously created 1D integer - dataset
   !!
   !! The necessary input is the value and the
   !! id of the dataset only.
   !!
   !! @author Moritz Reichert
   !! @date 28.01.21
   subroutine extend_1d_int_dataset(dset_id,value,in_size,chunksize)
      use global_class, only: net_size
      implicit none
      integer,intent(in)             :: chunksize!< Size of the amount of variables written
      integer,dimension(chunksize),intent(in)    :: value    !< New value that will get appended
      integer(HID_T),intent(in)      :: dset_id  !< ID of the dataset
      integer,intent(in)             :: in_size  !< New size of the array
      integer(HID_T)                 :: memspace !< Memory dataspace identifier
      integer                        :: i_stat   !< Reading status
      integer(HSIZE_T), dimension(1) :: offset   !< Offset, dont overwrite the old data
      integer(HSIZE_T), dimension(1) :: count    !< New line
      integer(HSIZE_T), dimension(1) :: new_size !< New total size of the dataset

      !Extend the dataset (new dimensions)
      new_size(1)   = in_size+chunksize
      dims_1d = (/chunksize/)
      call h5dset_extent_f(dset_id, new_size, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to extend the dataset.",&
                              "extend_1d_int_dataset",240014)
      end if

      call h5screate_simple_f (1, dims_1d, memspace, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to create new memspace.",&
                              "extend_1d_int_dataset",240015)
      end if

      ! Offset one row
      offset(1) = in_size
      ! write a new row
      count(1)  = chunksize

      ! Create space, select hyperslab
      call h5dget_space_f(dset_id, dataspace, i_stat)
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          offset, count, i_stat)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to select hyperslab.",&
                              "extend_1d_int_dataset",240016)
      end if

      ! Write to the new row
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, value, dims_1d, i_stat, &
         memspace, dataspace)
      if (i_stat .ne. 0) then
         call raise_exception("Unable to write to the file.",&
                              "extend_1d_int_dataset",240013)
      end if
   end subroutine extend_1d_int_dataset



   !>
   !! @brief Subroutine to write the average timescales (1/s) into the hdf5
   !!
   !! @author Moritz Reichert
   !! @date 5.02.21
   subroutine extend_av_timescales(time,tau_ga, tau_ag,tau_ng, tau_gn, tau_pg, tau_gp, tau_np, &
                                   tau_pn, tau_an, tau_na, tau_ap, tau_pa, tau_beta, tau_alpha, &
                                   tau_nfiss, tau_sfiss, tau_bfiss)
      implicit none
      real(r_kind),intent(in)  :: time                                        !< Time [s]
      real(r_kind),intent(in)  :: tau_ng, tau_gn, tau_pg, tau_gp, tau_np      !< Timescales [1/s]
      real(r_kind),intent(in)  :: tau_ga, tau_ag, tau_ap, tau_pa              !< Timescales [1/s]
      real(r_kind),intent(in)  :: tau_pn, tau_an, tau_na, tau_beta, tau_alpha !< Timescales [1/s]
      real(r_kind),intent(in)  :: tau_nfiss, tau_sfiss, tau_bfiss             !< Timescales [1/s]
      !
      integer                   :: i  !< Loop variable

      ! Count how large the space has to be
      chunk_counter_ts = chunk_counter_ts + 1
      !< Save the data temporary in an array
      chunk_store_ts(chunk_counter_ts,:) = (/time, tau_ga, tau_ag, tau_ng, tau_gn, tau_pg, tau_gp, tau_np, &
                                             tau_pn, tau_an, tau_na, tau_ap, tau_pa, tau_beta, tau_alpha, &
                                             tau_nfiss, tau_sfiss, tau_bfiss/)
      ! Write the chunk to the file
      if (chunk_counter_ts .eq. chunk_size_ts) then
         call write_av_timescales()
      end if

   end subroutine extend_av_timescales

   !>
   !! @brief Write average timescales and time into the hdf5 file
   !!
   !! These values are stored in a separate group (called timescales) that is initialized
   !! in \ref init_hdf5_module.
   !!
   !! @author Moritz Reichert
   !! @date 5.02.21
   subroutine write_av_timescales()
      implicit none
      integer :: i !< Loop variable

      if (chunk_counter_ts .lt. 1) return
      ! Save the time
      call extend_1d_dataset(ts_dset_id_t,chunk_store_ts(1:chunk_counter_ts,1),iter_ts,chunk_counter_ts)
      ! Save timescales
      do i=1, 17
         call extend_1d_dataset(ts_reac_ids(i),chunk_store_ts(1:chunk_counter_ts,i+1),iter_ts,chunk_counter_ts)
      end do
      iter_ts = iter_ts+chunk_counter_ts
      chunk_counter_ts = 0
   end subroutine write_av_timescales

   !>
   !! @brief Write tracked abundances and time into the hdf5 file
   !!
   !! These values are stored in a separate group (called tracked_nuclei) that is initialized
   !! in \ref init_hdf5_module.
   !!
   !! @author Moritz Reichert
   !! @date 5.02.21
   subroutine extend_track_nuclei(time,Y_array)
      use global_class, only: track_nuclei_nr,track_nuclei_indices,net_size
      implicit none
      real(r_kind),dimension(net_size),intent(in) :: Y_array !< Abundance to store
      real(r_kind),intent(in)                     :: time    !< Time [s]
      integer                                     :: i       !< Loop variable

      ! Count the amount of entries in an intermediate array
      chunk_counter_track = chunk_counter_track+1

      ! Store the time
      chunk_store_track(1,chunk_counter_track) = time
      ! Store abundances
      do i=1, track_nuclei_nr
         chunk_store_track(i+1,chunk_counter_track) = Y_array(track_nuclei_indices(i))
      end do

      if (chunk_counter_track .eq. chunk_size_track) then
         ! Extend the datasets
         call write_track_data()
      end if

   end subroutine extend_track_nuclei



   !> Write the abundances of track nuclei into a file
   !!
   !! When this subroutine is called depends on the chunk size.
   !!
   !! @author Moritz Reichert
   !! @date 09.02.21
   subroutine write_track_data()
      use global_class, only:track_nuclei_nr
      implicit none
      integer :: i !< Loop variable

      ! Do nothing if it recently got written already
      if (chunk_counter_track .lt. 1) return

      ! write time array by "chunk" timesteps
      call extend_1d_dataset(tracked_dset_id_t,chunk_store_track(1,1:chunk_counter_track),iter_tracked,chunk_counter_track)
      ! Write abundances
      call extend_hdf5_Y(tracked_abu_id,chunk_store_track(2:track_nuclei_nr+1,1:chunk_counter_track),track_nuclei_nr,iter_tracked,chunk_counter_track)
      iter_tracked = iter_tracked + chunk_counter_track
      chunk_counter_track = 0
   end subroutine write_track_data




   !>
   !! Write abundances and time into the hdf5 file
   !!
   !! These values are stored in a separate group (called snapshots) that is initialized
   !! in \ref init_hdf5_module.
   !!
   !! @author Moritz Reichert
   !! @date 5.02.21
   subroutine extend_snaps(time,Y_array)
      use global_class, only: net_size
      implicit none
      real(r_kind),dimension(net_size),intent(in) :: Y_array !< Abundance to store
      real(r_kind),intent(in)                     :: time    !< Time [s]

      ! Store the time
      chunk_counter_snaps = chunk_counter_snaps + 1
      chunk_store_snaps_time(chunk_counter_snaps) = time
      ! Save the abundances in intermediate arrays
      chunk_store_snaps_Y(:,chunk_counter_snaps) = Y_array

      ! Write to the file
      if (chunk_counter_snaps .eq. chunk_size_snaps) then
         call write_snaps()
      end if
   end subroutine extend_snaps


   !> Write the abundances of all nuclei into a file
   !!
   !! When this subroutine is called depends on the chunk size.
   !!
   !! @author Moritz Reichert
   !! @date 09.02.21
   subroutine write_snaps()
      use global_class, only: net_size
      implicit none

      ! Do nothing if it was recently called already
      if (chunk_size_snaps .lt. 1) return

      ! extend time array by one timestep
      call extend_1d_dataset(snaps_dset_id_t,chunk_store_snaps_time(1:chunk_counter_snaps),iter_snaps,chunk_counter_snaps)
      call extend_hdf5_Y(dset_id_Y,chunk_store_snaps_Y(:,1:chunk_counter_snaps),net_size,iter_snaps,chunk_counter_snaps)
      iter_snaps = iter_snaps + chunk_counter_snaps
      chunk_counter_snaps = 0 !< reset the chunk counter

   end subroutine write_snaps



   !> Write the mainout data to the hdf5
   !!
   !! This will write the iteration, time (s),
   !! temperature (GK), density (g/ccm), entropy(kB/nuc),
   !! radius (km), Ye, Yn, Yp, Ya, Ylight, Yheavy,
   !! abar, and zbar.
   !!
   !! @author Moritz Reichert
   !! @date 5.02.21
   subroutine extend_mainout(cnt,time,temp,dens,entr,rad,Y)
      use global_class, only: net_size,ipro,ihe4,ineu,isotope
      implicit none
      integer,intent(in)                          :: cnt     !< Iteration count
      real(r_kind),intent(in)                     :: time    !< Time [s]
      real(r_kind),intent(in)                     :: temp    !< Temperature [GK]
      real(r_kind),intent(in)                     :: dens    !< Density [g/ccm]
      real(r_kind),intent(in)                     :: entr    !< Entropy [kB/nuc]
      real(r_kind),intent(in)                     :: rad     !< Radius [km]
      real(r_kind),dimension(net_size),intent(in) :: Y       !< Abundances
      ! Helper variables
      real(r_kind)      :: Ye,Yneu,Ypro,Yhe4,ylight,yheavies,ysum,yasum,yzsum,abar,zbar

      ! Count the amount of entries in an intermediate array
      chunk_counter_mainout = chunk_counter_mainout+1

      ! output diagnostics
      Ye = sum(isotope(1:net_size)%p_nr*Y(1:net_size))
      ylight=   max(1.0d-99,sum(Y(ipro+1:ihe4-1)))
      yheavies= 0.0
      if (ihe4.ge.1) yheavies= max(1.0d-99,sum(Y(ihe4+1:net_size)))
      ! Calculate abar and zbar
      ysum  = sum(Y(1:net_size))
      yasum = sum(Y(1:net_size)*isotope(1:net_size)%mass)
      yzsum = sum(Y(1:net_size)*isotope(1:net_size)%p_nr)
      abar  = yasum/ysum
      zbar  = yzsum/ysum
      ! Get abundance of neutrons, protons, and alphas if they exist in the network
      Yneu= 0.0
      if (ineu.ge.1) Yneu= max(1.0d-99,Y(ineu))

      Ypro= 0.0
      if (ipro.ge.1) Ypro= max(1.0d-99,Y(ipro))

      Yhe4= 0.0
      if (ihe4.ge.1) Yhe4= max(1.0d-99,Y(ihe4))

      ! Store data in intermediate arrays to prevent to much IO going on
      chunk_store_mainout(chunk_counter_mainout,1) = time
      chunk_store_mainout(chunk_counter_mainout,2) = temp
      chunk_store_mainout(chunk_counter_mainout,3) = dens
      chunk_store_mainout(chunk_counter_mainout,4) = entr
      chunk_store_mainout(chunk_counter_mainout,5) = rad
      chunk_store_mainout(chunk_counter_mainout,6) = Ye
      chunk_store_mainout(chunk_counter_mainout,7) = Yneu
      chunk_store_mainout(chunk_counter_mainout,8) = Ypro
      chunk_store_mainout(chunk_counter_mainout,9) = Yhe4
      chunk_store_mainout(chunk_counter_mainout,10)= abar
      chunk_store_mainout(chunk_counter_mainout,11)= zbar
      chunk_store_mainout(chunk_counter_mainout,12)= ylight
      chunk_store_mainout(chunk_counter_mainout,13)= yheavies
      ! Same for iterations
      chunk_store_int_mainout(chunk_counter_mainout)= cnt

      if (chunk_counter_mainout .eq. chunk_size_mainout) then
         ! Extend the datasets
         call write_mainout_data()
      end if
   end subroutine extend_mainout



   !> Write the mainout data into a file
   !!
   !! When this subroutine is called depends on the chunk size.
   !!
   !! @author Moritz Reichert
   !! @date 09.02.21
   subroutine write_mainout_data()
      implicit none
      ! Do nothing if it recently got written already
      if (chunk_counter_mainout .lt. 1) return
      ! Extend the datasets
      call extend_1d_int_dataset(mainout_dset_id_iter  ,&
                                 chunk_store_int_mainout(1:chunk_counter_mainout),&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_t     ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,1)  ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_temp  ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,2)  ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_dens  ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,3)  ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_entr  ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,4)  ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_rad   ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,5)  ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_ye    ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,6)  ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_yn    ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,7)  ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_yp    ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,8)  ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_ya    ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,9)  ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_abar  ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,10) ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_zbar  ,&
                                 chunk_store_mainout(1:chunk_counter_mainout,11) ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_ylight,&
                                 chunk_store_mainout(1:chunk_counter_mainout,12) ,&
                                 iter_mainout, chunk_counter_mainout)
      call extend_1d_dataset(    mainout_dset_id_yheavy,&
                                 chunk_store_mainout(1:chunk_counter_mainout,13) ,&
                                 iter_mainout, chunk_counter_mainout)
      iter_mainout = iter_mainout + chunk_size_mainout
      chunk_counter_mainout = 0 !< reset the chunk counter
   end subroutine write_mainout_data



   !> @brief Extend the Y dataset for the net timestep.
   !!
   !! The dimension of this is given by global_class::net_size and the current
   !! iteration \ref iter.
   !!
   !! @see extend_1d_dataset, extend_1d_int_dataset
   !! @author Moritz Reichert
   !! @date 28.01.21
   subroutine extend_hdf5_Y(dset_id,Y_array,nr_nuc,old_size,chunksize)
      use global_class, only: net_size
      implicit none
      integer(HID_T),intent(in)                           :: dset_id  !< Dataset identifier
      integer,intent(in)                                  :: nr_nuc   !< Amount of abundances
      integer,intent(in)                                  :: old_size !< Old length of arrays
      integer,intent(in)                                  :: chunksize!< Size of the chunk
      real(r_kind),dimension(nr_nuc,chunksize),intent(in) :: Y_array  !< Abundance array to store
      integer(HSIZE_T),dimension(2)                       :: new_size !< New size of dataset
      integer(HSIZE_T),dimension(1:2)                     :: offset   !< offset of the data in the file
      integer(HSIZE_T),dimension(1:2)                     :: count    !< Size of the new data
      integer                                             :: i_stat   !< status variable
      integer(HID_T)                                      :: memspace !< Memory dataspace identifier

      !Extend the dataset.
      new_size(1:2)   = (/nr_nuc,old_size+chunksize/)
      ! Set the dimensions
      dims_Y = (/nr_nuc,chunksize/)

      call h5dset_extent_f(dset_id, new_size, i_stat)
      call h5screate_simple_f (2, dims_Y, memspace, i_stat)

      ! Offset by the old size
      offset(1:2) = (/0,old_size/)
      ! write a couple of new rows
      count(1:2)  = (/nr_nuc,chunksize/)

      ! Create space, select hyperslab
      call h5dget_space_f(dset_id, dataspace, i_stat)
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          offset, count, i_stat)
      ! Write to the new row
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Y_array, dims_Y, i_stat, &
         memspace, dataspace)

   end subroutine extend_hdf5_Y


   !> Finalize the hdf5 module
   !!
   !! This subroutine closes all datasets, dataspaces, files
   !! and other things.
   !!
   !! @author Moritz Reichert
   !! @date 28.01.21
   subroutine hdf5_module_finalize()
      use parameter_class, only: h_snapshot_every,h_mainout_every,&
                                 h_custom_snapshots,h_track_nuclei_every,&
                                 h_timescales_every,h_flow_every,&
                                 h_engen_every,h_engen_detailed,&
                                 h_nu_loss_every
      use global_class,    only: track_nuclei_nr
      implicit none
      integer :: i_stat        !< Status variable
      integer :: i             !< Loop variable

      ! Don't close anything if there is no file
      if (.not. hdf5_output) return

      ! Cleanup snapshots
      if ((h_snapshot_every .gt. 0) .or. (h_custom_snapshots)) then
         call write_snaps() ! Write the last data into the file
         call h5gclose_f(snaps_group_id,i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close snapshot group.",&
                                                 "hdf5_module_finalize",240017)
         ! Close all datasets
         call h5dclose_f(dset_id_Y   , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close abundance dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(snaps_dset_id_t   , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close time dataset.",&
                                                 "hdf5_module_finalize",240007)
      end if

      ! Cleanup mainout
      if (h_mainout_every .gt. 0) then
         call write_mainout_data() ! Write the last data into the file
         call h5gclose_f(mainout_group_id      , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close mainout group.",&
                                                 "hdf5_module_finalize",240017)
         call h5dclose_f(mainout_dset_id_t     , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close time dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_temp  , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close temperature dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_dens  , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close density dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_entr  , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close entropy dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_rad   , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close radius dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_ye    , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close electron fraction dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_iter  , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close iteration dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_ylight, i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close ylight dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_yheavy, i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close yheavy dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_yn    , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close yn dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_yp    , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close yp dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(mainout_dset_id_ya    , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close ya dataset.",&
                                                 "hdf5_module_finalize",240007)
      end if

      ! Cleanup tracked nuclei
      if (h_track_nuclei_every .gt. 0) then
         call write_track_data() ! Write the last datapoints
         call h5gclose_f(track_group_id    , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close track_nuclei group.",&
                                                 "hdf5_module_finalize",240017)
         call h5dclose_f(tracked_dset_id_t , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close time dataset.",&
                                                 "hdf5_module_finalize",240007)
         call h5dclose_f(tracked_abu_id , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close tracked abundance dataset.",&
                                                 "hdf5_module_finalize",240007)
      end if

      ! Cleanup neutrino loss/gain
      if ((h_nu_loss_every .gt. 0)) then
        call write_nu_loss_gain()
        call h5gclose_f(nuloss_group_id  , i_stat)
        if (i_stat .ne. 0) call raise_exception("Unable to close nuloss group.",&
                                                "hdf5_module_finalize",240017)
        call h5dclose_f(nuloss_dset_id_t , i_stat)
        if (i_stat .ne. 0) call raise_exception("Unable to close time dataset.",&
                                                "hdf5_module_finalize",240007)
        call h5dclose_f(nuloss_dset_id_bet , i_stat)
        if (i_stat .ne. 0) call raise_exception("Unable to close beta dataset.",&
                                                "hdf5_module_finalize",240007)
        call h5dclose_f(nuloss_dset_id_dens , i_stat)
        if (i_stat .ne. 0) call raise_exception("Unable to close density dataset.",&
                                                "hdf5_module_finalize",240007)
        call h5dclose_f(nuloss_dset_id_heat , i_stat)
        if (i_stat .ne. 0) call raise_exception("Unable to close heat dataset.",&
                                                "hdf5_module_finalize",240007)
        call h5dclose_f(nuloss_dset_id_nut , i_stat)
        if (i_stat .ne. 0) call raise_exception("Unable to close nu total dataset.",&
                                                "hdf5_module_finalize",240007)
        call h5dclose_f(nuloss_dset_id_rad , i_stat)
        if (i_stat .ne. 0) call raise_exception("Unable to close radius dataset.",&
                                                "hdf5_module_finalize",240007)
        call h5dclose_f(nuloss_dset_id_temp , i_stat)
        if (i_stat .ne. 0) call raise_exception("Unable to close temperature dataset.",&
                                                "hdf5_module_finalize",240007)
        call h5dclose_f(nuloss_dset_id_the , i_stat)
        if (i_stat .ne. 0) call raise_exception("Unable to close thermal dataset.",&
                                                "hdf5_module_finalize",240007)
      end if

      ! Cleanup average timescales
      if (h_timescales_every .gt. 0) then
         call write_av_timescales() ! Write the last datapoints
         call h5gclose_f(ts_group_id  , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close timescale group.",&
                                                 "hdf5_module_finalize",240017)
         call h5dclose_f(ts_dset_id_t , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close time dataset.",&
                                                 "hdf5_module_finalize",240007)
         do i=1,17
            call h5dclose_f(ts_reac_ids(i) , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close timescale dataset.",&
                                                    "hdf5_module_finalize",240007)
         end do
      end if

      ! Cleanup energy generation output
      if (h_engen_every .gt. 0) then
         call write_engen(h_engen_detailed) ! Write the last datapoints
         call h5gclose_f(en_group_id  , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close energy group.",&
                                                 "hdf5_module_finalize",240017)
         call h5dclose_f(en_dset_id_t , i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close time dataset.",&
                                                 "hdf5_module_finalize",240007)
         do i=1,10
            call h5dclose_f(en_reac_ids(i) , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close engen dataset.",&
                                                    "hdf5_module_finalize",240007)
         end do
         if (h_engen_detailed .eqv. .true.) then
            call h5dclose_f(en_det_bet_id , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close decay energy dataset.",&
                                                    "hdf5_module_finalize",240007)
            call h5dclose_f(en_det_ag_id , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close (a,g) energy dataset.",&
                                                    "hdf5_module_finalize",240007)
            call h5dclose_f(en_det_pg_id , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close (p,g) energy dataset.",&
                                                    "hdf5_module_finalize",240007)
            call h5dclose_f(en_det_ng_id , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close (n,g) energy dataset.",&
                                                    "hdf5_module_finalize",240007)
            call h5dclose_f(en_det_pn_id , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close (p,n) energy dataset.",&
                                                    "hdf5_module_finalize",240007)
            call h5dclose_f(en_det_ap_id , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close (a,p) energy dataset.",&
                                                    "hdf5_module_finalize",240007)
            call h5dclose_f(en_det_an_id , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close (a,n) energy dataset.",&
                                                    "hdf5_module_finalize",240007)
            call h5dclose_f(en_det_other_id , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close other energy dataset.",&
                                                    "hdf5_module_finalize",240007)
            call h5dclose_f(en_det_fiss_id , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close fiss energy dataset.",&
                                                    "hdf5_module_finalize",240007)
            call h5dclose_f(en_det_tot_id , i_stat)
            if (i_stat .ne. 0) call raise_exception("Unable to close tot energy dataset.",&
                                                    "hdf5_module_finalize",240007)
         end if

      end if

      ! Cleanup flow
      if (h_flow_every .gt. 0) then
         call h5gclose_f(flow_group_id, i_stat)
         if (i_stat .ne. 0) call raise_exception("Unable to close flow group.",&
                                                 "hdf5_module_finalize",240017)
      end if

      ! Close all other helper variables
      call h5fclose_f(file_id  , i_stat)
      if (i_stat .ne. 0) call raise_exception("Unable to close hdf5 file.",&
                                              "hdf5_module_finalize",240018)

      call h5pclose_f(crp_list , i_stat)
      if (i_stat .ne. 0) call raise_exception("Unable to close property list.",&
                                              "hdf5_module_finalize",240019)
   end subroutine hdf5_module_finalize

#endif

end module hdf5_module
