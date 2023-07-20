!> @file flow_module.f90
!!
!! The error file code for this file is ***W20***.
!! @brief Module \ref flow_module for calculating reaction flows
!!

!> Provides subroutines to calculate reaction flows
!!
!! @author  Christian Winteler
!! @date    07.10.10
!!
!! \b Edited:
!!           - 03.04.18, M. Jacobi   , Rewrote the module, the flows are now
!!                                     calculated with the help of the Jacobian
!!           - 22.01.21, M. Reichert , added more comments
!! .
#include "macros.h"
module flow_module
  use global_class,    only: net_size, ihe4, ineu, ipro, flow_vector
  use pardiso_class,   only: jind, vals
  implicit none



  type(flow_vector),dimension(:),allocatable, public  :: flows

  !
  ! Public and private fields and methods of the module
  !
  public:: &
      flow_init, flowcalc, flowprint

  private:: &
      flowsort

contains



!>
!! Initialise flow subroutine
!!
!! This subroutine counts the number of possible flows
!! and allocates the \ref flows array.
!!
!! \b Edited:
!!          - 11.01.14
!!          - 03.04.18, M. Jacobi
subroutine flow_init()
   implicit none

   integer                               :: ind, i, j

   INFO_ENTRY("flow_init")

   ind=0

   ! Cycle through Jacobian and look at J_ij and J_ji
   ! Get amount of flow instances
   do i=1,net_size
      ! Ignore hydrogen, neutrons and helium
      if (i.eq.ihe4) cycle
      if (i.eq.ipro) cycle
      if (i.eq.ineu) cycle

      do j=1,net_size
         if (j.eq.ihe4) cycle
         if (j.eq.ipro) cycle
         if (j.eq.ineu) cycle

         ! Ignore diagonal
         if (i.eq.j)    cycle
         ! Account for entries that are only present at one of the diagonals
         if (jind(i,j).eq.0) cycle
         if (i .lt. j) then
            if (.not. jind(j,i).eq.0) cycle
         end if

         ! Count amount of possible flows
         ind = ind +1

      end do
   end do

   allocate (flows(ind))

   INFO_EXIT("flow_init")

   return

end subroutine flow_init



!>
!! Flow calculation from jacobian. It is calculated with the help of the Jacobian.
!! \f[
!! F_{ij} = |(1/h - J_{ij}) \times Y_j - (1/h - J_{ji}) \times Y_i|
!! \f]
!!
!! @note Using the jacobian directly has the advantage
!!       that the flow will be correct if the calculation
!!       is correct. In previous versions, the flow
!!       was not calculated by using the jacobian.
!!
!! \b Edited:
!!          - 03.04.18, M. Jacobi
!! .
subroutine flowcalc(Y)
   implicit none

   real(r_kind),dimension(:),intent(in)  :: Y         !< abundance
   integer                               :: ind, i, j

   INFO_ENTRY("flowcalc")

   ind=0
   ! Cycle through values of the Jacobian and calculate flows
   do i=1,net_size
      ! Ignore hydrogen, neutrons and helium
      if (i.eq.ihe4) cycle
      if (i.eq.ipro) cycle
      if (i.eq.ineu) cycle

      do j=1,net_size
         if (j.eq.ihe4) cycle
         if (j.eq.ipro) cycle
         if (j.eq.ineu) cycle
         ! Ignore diagonal
         if (i.eq.j)    cycle
         ! Account for entries that are only present at one of the triangle matrix
         if (jind(i,j).eq.0) cycle
         if (i .lt. j) then
            if (.not. jind(j,i).eq.0) cycle
         end if

         ind = ind +1

         ! Store the flow i->j
         flows(ind)%iin  = i
         flows(ind)%iout = j
         flows(ind)%fwd  = -vals(jind(i,j))*Y(i)

         ! Store the flow j->i
         if (jind(j,i).eq.0) then
            flows(ind)%bwd = 0
         else
            flows(ind)%bwd = -vals(jind(j,i))*Y(j)
         end if
      end do
   end do

   ! Swap fwd and bwd flows if necessary
   call flowsort()

   INFO_EXIT("flowcalc")

   return

end subroutine flowcalc


!> Bring flows to correct format
!!
!! Make the flows positive and swap ingoing neutron and proton numbers
!! in case of negative flows.
subroutine flowsort
   implicit none
   integer      :: i    !< Loop variable
   integer      :: tmp  !< Temporary helper variable for an index
   real(r_kind) :: ftmp !< Temporary helper variable for a flow

   do i = 1,size(flows)
      if(flows(i)%bwd.gt.flows(i)%fwd) then
         ftmp = flows(i)%fwd
         flows(i)%fwd = flows(i)%bwd
         flows(i)%bwd = ftmp
         tmp = flows(i)%iin
         flows(i)%iin = flows(i)%iout
         flows(i)%iout = tmp
      end if

   end do

end subroutine flowsort

!>
!! Output reaction flows to a file
!!
!! An example of this file may look like:
!!\file{
!! time    temp    dens
!! 1.03895957612263E-01   7.19136097013393E+00   1.40977753502083E+06
!!  nin     zin     yin    nout    zout    yout    flow
!! 2   1   4.81807892321990E-08   1   1   2.13994533749120E-06   0.00000000000000E+00
!! 1   2   1.26489216252989E-09   1   1   2.13994533749120E-06   0.00000000000000E+00
!! 1   2   1.26489216252989E-09   2   1   4.81807892321990E-08   1.58426675189734E-10
!! 4   2   9.86495465952053E-13   3   3   2.15833022688002E-11   8.53665754802155E-13
!! ...}
!!
!! \b Edited:
!!         - 11.01.14
!!         - 03.04.18, M. Jacobi
!! .
subroutine flowprint(t,t9,dens,abu,cnt)
   use global_class, only: isotope_type, isotope
   use file_handling_class
   implicit none

   real(r_kind),intent(in)  :: t                !< time [s]
   real(r_kind),intent(in)  :: t9               !< temperature [GK]
   real(r_kind),intent(in)  :: dens             !< density [g/cm3]
   real(r_kind),dimension(:),intent(in)  :: abu !< abundances
   integer,intent(in)       :: cnt              !< flow snapshot counter
   !
   type(isotope_type)    :: nucin,nucout
   integer               :: i
   integer               :: ini,ino
   integer               :: flowunit
   character*50          :: flowfile
   real(r_kind)          :: Y_ini, Y_ino, flow_diff

   INFO_ENTRY("flowprint")

   select case (cnt)
   case(:9)
      write(flowfile,'(a16,i1,a4)')'flow/flow_000',cnt,'.dat'
   case(10:99)
      write(flowfile,'(a15,i2,a4)')'flow/flow_00',cnt,'.dat'
   case(100:999)
      write(flowfile,'(a14,i3,a4)')'flow/flow_0',cnt,'.dat'
   case default
      write(flowfile,'(a13,i4,a4)')'flow/flow_',cnt,'.dat'
   end select
   flowunit= open_outfile (adjustl(flowfile))

   write(flowunit,'(3a8)')'time','temp','dens'
   write(flowunit,'(3es23.14)')t,t9,dens
   write(flowunit,'(7(a8))')'nin','zin','yin','nout','zout','yout',&
                              'flow'
   do i = 1,size(flows)
      ! Declare helper variables
      ini = flows(i)%iin
      ino = flows(i)%iout
      nucin = isotope(ini)
      nucout = isotope(ino)

      ! Make a matching format and prevent something like 1.89-392 instead of 1.89E-392
      flow_diff = flows(i)%fwd - flows(i)%bwd
      Y_ino     = abu(ino)
      Y_ini     = abu(ini)
      if (abu(ini) .ne.  0.) Y_ini     = max(1.d-99,abu(ini))
      if (abu(ino) .ne.  0.) Y_ino     = max(1.d-99,abu(ino))
      if (flow_diff .ne. 0.) flow_diff = max(1.d-99,flows(i)%fwd - flows(i)%bwd)

      if (flow_diff .ne. 0) then
         write(flowunit,'(2(2i4,es23.14),es23.14)')        &
              nucin%n_nr,nucin%p_nr,Y_ini,                 &
              nucout%n_nr,nucout%p_nr,Y_ino,               &
              flow_diff
      end if
   end do

   call close_io_file(flowunit,flowfile)

   INFO_EXIT("flowprint")

   return

end subroutine flowprint

end module flow_module
