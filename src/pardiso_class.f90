!> @file pardiso_class.f90
!!
!! The error file code for this file is ***W35***.
!!
!! @brief Module \ref pardiso_class with the sparse solver
!!

!> Contains subroutines for sparse matrix assembly and
!! the solver core
!!
!! @note this module used to be split
!! into two: cscf and pardiso_class
!!
!! @author Christian Winteler
!! @date   27.10.2009
!!
!! \b Edited:
!!    - 14.01.14, Oleg Korobkin
!!    -  7.11.16, Oleg Korobkin
!!    -  6.04.18, M. Jacobi
!!    - 24.07.22, M. Reichert
!! .
#include "macros.h"
module pardiso_class
  use error_msg_class, only: raise_exception,num_to_str,int_to_str
  use global_class, only: net_size
  implicit none

  ! Variables are used in the jacobian class and the flow_module
  real(r_kind),dimension(:),allocatable, public  :: vals !< contains non-zero matrix elements (of the Jacobian)
  integer,dimension(:),allocatable     , public  :: rows !< rows(i) is the number of the row in the matrix that contains the i-th value in vals
  integer,dimension(:),allocatable     , public  :: pt_b !< pt_b(j) gives the index into vals that contains the FIRST non-zero element of column j of the matrix
  integer,dimension(:),allocatable     , public  :: pt_e !< pt_e(j)-pt_b(j) gives the index into vals that contains the LAST non-zero element of column j of the matrix
  integer,dimension(:),allocatable     , public  :: dia  !< dia(j) gives the index into vals that contains the (j,j) diagonal element of the matrix
  integer,dimension(:,:),allocatable   , public  :: jind !< jind(i,j) contains the index of jacobian entry in vals

  !
  ! Public and private fields and methods of the module.
  !
  public:: &
       netsolve, sparse
  ! private:: &
  !      nothing here

contains



!> Determines the position of jacobian entries in the cscf value array.
!!
!! This subroutine creates the basic sparse format of the jacobian. For this
!! it writes a dummy matrix with entries = 1 at places in the matrix
!! that are connected via nuclear reaction rates. Afterwards the matrix
!! is stored in compressed sparse column format (cscf) into the arrays
!! \ref vals, \ref rows, and \ref pt_b.
!! In cscf format, the matrix:
!! \f[
!! X= \left( \begin{array}{cccc}
!!          9 & 8 & 0 & 2\\
!!          0 & 0 & 0 & 5\\
!!          0 & 3 & 0 & 0\\
!!          5 & 0 & 4 & 9\end{array} \right).
!! \f]
!! will be stored like
!! \f[ \mathrm{vals} =[9,5,8,3,4,2,5,9] \f]
!! with vals storing the non-zero values of X,
!! \f[ \mathrm{rows} = [1,4,1,3,4,1,2,4] \f]
!! rows storing the position of a new row in vals, and
!! \f[ \mathrm{pt\_b} = [1,3,5,6,9] \f].
!! the position of the columns of the non-zero entries.
!!
!! @see [Winteler 2013](https://edoc.unibas.ch/29895/), section 2.4.2
!!
!! \b Edited:
!!     - 11.01.14
!!     - 28.07.22, MR: introduced new chapters
!! .
subroutine sparse
   use global_class, only: reactionrate_type, rrate, nreac, ineu
   use fission_rate_module, only: nfiss, fissrate
   use parameter_class, only: fissflag
   implicit none

   integer                 :: i,j,n   !< loop variables
   type(reactionrate_type) :: rr_tmp  !< Temporary reaction rate variable
   integer, dimension(net_size*net_size) :: ja_t
   integer, dimension(:,:),allocatable :: cscf_tmp_fiss
   integer, dimension(4,6) :: cscf_tmp
   integer                 :: cnt

   INFO_ENTRY("sparse")
   n = net_size
   allocate(pt_b(n+1))
   allocate(pt_e(n))
   allocate(dia(n))
   allocate(jind(n,n))
   jind = 0

!---- Write dummy jacobian
   do i = 1,n
      jind(i,i) = 1
   end do

   do i=1, nreac
      rr_tmp = rrate(i)
      select case (rr_tmp%group)
      case(1:3,11)
         if ((rr_tmp%reac_type.eq.rrt_sf).or.(rr_tmp%reac_type.eq.rrt_bf)) cycle  ! fission
         do j=1,5
            if (rr_tmp%parts(j).eq.0) exit
            jind(rr_tmp%parts(1),rr_tmp%parts(j)) = 1
         end do
      case(4:7)
         if (rr_tmp%reac_type.eq.rrt_nf) cycle           ! n-induced fission reaction
         do j=1,6
            if (rr_tmp%parts(j).eq.0) exit
            jind(rr_tmp%parts(1),rr_tmp%parts(j)) = 1
            jind(rr_tmp%parts(2),rr_tmp%parts(j)) = 1
         end do
      case(8:9)
         do j=1,6
            if (rr_tmp%parts(j).eq.0) exit
            jind(rr_tmp%parts(1),rr_tmp%parts(j)) = 1
            jind(rr_tmp%parts(2),rr_tmp%parts(j)) = 1
            jind(rr_tmp%parts(3),rr_tmp%parts(j)) = 1
         end do
      case(10)
         do j=1,6
            if (rr_tmp%parts(j).eq.0) exit
            jind(rr_tmp%parts(1),rr_tmp%parts(j)) = 1
            jind(rr_tmp%parts(2),rr_tmp%parts(j)) = 1
            jind(rr_tmp%parts(3),rr_tmp%parts(j)) = 1
            jind(rr_tmp%parts(4),rr_tmp%parts(j)) = 1
         end do
      end select
   end do

   if (fissflag.ne.0) then
      do i=1,nfiss
         do j=1,fissrate(i)%dimens
            jind(fissrate(i)%fissnuc_index,fissrate(i)%fissparts(j)) = 1
            if (fissrate(i)%reac_type.eq. rrt_nf) then               ! neutron-induced fission
               jind(ineu,fissrate(i)%fissparts(j)) = 1
            end if
         end do
      end do
   end if

!---- convert dummy jacobian to compressed sparse column format and
!---- write the index in a for a combination (i,j) into the jacobian
!---- at the corresponding position
   ja_t = 0
   pt_b = -1
   cnt = 0
   do i=1,n
      do j=1,n
         if (jind(j,i).eq.0) cycle
         cnt = cnt + 1
         if (j.eq.i) dia(i)=cnt
         jind(j,i) = cnt
         ja_t(cnt) = j
         if (pt_b(i).lt.0)pt_b(i)=cnt
      end do
      pt_e(i) = cnt + 1
   end do
   pt_b(n+1) = cnt + 1

   allocate(vals(cnt))
   allocate(rows(cnt))
   rows = ja_t(1:cnt)

!---- write the index in the dummy jacobian into cscf_ind of the
!---- corresponding reactionrates
   do i=1, nreac
      rr_tmp = rrate(i)
      cscf_tmp = 0
      select case (rr_tmp%group)
      case(1:3,11)
         if ((rr_tmp%reac_type.eq.rrt_sf).or.(rr_tmp%reac_type.eq.rrt_bf)) cycle      ! fission
         do j=1,5
            if (rr_tmp%parts(j).eq.0) exit
            cscf_tmp(1,j)=jind(rr_tmp%parts(1),rr_tmp%parts(j))
         end do
      case(4:7)
         if (rr_tmp%reac_type.eq.rrt_nf) cycle                    ! n-induced fission
         do j=1,6
            if (rr_tmp%parts(j).eq.0) exit
            cscf_tmp(1,j)=jind(rr_tmp%parts(1),rr_tmp%parts(j))
            cscf_tmp(2,j)=jind(rr_tmp%parts(2),rr_tmp%parts(j))
         end do
      case(8:9)
         do j=1,6
            if (rr_tmp%parts(j).eq.0) exit
            cscf_tmp(1,j)=jind(rr_tmp%parts(1),rr_tmp%parts(j))
            cscf_tmp(2,j)=jind(rr_tmp%parts(2),rr_tmp%parts(j))
            cscf_tmp(3,j)=jind(rr_tmp%parts(3),rr_tmp%parts(j))
         end do
      case(10)
         do j=1,6
            if (rr_tmp%parts(j).eq.0) exit
            cscf_tmp(1,j)=jind(rr_tmp%parts(1),rr_tmp%parts(j))
            cscf_tmp(2,j)=jind(rr_tmp%parts(2),rr_tmp%parts(j))
            cscf_tmp(3,j)=jind(rr_tmp%parts(3),rr_tmp%parts(j))
            cscf_tmp(4,j)=jind(rr_tmp%parts(4),rr_tmp%parts(j))
         end do
      end select

      rrate(i)%cscf_ind = cscf_tmp
   end do

   if (fissflag.ne.0) then
      do i=1,nfiss
        ! Size to correct dimensions, take care that it doesnt get allocated twice
        if (allocated(cscf_tmp_fiss)) deallocate(cscf_tmp_fiss)
        allocate(cscf_tmp_fiss(2,fissrate(i)%dimens))

         do j=1,fissrate(i)%dimens
            if (fissrate(i)%reac_type .eq. rrt_nf) then
               cscf_tmp_fiss(1,j) = jind(ineu,fissrate(i)%fissparts(j))
               cscf_tmp_fiss(2,j) = jind(fissrate(i)%fissnuc_index,fissrate(i)%fissparts(j))
            else
               cscf_tmp_fiss(1,j) = jind(fissrate(i)%fissnuc_index,fissrate(i)%fissparts(j))
            end if
         end do
         allocate(fissrate(i)%cscf_ind(2,fissrate(i)%dimens))
         fissrate(i)%cscf_ind = cscf_tmp_fiss
      end do
   end if

   INFO_EXIT("sparse")

end subroutine sparse


!>
!! The solver core
!!
!! This subroutine calls the sparse PARDISO solver
!! [intel pardiso](https://software.intel.com/en-us/node/470282) to
!! solve equations in the following form:
!! \f[ J\cdot \mathrm{res} = \mathrm{rhs}  \f]
!! with res being the unknown. In addition, this subroutine
!! checks for zero entries in vals and adjusts the sparse
!! format accordingly.
!!
subroutine netsolve(rhs, res)
   use parameter_class, only: solver
   implicit none

   real(r_kind),dimension(net_size),intent(inout)  :: rhs      !< right-hand sides of the system
   real(r_kind),dimension(net_size),intent(inout)  :: res      !< resulting abundances

!!!      definitions for sparse matrix solver pardiso
   integer*8                                       ::pt(64)
!!!     all other variables
   integer                                         :: i, j
   integer                                         :: maxfct, mnum, mtype, phase, nrhs, error, msglvl
   integer                                         :: iparm(64)
   integer                                         :: ia_par(net_size+1)
   integer, dimension(:), allocatable              :: ja_par
   real(r_kind), dimension(:), allocatable         :: a_par
   real(r_kind)                                    :: el(net_size)
   integer                                         :: idum(1)
   real(r_kind)                                    :: ddum(1)
   integer                                         :: cnt
   integer                                         :: nthreads
   character*3                                     :: omp_env
   integer                                         :: nz_cnt

   data nrhs /1/, maxfct /1/, mnum /1/

   INFO_ENTRY("netsolve")

   nz_cnt = count(vals.ne.0.d0)
   allocate(ja_par(nz_cnt))
   allocate(a_par(nz_cnt))
   ! Check if entries are zero and adjust
   ! sparse format
   cnt = 0
   ia_par = -1
   do i=1,net_size
      do j=pt_b(i),pt_e(i)-1
         if(vals(j).eq.0.d0) cycle
         cnt = cnt + 1
         a_par(cnt) = vals(j)
         ja_par(cnt) = rows(j)
         if (ia_par(i).lt.0) ia_par(i) = cnt
      end do
   end do
   ia_par(net_size+1) = cnt+1

   el = rhs

! setup pardiso control parameters and initialize the solvers
! internal adress pointers. this is only necessary for the first
! call of the pardiso solver.
!
   mtype     = 11      ! unsymmetric matrix
   call pardisoinit(pt, mtype, iparm)
! numbers of processors ( value of omp_num_threads )
   call getenv('OMP_NUM_THREADS',omp_env)
   read(omp_env,'(i3)')nthreads
   iparm(3)  =  nthreads
!    iparm(4)  = 61 ???
!    iparm(11) =  1 ???
!    iparm(13) =  1 ???
   msglvl    =  0       ! without statistical information
   iparm(8)  = 10       ! max numbers of iterative refinement steps
   iparm(10) = 13
   phase     = 13  ! analysis, numerical factorization, solve, iterative refinement

   idum = 0
   call pardiso (pt, maxfct, mnum, mtype, phase, net_size, a_par, ia_par, ja_par, &
        idum, nrhs, iparm, msglvl, rhs, el, error)

! termination and release of memory
   phase     = -1           ! release internal memory
   call pardiso (pt, maxfct, mnum, mtype, phase, net_size, ddum, idum, idum,&
        idum, nrhs, iparm, msglvl, ddum, ddum, error)

! check solution for NaNs
   do i=1, net_size
      if (el(i).ne.el(i)) then
         call raise_exception("el(i) is NaN","netsolve",350003)
      endif
   end do

! form output vectors
   select case (solver)
   case(0) ! implicit Euler solver
      do i=1, net_size
         if((i.gt.2).and.(el(i) .lt. 1.d-25)) el(i) = 0.d0
         rhs(i) = el(i) - res(i)
         res(i) = el(i)
      end do
   case(1) ! Gear solver
      res(1:net_size) = el(1:net_size)
   case default
      call raise_exception("Solver flag ("//trim(adjustl(int_to_str(solver)))//") not known."//&
                           NEW_LINE("A")//'Choose either "1" or "2".',"netsolve",350004)
   endselect

! end call pardiso
   deallocate(a_par)
   deallocate(ja_par)
   INFO_EXIT("netsolve")
   return

end subroutine netsolve

end module pardiso_class
