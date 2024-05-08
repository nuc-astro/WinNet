!> @file mergesort_class.f90
!!
!! The error file code for this file is ***W29***.
!! Contains the module \ref mergesort_module.

!>
!! Module \ref mergesort_module for merging arrays of rates
!!
!! A mergesort is a recursive algorithm to sort arrays.
!! It works as illustrated in these two figures (Source: en.wikipedia)
!! <div class="row">
!! @image html https://upload.wikimedia.org/wikipedia/commons/e/e6/Merge_sort_algorithm_diagram.svg "Illustration of a mergesort" height=280
!! @image html https://upload.wikimedia.org/wikipedia/commons/c/c5/Merge_sort_animation2.gif "Illustration of a mergesort" height=280 </div>
!!
#include "macros.h"
module mergesort_module
   use error_msg_class, only: raise_exception, int_to_str
  implicit none

  integer, private        :: debug_tabl    !< File ID for debugging tabulated rate class

  !
  ! Public and private fields and methods of the module
  !
  public:: &
      mergesort_init, nurate_ms, rrate_ms, mergesort_finalize, rrate_qs_replace, &
      rrate_sort, bubblesort, reorder, reorder_int, quicksort
  private:: &
      QSort

contains




!>
!! Initialize the mergesort module and open files for debugging.
!!
!! The file _debug_tabln.dat_ is created in case of
!! verbose_level>2.
!!
!! @author Moritz Reichert
!!
!! \b Edited:
!!       - 06.06.17
!! .
subroutine mergesort_init()
   use file_handling_class
   use parameter_class, only: use_tabulated_rates
   implicit none

   ! Debugging
   if (VERBOSE_LEVEL .gt. 2) then

       if (use_tabulated_rates) then
         ! Open debug file and write a header for the tabulated rates
         debug_tabl = open_outfile('debug_tabln.dat')
         write(debug_tabl,'(A)') '# This file contains data to debug the subroutine rrate_qs_replace, contained in mergesort_module.f90'
      end if
   end if

end subroutine mergesort_init



!>
!! Merges two tables of neutrino rates (of type global_class::nurate_type).
!!
!! Parameter r is a flag with the following meaning:
!!
!! <table>
!! <caption id="multi_row">Explanation of the flag "r"</caption>
!! <tr><th> r mode    <th> Behavior
!! <tr><td> 1         <td> Values of y replace corresponding values of x
!! <tr><td> 2         <td> Values of x replace corresponding values of y
!! <tr><td> other     <td> Both values are stored one after the other\
!! </table>
!!
!! @remark At the moment, this function is only called with \f$ r=0 \f$
!!         i.e., case "other".
!!
!! @see rrate_ms, network_init_module::network_init
!!
!! \b Edited:
!!        - 11.01.14
!! .
subroutine nurate_ms(x,xs,y,ys,r)
   use global_class, only: nurate_type, nurate

   implicit none

   integer, intent(in)                        :: xs !< Length of x
   integer, intent(in)                        :: ys !< Length of y
   type(nurate_type),dimension(xs),intent(in) :: x  !< the first nurates table
   type(nurate_type),dimension(ys),intent(in) :: y  !< the second nurates table
   integer, intent(in)                        :: r  !< flag
   !
   type(nurate_type), dimension(:), allocatable   :: tmp
   integer    :: x1,x2,y1,y2,n
   integer    :: ptx, pty, ptz
   integer    :: fin
   integer    :: ptb, pte
   integer    :: i
   integer    :: istat

   n = xs+ys
   ptx = 1
   pty = 1
   ptz = 1
   fin = 0

   allocate(tmp(n),stat=istat)
   if (istat /= 0) call raise_exception('Allocation of "tmp" failed.',&
                                        "nurate_ms",290001)

   do while (fin.eq.0)
      x1 = x(ptx)%parts(1)
      x2 = x(ptx)%parts(2)

      y1 = y(pty)%parts(1)
      y2 = y(pty)%parts(2)

      select case(x1 - y1)
      case(:-1)
         tmp(ptz) = x(ptx)
         ptx = ptx+1
         ptz = ptz+1
      case(1:)
         tmp(ptz) = y(pty)
         pty = pty+1
         ptz = ptz+1
      case(0)
         select case(x2 - y2)
         case(:-1)
            tmp(ptz) = x(ptx)
            ptx = ptx+1
            ptz = ptz+1
         case(1:)
            tmp(ptz) = y(pty)
            pty = pty+1
            ptz = ptz+1
         case(0)
            select case(r)
            case(1)
               tmp(ptz) = y(pty)
               pty = pty+1
               ptx = ptx+1
               ptz = ptz+1
            case(2)
               tmp(ptz) = x(ptx)
               pty = pty+1
               ptx = ptx+1
               ptz = ptz+1
            case default
               tmp(ptz) = x(ptx)
               ptx = ptx+1
               ptz = ptz+1
               tmp(ptz) = y(pty)
               pty = pty+1
               ptz = ptz+1
            end select
         end select
      end select
      if (ptx.gt.xs) fin = fin+1
      if (pty.gt.ys) fin = fin+2
   end do

   if (fin .eq. 1) then
      ptb = pty
      pte = ys
   else if (fin.eq.2) then
      ptb = ptx
      pte = xs
   else
      ptz = ptz -1
      if(allocated(nurate)) deallocate(nurate)
      allocate(nurate(ptz))
      nurate = tmp(1:ptz)
      deallocate(tmp)
      return
   end if

   do i = ptb,pte
      if (fin.eq.1) tmp(ptz) = y(i)
      if (fin.eq.2) tmp(ptz) = x(i)
      ptz = ptz+1
   end do

   ptz = ptz -1
   if(allocated(nurate)) deallocate(nurate)
   allocate(nurate(ptz))
   nurate = tmp(1:ptz)
   deallocate(tmp)

   return
end subroutine nurate_ms


!>
!! Merges two tables of rates (of type global_class::reactionrate_type).
!!
!! Parameter r is a flag with the following meaning:
!!
!! <table>
!! <caption id="multi_row">Evolution mode switches</caption>
!! <tr><th> r mode    <th> Behavior
!! <tr><td> 1         <td> Values of y replace corresponding values of x
!! <tr><td> 2         <td> Values of x replace corresponding values of y
!! <tr><td> other     <td> Both values are stored one after the other\
!! </table>
!!
!! In case of theoretical weak rates that get merged, this function counts
!! the amount of rates that get replaced.
!!
!! @remark The function is used to merge the rate array with theoretical weak rates
!!         (with r=1) and neutrino rates (r=0).
!!
!!
!! @see nurate_ms, network_init_module::network_init
!!
!! \b Edited:
!!      - 17.02.14
!!      - 27.01.21, MR - introduced "rate_out" instead of writing directly to rrate
!!      - 04.08.22, MR - removed hardcoded " ffn" flag and replaced it with "rrs_twr"
!! .
  subroutine rrate_ms(x,xs,y,ys,r,ptz,rate_out)
    use global_class,    only: rrate, reactionrate_type
    use global_class,    only: common_weak_rates, only_theo_weak_rates
    use global_class,    only: rrate_weak_exp
    use error_msg_class, only: raise_exception
    implicit none

    integer, intent(in)                              :: xs !< Length of x
    integer, intent(in)                              :: ys !< Length of y
    type(reactionrate_type),dimension(xs),intent(in) :: x  !< the first nurates table
    type(reactionrate_type),dimension(ys),intent(in) :: y  !< the second nurates table
    integer, intent(in)                              :: r  !< flag
    integer, intent(out)                             :: ptz!< length of new (merged) array
    type(reactionrate_type), dimension(:), &
         allocatable,optional, intent(out) :: rate_out !< returns the merged
                                                      !< reactions in this array, if present
    !
    type(reactionrate_type), dimension(:), allocatable :: tmp
    type(reactionrate_type), dimension(:), allocatable :: tmp2
    integer              :: x1,x2,y1,y2,n
    integer              :: ptx, pty
    integer              :: fin
    integer              :: ptb, pte
    integer              :: alloc_stat !< Allocation status flag
    integer              :: i
    integer              :: rp

    n = xs+ys
    ptx = 1
    pty = 1
    ptz = 1
    fin = 0
    rp  = 0

    allocate(tmp(n))

    !common_weak_rates     = 0
    !only_theo_weak_rates  = 0
    if (y(pty)%reac_src .eq. rrs_twr) then
       if (.not. allocated(tmp2)) allocate(tmp2(ys))
    end if

    do while (fin.eq.0)
       x1 = x(ptx)%parts(1)
       x2 = x(ptx)%parts(2)

       y1 = y(pty)%parts(1)
       y2 = y(pty)%parts(2)

       select case(r)
       case(1)
          if(ptz.gt.1) then
             if(rp.eq.1) then
                if((x1.eq.tmp(ptz-1)%parts(1)).and.(x2.eq.tmp(ptz-1)%parts(2))) then
                   ptx=ptx+1
                   cycle
                end if
             end if
          end if
       case(2)
          if(ptz.gt.1) then
             if(rp.eq.1) then
                if((y1.eq.tmp(ptz-1)%parts(1)).and.(y2.eq.tmp(ptz-1)%parts(2))) then
                   pty=pty+1
                   cycle
                end if
             end if
          end if
       end select
       rp = 0
       select case(x1 - y1)
       case(:-1)
          tmp(ptz) = x(ptx)
          ptx = ptx+1
          ptz = ptz+1
       case(1:)
          tmp(ptz) = y(pty)
          if (y(pty)%reac_src .eq. rrs_twr) then
             only_theo_weak_rates = only_theo_weak_rates + 1
          end if
          pty = pty+1
          ptz = ptz+1
       case(0)
          select case(x2 - y2)
          case(:-1)
             tmp(ptz) = x(ptx)
             ptx = ptx+1
             ptz = ptz+1
          case(1:)
             tmp(ptz) = y(pty)
             if (y(pty)%reac_src .eq. rrs_twr) then
                only_theo_weak_rates = only_theo_weak_rates + 1
             end if
             pty = pty+1
             ptz = ptz+1
          case(0)
             select case(r)
             case(1)
                tmp(ptz) = y(pty)
                pty = pty+1
                ptx = ptx+1
                ptz = ptz+1
                rp = 1
             case(2)
                tmp(ptz) = x(ptx)
                pty = pty+1
                ptx = ptx+1
                ptz = ptz+1
                rp = 1
             case default
                tmp(ptz) = x(ptx)
                ptx = ptx+1
                ptz = ptz+1
                tmp(ptz) = y(pty)
                pty = pty+1
                ptz = ptz+1
             end select
             if (y(pty-1)%reac_src .eq. rrs_twr) then
                common_weak_rates = common_weak_rates + 1
                tmp2(common_weak_rates) = x(ptx-1)
             end if
          end select
       end select
       if (ptx.gt.xs) fin = fin+1
       if (pty.gt.ys) fin = fin+2
    end do

    if (fin .eq. 1) then
       ptb = pty
       pte = ys
    else if (fin.eq.2) then
       ptb = ptx
       pte = xs
    else
       ptz = ptz -1
       if (present(rate_out))then
         if (allocated(rate_out)) deallocate(rate_out)
         allocate(rate_out(ptz),stat=alloc_stat)
         rate_out(1:ptz) = tmp(1:ptz)
       else
         rrate(1:ptz) = tmp(1:ptz)
       end if
       deallocate(tmp)
       return
    end if

    do i = ptb,pte
       if (fin.eq.1) then
          tmp(ptz) = y(i)
          if (y(i)%reac_src .eq. rrs_twr) then
             only_theo_weak_rates = only_theo_weak_rates + 1
          end if
       end if
       if (fin.eq.2) tmp(ptz) = x(i)
       ptz = ptz+1
    end do

    ptz = ptz -1
    if (present(rate_out)) then
      if (allocated(rate_out)) deallocate(rate_out)
      allocate(rate_out(ptz),stat=alloc_stat)
      if ( alloc_stat /= 0) call raise_exception('Allocation of "rrate_array" failed.',&
                                                 "merge_neutrino_rates",290001)
      rate_out(1:ptz) = tmp(1:ptz)
    else
      rrate(1:ptz) = tmp(1:ptz) !< store directly into rrate and hope for the best
    end if

    deallocate(tmp)

    if (y(1)%reac_src .eq. rrs_twr) then
       if (.not. allocated(rrate_weak_exp)) allocate(rrate_weak_exp(common_weak_rates))
       rrate_weak_exp(1:common_weak_rates) = tmp2(1:common_weak_rates)
       deallocate(tmp2)
    end if

    return
  end subroutine rrate_ms



!> Wrapper around the quicksort subroutine
!!
!! This subroutine is able to replace rates that appear with
!! multiple entries in the REACLIB (e.g., because it has a resonance)
!! with a tabulated rate. In debug mode, the file debug_tabln.dat
!! is written to debug this subroutine.
!!
!! @see qsort
!!
!! @author D. Martin
!!
!! \b Edited:
!!     - 15.04.15
!!     - 07.06.17
!!     - 04.08.22, MR: removed hardcoded "tabl" flag and replaced it with "rrs_tabl"
!!     - 05.08.22, MR: fixed bug and gave a new array as output (out_array)
!! .
  subroutine rrate_qs_replace(x,xs,y,ys,r,out_array,ptz)
    use global_class, only: reactionrate_type,net_names

    implicit none

    integer, intent(in)                                :: r,xs,ys !< replacement rule, lenghts of arrays
    type(reactionrate_type), dimension(xs), intent(in) :: x       !< array 1 (rrate)
    type(reactionrate_type), dimension(ys), intent(in) :: y       !< array 2 (tabl)
    type(reactionrate_type), dimension(:),allocatable, intent(out) :: out_array  !< new merged array
    integer, intent(out)                               :: ptz     !y final length of merged array
    !
    type(reactionrate_type), dimension(:), allocatable :: tmp
    integer                                            :: n, i, j, k
    integer                                            :: nrDuplicates
    logical                                            :: duplicate
    type(reactionrate_type)                            :: tmp_find_duplicates
    integer                                            :: count_duplicates
    integer                                            :: tabl_position
    integer                                            :: new_rate_length
    character(5),dimension(6)                          :: part_names       ! For debugging

    duplicate = .false.


    ! prepare array to be sorted
    n = xs+ys
    allocate(tmp(n))
    tmp(1:xs)   = x
    tmp(xs+1:n) = y

    new_rate_length = n

    ! sort array
    call QSort(tmp,n)
    ! remove duplicates and create array without duplicates, if r!=0

    j = 1
    nrDuplicates = 0
    if(r .ne. 0) then
      do i=1,n-1
        ! Prepare loops to identify the number of the duplicates and the position of the tabulated rate
        tmp_find_duplicates = tmp(i)
        count_duplicates = 1
        if (tmp_find_duplicates%reac_src .eq. rrs_tabl) then
           tabl_position = i
        else
           tabl_position = 0
        end if

        do k =i, new_rate_length-1
           if (all(tmp_find_duplicates%parts .eq. tmp(k+1)%parts) .and. &
              (tmp_find_duplicates%group .eq. tmp(k+1)%group)) then
             ! Count the amount of duplicates
              count_duplicates= count_duplicates + 1
              ! Check the position of the tabulated rate
              if (tmp(k+1)%reac_src .eq. rrs_tabl) then
                 tabl_position = k+1
              end if
           else
             exit
           end if
        end do


        ! If there was a tabulated rate included, replace the first rate with that one
        if (tabl_position .ne. 0) then
           ! Write debug file
           if (VERBOSE_LEVEL .gt. 2) then
              ! Convert parts to corresponding names
              do k=1, 6
                 if (tmp_find_duplicates%parts(k) .ne. 0) then
                    part_names(k) = net_names(tmp_find_duplicates%parts(k))
                 else
                    part_names(k) = '     '
                 end if
              end do
              ! Write an output
              write(debug_tabl,"('Replacing rate ',6a,' found ',i10,' reactions in Reaclib')") part_names(:),count_duplicates-1
           end if

           ! Replace the first rate with the tabulated one
           tmp(i) = tmp(tabl_position)
           do k=i,new_rate_length - count_duplicates
             tmp(k+1)=tmp(k+count_duplicates)
           end do
           new_rate_length = (new_rate_length - count_duplicates) + 1
        end if

        ! Exit condition
        if (i .ge. new_rate_length) then
           exit
        end if
      enddo


      if (VERBOSE_LEVEL .gt. 2) then
         ! Write an output
         write(debug_tabl,*)
         write(debug_tabl,'(A,i10,A)') 'Merged in ', n-xs, ' rates.'
      end if

    endif
    if (VERBOSE_LEVEL .gt. 2)  print *, "duplicates: ", new_rate_length

    ptz = new_rate_length

    allocate(out_array(new_rate_length))
    ! TODO throw exception
    ! assign rates to rrate
    out_array(1:new_rate_length) = tmp(1:new_rate_length)
    deallocate(tmp)

    return
  end subroutine rrate_qs_replace



!> Modified quicksort algorithm
!!
!! A quicksort is a recursive sorting algorithm. An illustration of it is shown
!! in the following (Source en.wikipedia):
!! @image html https://upload.wikimedia.org/wikipedia/commons/6/6a/Sorting_quicksort_anim.gif "Illustration of the quicksort" width=300
!! In this gif, the red indicated number is the so called pivot-element.
!!
!! @see [Rosettacode](http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran),
!!      rateSortKey, rrate_qs_replace
!!
!! \b Edited:
!!       -  15.04.15
!!       -  25.4.19, M.R.
!! .
  recursive subroutine QSort(A,nA)
    use global_class, only: reactionrate_type

    implicit none

    ! DUMMY ARGUMENTS
    type (reactionrate_type), dimension(nA), intent(in out) :: A !< Reaction rate array
    integer, intent(in) :: nA  !< Length of the investigated reaction rate array.
                               !! The recursion anchor of the routine is given
                               !! for nA \f$\le =1 \f$.

    ! LOCAL VARIABLES
    character(31)                          :: rateSortKey
    integer                                :: left, right
    real                                   :: random
    character(31)                          :: pivot
    type (reactionrate_type)               :: temp
    integer                                :: marker, randI, clock
    integer                                :: n
    integer,dimension(:) ,allocatable      :: variable_put
    integer                                :: istat


    if (nA > 1) then

        call system_clock(count=clock)

        ! Get the size in order to call put properly and system independent
        call random_seed(size=n)
        allocate(variable_put(n),stat=istat)
        if (istat /= 0) call raise_exception('Allocation failed.',"QSort",290001)
        variable_put(:) = clock

        call random_seed(put = variable_put)
        call random_number(random)
        randI = int(random*real(nA-1))+1
        pivot = rateSortKey(A(randI)%group, A(randI)%parts)   ! random pivot (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1

        do while (left < right)
            right = right - 1
            do while (LGT(rateSortKey(A(right)%group, A(right)%parts), pivot))
                right = right - 1
            end do
            left = left + 1
            do while (LLT(rateSortKey(A(left)%group, A(left)%parts), pivot))
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
            end if
        end do

        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if

        call QSort(A(:marker-1),marker-1)
        call QSort(A(marker:),nA-marker+1)

    end if

  end subroutine QSort




  !> Close files and finalize this module
  !!
  !! This subroutine closes the debug file "debug_tabln.dat"
  !!
  !! @author: Moritz Reichert
  !!
  !! \b Edited:
  !!     - 07.06.17
  !! .
  subroutine mergesort_finalize()
    use file_handling_class
    use parameter_class, only: use_tabulated_rates
    implicit none

    ! Debugging
    if (VERBOSE_LEVEL .gt. 2) then
       if (use_tabulated_rates) then
          ! Close the debug file again
          call close_io_file(debug_tabl,'debug_tabln.dat')
      end if
    end if

  end subroutine mergesort_finalize

  !> Sorts chapter one of the rate array
  !!
  !! The sorting is based on the second nucleus in the rate array.
  !!
  !! ### Example
  !! <pre>
  !! Rate 1: cl34 (Z=17) ->  s34 (Z=16)
  !! Rate 2: cl36 (Z=17) -> ar36 (Z=18)
  !! Rate 3: cl36 (Z=17) ->  s36 (Z=16)
  !! Rate 4: cl38 (Z=17) -> ar38 (Z=18)
  !! </pre>
  !! will switch the order of the rates to:
  !! <pre>
  !! Rate 1: cl34 (Z=17) ->  s34 (Z=16)
  !! Rate 3: cl36 (Z=17) ->  s36 (Z=16)
  !! Rate 2: cl36 (Z=17) -> ar36 (Z=18)
  !! Rate 4: cl38 (Z=17) -> ar38 (Z=18)
  !! </pre>
  !! In case the sunet is ordered by proton number.
  !!
  !! @todo Expand and change this subroutine
  !!
  !! @warning This subroutine assumes a certain ordering of the
  !!          sunet!
  subroutine rrate_sort(num,rate_array)
    use global_class, only: reactionrate_type

    implicit none

    integer,intent(in)                                    :: num        !< Number of rates in chapter 1
    type(reactionrate_type),dimension(num),intent(inout)  :: rate_array !< Rate array that will get sorted
    type(reactionrate_type),dimension(1)                  :: temp       !< Temporary reaction rate
    integer     :: i

    INFO_ENTRY("rrate_sort")
    i=1

    do i=1,num-1
       if (rate_array(i+1)%group.ne.1)exit
       if ((rate_array(i)%parts(1).eq.rate_array(i+1)%parts(1)) .and. &
            (rate_array(i)%parts(2).gt.rate_array(i+1)%parts(2))) then
          if (VERBOSE_LEVEL .ge. 2) then
             print '(a10,4i6)','switch ' ,rate_array(i)%parts(1:2), rate_array(i+1)%parts(1:2)
          end if
          temp(1) = rate_array(i+1)
          rate_array(i+1) = rate_array(i)
          rate_array(i) = temp(1)
       end if
    end do

    INFO_EXIT("rrate_sort")

    return
  end subroutine rrate_sort



  !> Bubblesort of array.
  !!
  !! This subroutine performs a bubblesort of an array. See the following
  !! illustration (from wikipedia):
  !! @image html https://upload.wikimedia.org/wikipedia/commons/c/c8/Bubble-sort-example-300px.gif "Illustration of a bubblesort" width=300
  !! The subroutine is able to run in 2 modes,
  !! - mode 0 : Ascending order
  !! - mode 1 : Descending order
  !! .
  !! Furthermore, a second array can be given which will return the indices of
  !! the changes within the first array.
  !!
  !! @see reorder, tw_rate_module::readweak_logft
  !!
  !! @author M. Reichert
  subroutine bubblesort(mode,length,arr1,arr2)
    use error_msg_class, only: int_to_str,raise_exception
    implicit none
    real(r_kind),dimension(:),intent(inout)          :: arr1       !< Array to get sorted
    integer,dimension(:),intent(out),optional        :: arr2       !< Indices of the sort
    integer,intent(in)                               :: length     !< Length of the arrays
    integer,intent(in)                               :: mode       !<0 descending; 1 ascending
    ! Internal helper variables
    logical      :: unsorted
    integer      :: i
    real(r_kind) :: h

    !Create indices in array2
    if (present(arr2)) then
      do i=1, length
        arr2(i) = i
      end do
    end if

    unsorted = .True.

    do while(unsorted)
      unsorted = .False.

      do i = 1 , length-1
        if (mode .eq. 1) then
          ! descending order
          if (arr1(i)>arr1(i+1)) then
            h         = arr1(i)
            arr1(i)   = arr1(i+1)
            arr1(i+1) = h
            ! Also sort arr2 if present
            if (present(arr2)) then
              h         = arr2(i)
              arr2(i)   = arr2(i+1)
              arr2(i+1) = int(h)
            end if
            unsorted=.True.
          end if
        elseif (mode .eq. 0) then
          ! ascending order
          if (arr1(i)<arr1(i+1)) then
            h         = arr1(i)
            arr1(i)   = arr1(i+1)
            arr1(i+1) = h
            ! Also sort arr2 if present
            if (present(arr2)) then
              h         = arr2(i)
              arr2(i)   = arr2(i+1)
              arr2(i+1) = int(h)
            end if
            unsorted=.True.
          end if
        else
          call raise_exception("Unknown sort mode "//int_to_str(mode),&
                               "bubblesort",290003)
        end if
      end do
    end do

  end subroutine bubblesort



  !> Quicksort of array.
  !!
  !! This subroutine performs a quicksort of an array.
  !! A quicksort is a recursive sorting algorithm. An illustration of it is shown
  !! in the following (Source en.wikipedia):
  !! @image html https://upload.wikimedia.org/wikipedia/commons/6/6a/Sorting_quicksort_anim.gif "Illustration of the quicksort" width=300
  !! In this gif, the red indicated number is the so called pivot-element.
  !! In contrast to \ref QSort, this subroutine is able to sort an array that is not a derived type and comes in form of
  !! an array of doubles.
  !! The subroutine is able to run in 2 modes,
  !! - mode 0 : Ascending order
  !! - mode 1 : Descending order
  !! .
  !! Furthermore, a second array can be given which will return the indices of
  !! the changes within the first array.
  !!
  !! @see [Rosettacode](http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran),
  !!
  !! @see reorder
  !!
  !! @author M. Reichert
  !! @date 08.05.24
  recursive subroutine quicksort(mode, length, arr1, arr2, is_recursion)
    use error_msg_class, only: int_to_str, raise_exception
    implicit none
    real(r_kind), dimension(:), intent(inout) :: arr1       !< Array to get sorted
    integer, dimension(:), intent(out), optional :: arr2    !< Indices of the sort
    integer, intent(in) :: length                           !< Length of the arrays
    integer, intent(in) :: mode                             !< 0 descending; 1 ascending
    logical, intent(in), optional :: is_recursion           !< Flag to indicate if the subroutine is called recursively
    ! Internal helper variables
    real(r_kind) :: pivot_value
    integer :: pivot_index, i, j
    logical :: recursive_call
    real(r_kind) :: h

    !Create indices in array2
    if (present(is_recursion)) then
        recursive_call = is_recursion
    else
        recursive_call = .False.
    end if

    if (.not. recursive_call) then
        if (present(arr2)) then
            arr2 = [(i, i=1,length)]
        end if
    end if

    ! Check for base case
    if (length <= 1) return

    ! Choose pivot element
    pivot_index = length / 2
    pivot_value = arr1(pivot_index)

    ! Partition step
    i = 1
    j = length

    do while (i <= j)
      if (mode == 1) then
        do while (arr1(i) < pivot_value)
          i = i + 1
        end do
        do while (arr1(j) > pivot_value)
          j = j - 1
        end do
      else
        do while (arr1(i) > pivot_value)
          i = i + 1
        end do
        do while (arr1(j) < pivot_value)
          j = j - 1
        end do
      end if

      if (i <= j) then
        ! Swap elements
        h = arr1(i)
        arr1(i) = arr1(j)
        arr1(j) = h

        if (present(arr2)) then
          h = arr2(i)
          arr2(i) = arr2(j)
          arr2(j) = h
        end if

        i = i + 1
        j = j - 1
      end if
    end do

    ! Recursive calls
    if (present(arr2)) then
        call quicksort(mode, j, arr1, arr2, .True.)
        call quicksort(mode, length - i + 1, arr1(i:), arr2(i:),.True.)
    else
        call quicksort(mode, j, arr1, is_recursion=.True.)
        call quicksort(mode, length - i + 1, arr1(i:), is_recursion=.True.)
    end if

end subroutine quicksort




  !> Takes an array and an array with indices to reorder the first array
  !!
  !! @see bubblesort, tw_rate_module::readweak_logft
  !!
  !! @author M. Reichert
  subroutine reorder(arr1,arr2,length)
    use error_msg_class, only: int_to_str,raise_exception
    implicit none
    real(r_kind),dimension(:),intent(inout) :: arr1   !< Array to reorder
    integer,dimension(:),intent(in)         :: arr2   !< Array containing the indices
    integer,intent(in)                      :: length !< Length of the arrays

    real(r_kind),dimension(length)          :: helper !< Helper variable
    integer :: i

    do i=1,length
      ! Check for weird things in index array
      if ((arr2(i) .lt. 1) .or. (arr2(i) .gt. length)) then
        call raise_exception("Inconsistent index array! Got index "//int_to_str(arr2(i)),&
                             "reorder",290004)
      end if

      helper(i) = arr1(arr2(i))
    end do

    arr1(:) = helper(:)
  end subroutine reorder


  !> Takes an array and an array with indices to reorder the first array
  !!
  !! Same as \ref reorder, but for integer arrays.
  !!
  !! @see bubblesort, tw_rate_module::readweak_logft
  !!
  !! @author M. Reichert
  subroutine reorder_int(arr1,arr2,length)
    use error_msg_class, only: int_to_str,raise_exception
    implicit none
    integer,dimension(:),intent(inout)      :: arr1   !< Array to reorder
    integer,dimension(:),intent(in)         :: arr2   !< Array containing the indices
    integer,intent(in)                      :: length !< Length of the arrays

    integer,dimension(length)               :: helper !< Helper variable
    integer :: i

    do i=1,length
      ! Check for weird things in index array
      if ((arr2(i) .lt. 1) .or. (arr2(i) .gt. length)) then
        call raise_exception("Inconsistent index array! Got index "//int_to_str(arr2(i)),&
                             "reorder",290004)
      end if

      helper(i) = arr1(arr2(i))
    end do

    arr1(:) = helper(:)
  end subroutine reorder_int


end module mergesort_module





!>
!! Returns a key for sorting a rate array by concatenating isotope indices.
!!
!! Rule: order of elements depends on chapter. \n
!! The here produced string is compared for
!! being lexically greater or smaller in \ref qsort.
!!
!! @returns String with that determines the ordering of an isotope
!!
!! @see qsort
!!
!! \b Edited:
!!      - 15.04.15
!!      - 28.07.22, MR: Introduced new chapters
!! .
character(31) function rateSortKey(group,parts)
  !
  implicit none
  integer              :: group
  integer,dimension(6) :: parts

  select case (group)
    case (1:3,11)
       write (rateSortKey, "(I1,6I5)") group, parts
    case (4:7)
       write (rateSortKey, "(I1,6I5)") group, parts(2), parts(1), parts(3:)
    case (8:9)
       write (rateSortKey, "(I1,6I5)") group, parts(3), parts(2), parts(1), parts(4:)
    case (10)
       write (rateSortKey, "(I1,6I5)") group, parts(4), parts(3), parts(2), parts(1), parts(5:)
  end select

end function rateSortKey
