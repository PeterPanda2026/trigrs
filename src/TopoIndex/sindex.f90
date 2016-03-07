subroutine sindex(n, ra, ndx)
! Index sort: returns ndx such that ra(ndx(1)) <= ra(ndx(2)) <= ... <= ra(ndx(n)).
! The array ra is not modified.
!
! Method: heapsort on the index array with bottom-up heap
! construction.  A max-heap keyed on ra is
! built over ndx, then the largest remaining entry is repeatedly moved
! from the top of the heap to the end of the unsorted part.
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: ra(n)
  integer, intent(out) :: ndx(n)
  integer :: k, last, top

  do k = 1, n
    ndx(k) = k
  end do
  if (n < 2) return

  ! Build the heap: sift down every node that has a child, last parent first.
  do k = n / 2, 1, -1
    call sift_down(k, n)
  end do

  ! Repeatedly move the largest entry behind the shrinking heap.
  do last = n, 2, -1
    top = ndx(1)
    ndx(1) = ndx(last)
    ndx(last) = top
    call sift_down(1, last - 1)
  end do

contains

  subroutine sift_down(start, heap_size)
    ! Move the entry at position start down until neither child is larger.
    integer, intent(in) :: start, heap_size
    integer :: parent, child, moving
    real :: key

    moving = ndx(start)
    key = ra(moving)
    parent = start
    do
      child = 2 * parent
      if (child > heap_size) exit
      ! pick the larger child; the left one wins ties
      if (child < heap_size) then
        if (ra(ndx(child + 1)) > ra(ndx(child))) child = child + 1
      end if
      if (ra(ndx(child)) <= key) exit
      ndx(parent) = ndx(child)
      parent = child
    end do
    ndx(parent) = moving
  end subroutine sift_down

end subroutine sindex
