! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd
! Made double precision by John Baumgardner January 2007
! Rename module to "sort_module" and add heapsort subroutine 
!   by Wesley Brewer February 2013

module sort_module

implicit none
public :: QsortC
public :: heapsort
private :: Partition

contains

recursive subroutine QsortC(A)
  real*8, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real*8, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real*8 :: temp
  real*8 :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

subroutine heapsort(a,n)

!   This routine applies the heapsort algorithm to sort array a of
!   length n into ascending numerical order. Array a is replaced on
!   output by its sorted rearrangement.
!   Taken from Press, et al., Numerical Recipes, 1986, p. 231.

implicit none
real*8 a(*), ra
integer i, ir, j, l, n

l  = n/2 + 1
ir = n

!   The index l will be decremented from its initial value down to 1
!   during the hiring (heap creation) phase.  Once it reaches 1, the 
!   index ir will be decremented from its initial value down to 1
!   during the retirement-and-promotion (heap selection) phase.

10 continue

if(l > 1) then         ! Still in hiring phase.
   l  = l - 1
   ra = a(l)
else                   ! In retirement-and-promotion phase.
   ra = a(ir)          ! Clear a space at the end of the array.
   a(ir) = a(1)        ! Retire the top of the heap into it.
   ir = ir - 1         ! Decrease the size of the corporation.
   if(ir == 1) then    ! Completed the last promotion.
      a(1) = ra        ! Identify the least competent worker.
      return
   end if
end if

i = l                  ! Set up to sift element ra to its proper
j = l + l              ! level.

do while(j <= ir)

   if(j < ir) then
      if(a(j) < a(j+1)) j = j + 1  ! Compare to better underling.
   end if

   if(ra < a(j)) then  ! Demote ra.
      a(i) = a(j)
      i = j
      j = j + j
   else                ! This is ra's level. Set j to terminate
      j = ir + 1       ! sift-down.
   end if

end do
   
a(i) = ra              ! Put ra into its slot.

go to 10

end subroutine heapsort

end module sort_module
