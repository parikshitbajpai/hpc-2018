!------------------------Matrix Vector Multiplication--------------------!
!   This code executes matrix-vector multiplication using double loops.  !


module globals
! Globally used stuff
implicit none
integer :: n
end module globals

module auxiliary
! Useful stuff
use globals
implicit none
double precision, allocatable, dimension(:,:) :: A
double precision, allocatable, dimension(:) :: v,w
contains
subroutine make_random_matrix

call random_number(A)

end subroutine make_random_matrix

subroutine make_random_vector

call random_number(v)

end subroutine make_random_vector
end module auxiliary

program main
! Naive matrix-vector product computation
use globals
use auxiliary
implicit none
integer :: i,j,p
double precision :: wtime,begin,end

! Loop over matrix sizes
do p=9,14
   n=2**p

! Allocate and de-allocate inside the loop
   allocate(A(n,n),v(n),w(n))
   
   call make_random_matrix
   call make_random_vector

! GNU Fortran implementation of the CPU clock
   call cpu_time(begin)

   do i=1,n
      w(i)=0d0
      do j=1,n
         w(i)=w(i)+A(i,j)*v(j)
      end do
   end do

   call cpu_time(end)
   wtime=end-begin

   print *,n,wtime

   deallocate(A,v,w)
end do
end program main
