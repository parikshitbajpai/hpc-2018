!------------------------Matrix Vector Multiplication--------------------!
!    This code executes matrix-vector multiplication using double loops  !

!   Note: The matrix size in this program has been fixed to 16384*16384  !


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

! Fix matrix size to 2^14 = 16384. Change p to change matrix size 
p=14
n=2**p

! Allocate the variables.
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
end program main
