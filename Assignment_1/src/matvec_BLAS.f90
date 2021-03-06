!------------------------Matrix Vector Multiplication--------------------!
!       This code executes matrix-vector multiplication using BLAS.      !

!    Note: The matrix size in this program has been fixed to 16384*16384   !

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

! Set matrix size to 2^14 = 16384. To change matrix size edit p
   p=14
   n=2**p

! Allocate and de-allocate inside the loop
   allocate(A(n,n),v(n),w(n))
   
   call make_random_matrix
   call make_random_vector

! GNU Fortran implementation of the CPU clock
   call cpu_time(begin)

   call dgemv("N", n, n, 1.0D0, a, n, v, 1, 0.0D0, w, 1) ! Call to DGEMV routine in BLAS
   !now =yAu
   call cpu_time(end)
   wtime=end-begin

   print *,n,wtime

   deallocate(A,v,w)

end program main
