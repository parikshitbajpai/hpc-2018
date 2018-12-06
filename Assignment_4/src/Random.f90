include 'mkl_vsl.f90'                     ! For quasi-random number generation with the MKL VSL library
module globals
USE MKL_VSL_TYPE                          ! Routines for quasi-random 
USE MKL_VSL                               ! number generation
implicit none
include 'mpif.h'
integer :: n, Ni, INFO,i, j, flag, job_count, ith, bin, n_bin
integer, allocatable, dimension(:) :: Tally
double precision :: EV, bin_width
integer :: num_procs,proc_num,DEFAULT_TAG=1
end module globals

module auxiliary
use globals
implicit none
TYPE (VSL_STREAM_STATE) :: stream
contains
subroutine open_random_stream
! Setup of stream of pseudo-random numbers using MKL library.
integer :: un,seed,istat,ierr

! Read seed integer from /dev/urandom, a file filled with pseudo-random numbers generated from the state of the computer
open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
   read(un) seed
close(un)
!print '(a11,i14)','Random seed=',seed

! Use a Mersenne Twister pseudorandom number generator; the PRN sequence is identified by the variable "stream"
ierr=vslnewstream(stream,VSL_BRNG_MT19937,seed)
if(ierr.ne.0) then
   print *,'Trouble with the MKL RNG, returns flag ',ierr
end if

end subroutine open_random_stream

subroutine close_random_stream
! Close pseudo-random number stream
integer :: rnd_ierr

rnd_ierr=vsldeletestream( stream )

end subroutine close_random_stream
end module auxiliary


