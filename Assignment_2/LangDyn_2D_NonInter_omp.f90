module globals
use omp_lib
! Global variables
implicit none
integer :: n=1000         ! number of particles
double precision, parameter :: pi=2q0*asin(1q0) ! numerical constant
end module globals

module Langevin
! Initialization and update rule for Langevin particles
use globals
implicit none
double precision :: dt,kT,g,m                   ! time step size and physical parameters
double precision :: pref1,pref2                 ! auxiliary parameters
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy
contains
subroutine set_parameters

! Set time step and physical parameters
dt=0.01d0
kT=1d0
g=1d0
m=1d0

! Set auxiliary parameters
pref1=g
pref2=sqrt(24d0*kT*g/dt)

end subroutine set_parameters
subroutine initialize_particles
integer :: i
double precision :: ran1,ran2,gr1,gr2
! Give particles initial position and velocity
do i=1,n
   x(i)=0d0
   y(i)=0d0
   ax(i)=0d0
   ay(i)=0d0
   call random_number(ran1)
   call random_number(ran2)
   gr1=sqrt(kT/(m))*sqrt(-2*log(ran1))*cos(2*pi*ran2) ! Box-Mueller transform
   gr2=sqrt(kT/(m))*sqrt(-2*log(ran1))*sin(2*pi*ran2)
   vx(i)=gr1
   vy(i)=gr2
end do

end subroutine initialize_particles
end module Langevin

program main
use globals
use Langevin
implicit none
integer :: i,skip
double precision :: t,t_max,m1,m2,dt2
double precision :: wtime,t0,t1
double precision, allocatable, dimension(:) :: ran1, ran2,x_temp,y_temp
! Open files
open(11,file='trajectories')
open(12,file='means')
call omp_set_num_threads(2)

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),x_temp(n),y_temp(n))
allocate(ran1(n),ran2(n))
t=0d0
t_max=100d0

call set_parameters
call initialize_particles
dt2=dt/2d0
t0=omp_get_wtime()
! skip=50

!$omp parallel
do while(t.lt.t_max)
! begin omp sections
   x_temp=x
   y_temp=y
   !$omp sections firstprivate(t),lastprivate(t)
   !$omp section
   ! write to disk
   do i=1,n!,skip
      write(11,*) x_temp(i),y_temp(i)
   end do
   write(11,*) ''
   write(12,*) t,sqrt(sum(x_temp**2+y_temp**2)/real(n,8))
   !$omp section
   ! get random numbers
   call random_number(ran1)
   call random_number(ran2)
   ran1=pref2*(ran1-0.5d0)
   ran2=pref2*(ran2-0.5d0)
   vhx=vx+ax*dt2
   vhy=vy+ay*dt2
   x=x+vhx*dt
   y=y+vhy*dt
   ax=0d0                   ! Add forces here if any
   ay=0d0                   ! Add forces here if any
   ax=ax-pref1*vhx+ran1
   ay=ay-pref1*vhy+ran2
   vx=vhx+ax*dt2
   vy=vhy+ay*dt2
   t=t+dt
   !$omp end sections
end do
!$omp end parallel
t1=omp_get_wtime()

print *,'Wtime=',t1-t0

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,ran1,ran2,x_temp,y_temp)
! Close files
close(11)
close(12)

end program main
