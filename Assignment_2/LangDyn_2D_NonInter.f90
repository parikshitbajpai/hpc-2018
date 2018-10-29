module globals
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

module BC
  use globals
  use Langevin
  implicit none
contains
  subroutine implement_BC

  end subroutine implement_BC
end module BC

program main
use globals
use Langevin
use BC
implicit none
integer :: i
double precision :: t,t_max,ran1,ran2,m1,m2
double precision :: wtime,begin,end

! Open files
open(11,file='trajectories')
open(12,file='means')

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n))

t=0d0
t_max=100d0

call set_parameters
call initialize_particles

call cpu_time(begin)

do while(t.lt.t_max)
   do i=1,n
      vhx(i)=vx(i)+0.5*ax(i)*dt
      vhy(i)=vy(i)+0.5*ay(i)*dt
      x(i)=x(i)+vhx(i)*dt
      y(i)=y(i)+vhy(i)*dt

!     call impose_bc
      
      ax(i)=0d0                   ! Add forces here if any
      ay(i)=0d0                   ! Add forces here if any

      call random_number(ran1)
      ran1=ran1-0.5d0
      call random_number(ran2)
      ran2=ran2-0.5d0
      ax(i)=ax(i)-pref1*vhx(i)+pref2*ran1
      ay(i)=ay(i)-pref1*vhy(i)+pref2*ran2
      
      vx(i)=vhx(i)+0.5*ax(i)*dt
      vy(i)=vhy(i)+0.5*ay(i)*dt
   end do
   t=t+dt
   do i=1,n
      write(11,*) x(i),y(i)
   end do
   write(11,*) ''
   write(12,*) t,sqrt(sum(x**2+y**2)/real(n,8))
end do

call cpu_time(end)
print *,'Wtime=',end-begin

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay)
! Close files
close(11)
close(12)

end program main
