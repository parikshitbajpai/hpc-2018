module globals
! Global variables
implicit none
integer :: n=1000         ! number of particles
double precision :: L=100d0
double precision, parameter :: pi=2q0*asin(1q0) ! numerical constant
end module globals

module Langevin
! Initialization and update rule for Langevin particles
use globals
implicit none
double precision :: dt,kT,g,m                   ! time step size and physical parameters
double precision :: pref1,pref2                 ! auxiliary parameters
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy,x0,y0,vx0,vy0
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
   ! In preparation for the introduction of interactions, spread the particles out
   ! but remember their initial position to compute the rms displacament.
   call random_number(ran1)
   call random_number(ran2)
   x(i)=L*(ran1-0.5d0)
   x0(i)=x(i)
   y(i)=L*(ran2-0.5d0)
   y0(i)=y(i)
   ax(i)=0d0
   ay(i)=0d0
   call random_number(ran1)
   call random_number(ran2)
   gr1=sqrt(kT/(m))*sqrt(-2*log(ran1))*cos(2*pi*ran2) ! Box-Mueller transform
   gr2=sqrt(kT/(m))*sqrt(-2*log(ran1))*sin(2*pi*ran2)
   vx(i)=gr1
   vx0(i)=vx(i)
   vy(i)=gr2
   vy0(i)=vy(i)
end do

end subroutine initialize_particles
end module Langevin

module BC
  ! Subroutines related to the boundary conditions
  use globals
  use Langevin
  implicit none
contains
  subroutine impose_BC(i)
    integer :: i
    ! See if particle i is outside the box [-L/2,L/2] X [-L/2,L/2]
    ! if it is, relfect its position in the boundary and flip its velocity
    if(x(i).lt.-L/2d0) then
       x(i)=-L-x(i)
       vhx(i)=-vhx(i)
    end if
    if(x(i).gt.L/2d0) then
       x(i)=L-x(i)
       vhx(i)=-vhx(i)
    end if
    if(y(i).lt.-L/2d0) then
       y(i)=-L-y(i)
       vhy(i)=-vhy(i)
    end if
    if(y(i).gt.L/2d0) then
       y(i)=L-y(i)
       vhy(i)=-vhy(i)
    end if

  end subroutine impose_BC
end module BC

program main
use globals
use Langevin
use BC
implicit none
integer :: i, skip
double precision :: t,t_max,ran1,ran2,m1,m2
double precision :: wtime,begin,end

! Open files
open(11,file='trajectories')
open(12,file='means')

! Open additional file to save the velocity correlation function
open(13,file='velocity_correlation')

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),x0(n),y0(n),vx0(n),vy0(n))

t=0d0
t_max=100d0
skip=50

call set_parameters
call initialize_particles

call cpu_time(begin)

do while(t.lt.t_max)
   do i=1,n
      vhx(i)=vx(i)+0.5d0*ax(i)*dt
      vhy(i)=vy(i)+0.5d0*ay(i)*dt
      x(i)=x(i)+vhx(i)*dt
      y(i)=y(i)+vhy(i)*dt

      call impose_bc(i)
      
      ax(i)=0d0                   ! Add forces here if any
      ay(i)=0d0                   ! Add forces here if any

      call random_number(ran1)
      ran1=ran1-0.5d0
      call random_number(ran2)
      ran2=ran2-0.5d0
      ax(i)=ax(i)-pref1*vhx(i)+pref2*ran1
      ay(i)=ay(i)-pref1*vhy(i)+pref2*ran2
      
      vx(i)=vhx(i)+0.5d0*ax(i)*dt
      vy(i)=vhy(i)+0.5d0*ay(i)*dt
   end do
   t=t+dt
   do i=1,n,skip
      write(11,*) x(i),y(i)
   end do
   write(11,*) ''
   write(12,*) t,sqrt(sum((x-x0)**2+(y-y0)**2)/real(n,8))
   write(13,*) t,(sum(vx*vx0 + vy*vy0)/real(n,8))/(sum(vx*vx + vy*vy)/real(n,8))
end do

call cpu_time(end)
print *,'Wtime=',end-begin

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,x0,y0)
! Close files
close(11)
close(12)

end program main
